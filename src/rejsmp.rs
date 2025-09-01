//! Weighted MinHash via Efficient Rejection Sampling (ERS).
//!
//! Implements:
//!   - ERS (Li & Li 2021): K independent fixed-length sequences
//!     r_{j,1..L} per hash position j; take the first green if any, otherwise mark
//!     empty; then densify empties by rotating to a non-empty bucket with a
//!     per-j random offset (data-independent).
//!
//! Inputs: sparse weighted vector `&[(u64, f64)]` where id ∈ [0, D) and weight ≥ 0.
//! Randomness: purely via Tab32/Tab64 tabulation hashing (no stateful RNG required).
//!
//! IMPORTANT: Caps `m_i` are **real-valued** (`f64`) and should be set to the
//! *tight* per-dimension maxima across the dataset: `m_i = max_s x_i(s)`.
//! Using tight caps reduces total M = sum_i m_i, increases acceptance probability,
//! and lets you use much smaller L in ERS.

use tab_hash::{Tab32Simple, Tab64Simple};
use rand_core::RngCore;
use crate::hash_utils::{tab32_from_rng, tab64_from_rng, to_unit};
use crate::rng_utils::MtRng;

/// A single (id, rank) pair compatible with your DartMinHash plumbing.
pub type Dart = (u64, f64);

/// Continuous (real-valued) line partition for the red–green test.
/// Stores cumulative caps so we can map r∈[0,M) → component i by binary search.
#[derive(Clone)]
pub struct RedGreenIndex {
    /// cumulative sums: cum[0] = 0, cum[i+1] = cum[i] + m_i  (length = D+1)
    cum: Vec<f64>,
    d: usize,
    m_total: f64,
}

impl RedGreenIndex {
    /// Build from **real-valued** caps (m_i ≥ 0). Zeros are allowed.
    pub fn from_caps(m_per_dim: &[f64]) -> Self {
        let d = m_per_dim.len();
        let mut cum = Vec::with_capacity(d + 1);
        cum.push(0.0);
        let mut acc = 0.0f64;
        for &mi in m_per_dim {
            debug_assert!(mi >= 0.0, "caps must be non-negative");
            acc += mi;
            cum.push(acc);
        }
        Self {
            cum,
            d,
            m_total: acc,
        }
    }

    #[inline]
    pub fn d(&self) -> usize {
        self.d
    }

    #[inline]
    pub fn m_total(&self) -> f64 {
        self.m_total
    }

    /// Map r∈[0,M) to its component i and the left boundary cum[i].
    /// If r==M due to FP, clamp to nextafter(M, -∞).
    #[inline]
    pub fn comp_of(&self, mut r: f64) -> (usize, f64) {
        // Clamp r to [0, nextafter(M, -inf))
        if r >= self.m_total {
            // nextafter(M, -inf)
            r = f64::from_bits(self.m_total.to_bits() - 1);
        }
        // upper_bound: smallest j with cum[j] > r
        // then i = j - 1 so that cum[i] <= r < cum[i+1]
        let mut lo = 1usize;
        let mut hi = self.cum.len();
        while lo < hi {
            let mid = (lo + hi) >> 1;
            if self.cum[mid] > r {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        let i = lo - 1;
        (i, unsafe { *self.cum.get_unchecked(i) })
    }
}

/// Dense array of weights (length D) from a sparse vector.
#[inline]
fn dense_weights(d: usize, x: &[(u64, f64)]) -> Vec<f64> {
    let mut w = vec![0.0f64; d];
    for &(i, xi) in x {
        debug_assert!((i as usize) < d);
        if xi > 0.0 {
            w[i as usize] = xi;
        }
    }
    w
}

/// ERS (AAAI Algorithm 2): K independent fixed-length random sequences.
/// For each j in 0..K, scan r_{j,1},...,r_{j,L}; take first green. If none, mark E.
/// Then densify empties by rotating to the next non-empty bucket with a
/// per-j random offset (data-independent).
pub struct ErsWmh {
    index: RedGreenIndex,
    // tabulation generators
    t_u: Tab64Simple,   // U(0,1) for r_{j,t}
    t_id: Tab64Simple,  // ID from accepted draw r (via r.to_bits())
    t_rot: Tab32Simple, // offset for densification
    k: usize,
}

impl ErsWmh {
    /// `caps`: real-valued caps (tight upper bounds). `k`: number of hashes.
    pub fn new_mt(rng: &mut MtRng, caps: &[f64], k: u64) -> Self {
        let index = RedGreenIndex::from_caps(caps);
        let t_u = tab64_from_rng(rng);
        let t_id = tab64_from_rng(rng);
        let t_rot = tab32_from_rng(rng);
        Self {
            index,
            t_u,
            t_id,
            t_rot,
            k: k as usize,
        }
    }

    #[inline]
    fn is_green(&self, w_dense: &[f64], r: f64) -> bool {
        let (i, base) = self.index.comp_of(r);
        let xi = unsafe { *w_dense.get_unchecked(i) };
        r <= base + xi
    }

    /// `max_attempts` is interpreted as L (sequence length per hash position).
    /// If None, uses a moderate default (1024).
    pub fn sketch(&self, x: &[(u64, f64)], max_attempts: Option<u64>) -> Vec<Dart> {
        const L_DEFAULT: u32 = 1024;
        let l_per_hash: u32 = max_attempts.map(|v| v as u32).unwrap_or(L_DEFAULT);

        let w = dense_weights(self.index.d(), x);
        let m = self.index.m_total();

        // Degenerate: no mass or M==0 → deterministic, shared fallback
        let mass: f64 = w.iter().sum();
        if m == 0.0 || mass == 0.0 {
            let mut out = Vec::with_capacity(self.k);
            for j in 0..self.k {
                let fake = (self.t_rot.hash(j as u32) as u64) << 32 | (j as u64);
                out.push((fake, f64::INFINITY));
            }
            return out;
        }

        // One slot per hash position j
        let mut buckets: Vec<Option<(u64 /*id*/, u32 /*time*/)>> = vec![None; self.k];

        // Fixed-length sequences r_{j,t}; accept first green per j
        for j in 0..self.k {
            for t in 1..=l_per_hash {
                // key = (j, t)  → u ∈ [0,1) → r ∈ [0,M)
                let key = ((j as u64) << 32) ^ (t as u64);
                let mut u = to_unit(self.t_u.hash(key));
                if u >= 1.0 {
                    // guard against a theoretical 1.0 due to FP edge
                    u = f64::from_bits(0x3fefffffffffffff);
                }
                let r = m * u;

                if self.is_green(&w, r) {
                    // ID derived from the accepted *draw* r (ties collisions to the same r)
                    let id = self.t_id.hash(r.to_bits());
                    buckets[j] = Some((id, t));
                    break;
                }
            }
        }

        // If *all* buckets empty (very rare with decent L), fallback
        if buckets.iter().all(|b| b.is_none()) {
            let mut out = Vec::with_capacity(self.k);
            for j in 0..self.k {
                let fake = (self.t_rot.hash(j as u32) as u64) << 32 | (j as u64);
                out.push((fake, f64::INFINITY));
            }
            return out;
        }

        // Rotation densification: for each empty j, start from a per-j offset and
        // scan sequentially (mod k) until a non-empty bucket is found; copy it.
        for j in 0..self.k {
            if buckets[j].is_none() {
                // offset in {1,..,k-1}
                let offset = (self.t_rot.hash(j as u32) as usize % (self.k.saturating_sub(1)).max(1)) + 1;
                let mut idx = (j + offset) % self.k;

                // probe up to k-1 positions
                for _ in 0..(self.k - 1) {
                    if let Some(val) = buckets[idx] {
                        buckets[j] = Some(val);
                        break;
                    }
                    idx += 1;
                    if idx == self.k {
                        idx = 0;
                    }
                }

                // ultra-rare guard
                if buckets[j].is_none() {
                    let fake = (self.t_rot.hash(j as u32) as u64) << 32 | (j as u64);
                    buckets[j] = Some((fake, u32::MAX));
                }
            }
        }

        // Convert to (id, rank) = (hash_id, time as f64)
        let mut out = Vec::with_capacity(self.k);
        for j in 0..self.k {
            let (id, t) = buckets[j].unwrap();
            out.push((id, t as f64));
        }
        out
    }

    /// Uses default L (L_DEFAULT).
    #[inline]
    pub fn sketch_early_stop(&self, x: &[(u64, f64)]) -> Vec<Dart> {
        self.sketch(x, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng_utils::{mt_from_seed, MtRng};
    /// Generate a random weighted set with ids in [0, d)
    fn generate_weighted_set(d: usize, l0: u64, l1: f64, rng: &mut MtRng) -> Vec<(u64, f64)> {
        use std::collections::HashSet;
        let mut elements = HashSet::with_capacity(l0 as usize);
        while elements.len() < l0 as usize {
            let id = (rng.next_u64() as usize) % d;
            elements.insert(id as u64);
        }
        fn uniform01(rng: &mut MtRng) -> f64 {
            mt19937::gen_res53(rng)
        }
        let mut z: Vec<f64> = (0..(l0 - 1)).map(|_| uniform01(rng)).collect();
        z.push(1.0);
        z.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut prev = 0.0;
        let mut j = 0usize;
        let mut out: Vec<(u64, f64)> = Vec::with_capacity(l0 as usize);
        for idx in elements {
            let w = l1 * (z[j] - prev);
            out.push((idx, w.max(0.0)));
            prev = z[j];
            j += 1;
        }
        out.sort_by_key(|p| p.0);
        out
    }

    /// Generate Y from X with target overlap rel∈[0,1], ids in [0,d).
    /// Construction gives true weighted Jaccard J = rel / (2 - rel).
    fn generate_similar_weighted_set(
        d: usize,
        x: &[(u64, f64)],
        relative_overlap: f64,
        rng: &mut MtRng,
    ) -> Vec<(u64, f64)> {
        let free_id: u64 = loop {
            let cand = (rng.next_u64() as usize) % d;
            if x.binary_search_by_key(&(cand as u64), |p| p.0).is_err() {
                break cand as u64;
            }
        };
        let mut excess = 0.0;
        let mut y = Vec::with_capacity(x.len() + 1);
        for &(id, w) in x {
            let w_scaled = w * relative_overlap;
            excess += w - w_scaled;
            y.push((id, w_scaled.max(0.0)));
        }
        if excess > 0.0 {
            y.push((free_id, excess));
        }
        y.sort_by_key(|p| p.0);
        y
    }

    /// Build **real-valued** caps m_i that dominate all provided sets:
    /// m_i = max_s w_i(s)  (NO ceil, NO max(1)).
    fn caps_from_sets(d: usize, sets: &[&[(u64, f64)]]) -> Vec<f64> {
        let mut m = vec![0.0f64; d];
        for s in sets {
            for &(id, w) in *s {
                if w > 0.0 {
                    let idx = id as usize;
                    if w > m[idx] {
                        m[idx] = w;
                    }
                }
            }
        }
        m
    }

    #[test]
    fn ers_early_stop_fills_all_buckets() {
        let mut rng = mt_from_seed(1337);
        let d = 200_000usize;
        let k = 4096;

        // Base set
        let x = generate_weighted_set(d, 50_000, 10_000.0, &mut rng);

        // Build tight caps **from the data actually being sketched**
        let m = caps_from_sets(d, &[&x]);

        // ERS with data-consistent caps
        let ers = ErsWmh::new_mt(&mut rng, &m, k as u64);

        // Algorithm 2: per-hash sequence length L
        let l: u64 = 512;
        let sk = ers.sketch(&x, Some(l));

        assert_eq!(sk.len(), k);
    }

    #[test]
    fn ers_approximates_weighted_jaccard() {
        use crate::similarity::jaccard_similarity;

        let mut rng = mt_from_seed(8675309);
        let d = 200_000usize;
        let k = 4096;

        // Fixed per-hash sequence length (Algorithm 2)
        let l: u64 = 1024; // with tight caps, even smaller L often works

        let l0 = 50_000u64;
        let l1 = 10_000.0;
        let x = generate_weighted_set(d, l0, l1, &mut rng);
        let targets = [
            0.99, 0.96, 0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1,
            0.05, 0.01,
        ];

        for &rel in &targets {
            let y = generate_similar_weighted_set(d, &x, rel, &mut rng);
            let j_true = jaccard_similarity(&x, &y);
            println!("true weighted Jaccard: {:?}", j_true);
            // Tight caps must dominate BOTH vectors in the comparison
            let m_per_dim = caps_from_sets(d, &[&x, &y]);

            // Rebuild ERS for this pair with valid caps
            let ers = ErsWmh::new_mt(&mut rng, &m_per_dim, k as u64);

            // ERS (Alg.2): collision rate of per-bucket IDs
            let sk_x = ers.sketch(&x, Some(l));
            let sk_y = ers.sketch(&y, Some(l));

            let hits = sk_x.iter().zip(&sk_y).filter(|(a, b)| a.0 == b.0).count();
            let j_est = hits as f64 / k as f64;
            println!("estimated weighted Jaccard: {:?}", j_est);
            // σ-aware tolerance
            let sd = (j_true * (1.0 - j_true) / (k as f64)).sqrt();
            let tol = (3.2 * sd).max(1.25 / (k as f64).sqrt());
            let err = (j_true - j_est).abs();
            assert!(
                err <= tol,
                "rel={rel:.3}, true={j_true:.6}, est={j_est:.6}, err={err:.6}, tol={tol:.6}"
            );
        }
    }
}