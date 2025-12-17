//! Weighted MinHash via Efficient Rejection Sampling (ERS).
//!
//! Implements:
//!   - ERS (Li & Li 2021): k independent fixed-length sequences
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
//!
//! PERFORMANCE NOTE (this version):
//! - Removes the hot O(log D) binary search over prefix-sums for each draw by using
//!   a Walker alias table to sample the interval i in O(1).
//! - Keeps your original semantics for ID hashing: id = hash(r.to_bits()) where
//!   r = base[i] + off, off ~ Uniform(0, m_i).

use tab_hash::{Tab32Simple, Tab64Simple};
use rand_core::RngCore;

use crate::hash_utils::{tab32_from_rng, tab64_from_rng, to_unit};
use crate::rng_utils::MtRng;

use std::cell::RefCell;

/// A single (id, rank) pair compatible with your DartMinHash plumbing.
pub type Dart = (u64, f64);

/// Continuous (real-valued) line partition for the red–green test.
///
/// We store:
/// - base[i] = sum_{h<i} m_h  (left boundary)
/// - cap[i]  = m_i
/// And a Walker alias table to sample i with P(i)=m_i/M in O(1).
#[derive(Clone)]
pub struct RedGreenIndex {
    base: Vec<f64>,
    cap: Vec<f64>,
    d: usize,
    m_total: f64,

    // Walker alias table for discrete distribution p_i = cap[i]/m_total
    prob: Vec<f64>,  // in [0,1]
    alias: Vec<u32>, // in [0,d)
}

impl RedGreenIndex {
    /// Build from **real-valued** caps (m_i ≥ 0). Zeros are allowed.
    pub fn from_caps(m_per_dim: &[f64]) -> Self {
        let d = m_per_dim.len();

        // base prefix sums
        let mut base = Vec::with_capacity(d);
        let mut acc = 0.0f64;
        for &mi in m_per_dim {
            debug_assert!(mi >= 0.0, "caps must be non-negative");
            base.push(acc);
            acc += mi;
        }
        let m_total = acc;

        // copy caps
        let cap = m_per_dim.to_vec();

        let mut prob = vec![0.0f64; d];
        let mut alias = vec![0u32; d];

        // Degenerate cases
        if d == 0 || m_total == 0.0 {
            return Self { base, cap, d, m_total, prob, alias };
        }

        // Walker alias build
        // scaled probabilities: q_i = p_i * d = (cap[i]/m_total) * d
        let mut q: Vec<f64> = cap
            .iter()
            .map(|&mi| (mi / m_total) * (d as f64))
            .collect();

        let mut small = Vec::<usize>::new();
        let mut large = Vec::<usize>::new();
        for (i, &qi) in q.iter().enumerate() {
            if qi < 1.0 { small.push(i); } else { large.push(i); }
        }

        while let (Some(s), Some(l)) = (small.pop(), large.pop()) {
            prob[s] = q[s];      // < 1
            alias[s] = l as u32; // redirect

            q[l] = (q[l] + q[s]) - 1.0;
            if q[l] < 1.0 {
                small.push(l);
            } else {
                large.push(l);
            }
        }

        // leftovers
        for i in small.into_iter().chain(large.into_iter()) {
            prob[i] = 1.0;
            alias[i] = i as u32;
        }

        Self { base, cap, d, m_total, prob, alias }
    }

    #[inline]
    pub fn d(&self) -> usize { self.d }

    #[inline]
    pub fn m_total(&self) -> f64 { self.m_total }

    #[inline]
    pub fn base_of(&self, i: usize) -> f64 {
        unsafe { *self.base.get_unchecked(i) }
    }

    #[inline]
    pub fn cap_of(&self, i: usize) -> f64 {
        unsafe { *self.cap.get_unchecked(i) }
    }

    /// Sample an interval i with P(i)=cap[i]/M, plus offset off ∈ [0,cap[i]).
    /// Uses only tab-hash-derived uniforms (stateless).
    #[inline]
    pub fn sample_interval_and_offset(&self, t_u: &Tab64Simple, key: u64) -> (usize, f64) {
        debug_assert!(self.d > 0);

        // u0 chooses the column in [0,d)
        let mut u0 = to_unit(t_u.hash(key));
        if u0 >= 1.0 {
            u0 = f64::from_bits(0x3fefffffffffffff); // < 1.0
        }
        let mut col = (u0 * (self.d as f64)) as usize;
        if col >= self.d {
            col = self.d - 1;
        }

        // u1 decides alias/keep
        let mut u1 = to_unit(t_u.hash(key ^ 0x9e37_79b9_7f4a_7c15));
        if u1 >= 1.0 {
            u1 = f64::from_bits(0x3fefffffffffffff);
        }
        let i = if u1 < unsafe { *self.prob.get_unchecked(col) } {
            col
        } else {
            unsafe { *self.alias.get_unchecked(col) as usize }
        };

        // u2 chooses offset within interval i
        let mut u2 = to_unit(t_u.hash(key ^ 0xbf58_476d_1ce4_e5b9));
        if u2 >= 1.0 {
            u2 = f64::from_bits(0x3fefffffffffffff);
        }
        let off = self.cap_of(i) * u2;
        (i, off)
    }
}

/// Per-thread dense scratch to avoid allocating/zeroing a length-D vector per sample.
/// We only touch indices present in x, and only clear those indices afterward.
#[derive(Default)]
struct DenseScratch {
    w: Vec<f64>,
    touched: Vec<usize>,
}

impl DenseScratch {
    #[inline]
    fn ensure_len(&mut self, d: usize) {
        if self.w.len() < d {
            self.w.resize(d, 0.0);
        }
    }

    /// Populate dense weights from sparse `x` and return `mass` computed with the same
    /// semantics as a dense vector sum (handles duplicate ids by overwrite).
    #[inline]
    fn fill_from_sparse_and_mass(&mut self, d: usize, x: &[(u64, f64)]) -> f64 {
        self.ensure_len(d);
        self.touched.clear();

        let mut mass = 0.0f64;

        for &(i_u64, xi) in x {
            let idx = i_u64 as usize;
            debug_assert!(idx < d);

            if xi <= 0.0 {
                continue;
            }

            let old = self.w[idx];
            if old == 0.0 {
                self.touched.push(idx);
                self.w[idx] = xi;
                mass += xi;
            } else {
                // duplicate id: overwrite semantics
                self.w[idx] = xi;
                mass += xi - old;
            }
        }

        mass
    }

    #[inline]
    fn clear_touched(&mut self) {
        for &idx in &self.touched {
            self.w[idx] = 0.0;
        }
        self.touched.clear();
    }
}

thread_local! {
    static ERS_SCRATCH: RefCell<DenseScratch> = RefCell::new(DenseScratch::default());
}

/// ERS (AAAI Algorithm 2): k independent fixed-length random sequences.
/// For each j in 0..k, scan r_{j,1},...,r_{j,L}; take first green. If none, mark empty.
/// Then densify empties by rotating to a non-empty bucket with a per-j random offset.
pub struct ErsWmh {
    index: RedGreenIndex,
    // tabulation generators
    t_u: Tab64Simple,   // U(0,1) for draws
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
        Self { index, t_u, t_id, t_rot, k: k as usize }
    }

    #[inline]
    fn is_green_offset(&self, w_dense: &[f64], i: usize, off: f64) -> bool {
        // green iff off <= x_i (since r = base[i] + off and green region is [base, base + x_i])
        let xi = unsafe { *w_dense.get_unchecked(i) };
        off <= xi
    }

    /// `max_attempts` is interpreted as L (sequence length per hash position).
    /// If None, uses a moderate default (1024).
    pub fn sketch(&self, x: &[(u64, f64)], max_attempts: Option<u64>) -> Vec<Dart> {
        const L_DEFAULT: u32 = 1024;
        let l_per_hash: u32 = max_attempts.map(|v| v as u32).unwrap_or(L_DEFAULT);

        let d = self.index.d();
        let m = self.index.m_total();

        // One slot per hash position j
        let mut buckets: Vec<Option<(u64 /*id*/, u32 /*time*/)>> = vec![None; self.k];

        // Use per-thread scratch to avoid O(D) alloc/zero and O(D) mass sum.
        let mut out: Option<Vec<Dart>> = None;

        ERS_SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();

            // Fill dense vector (only touched indices) and compute mass with dense semantics.
            let mass = scratch.fill_from_sparse_and_mass(d, x);

            // Degenerate: no mass or M==0 → deterministic fallback
            if m == 0.0 || mass == 0.0 || d == 0 {
                let mut fallback = Vec::with_capacity(self.k);
                for j in 0..self.k {
                    let fake = (self.t_rot.hash(j as u32) as u64) << 32 | (j as u64);
                    fallback.push((fake, f64::INFINITY));
                }
                scratch.clear_touched();
                out = Some(fallback);
                return;
            }

            let w = &scratch.w;

            // Fixed-length sequences; accept first green per j
            for j in 0..self.k {
                for t in 1..=l_per_hash {
                    // key = (j, t)
                    let key = ((j as u64) << 32) ^ (t as u64);

                    // O(1) interval sample + offset
                    let (i, off) = self.index.sample_interval_and_offset(&self.t_u, key);

                    if self.is_green_offset(w, i, off) {
                        // Reconstruct r so ID hashing matches the previous definition.
                        let r = self.index.base_of(i) + off;
                        let id = self.t_id.hash(r.to_bits());
                        buckets[j] = Some((id, t));
                        break;
                    }
                }
            }

            // If *all* buckets empty (very rare with decent L), fallback
            if buckets.iter().all(|b| b.is_none()) {
                let mut fallback = Vec::with_capacity(self.k);
                for j in 0..self.k {
                    let fake = (self.t_rot.hash(j as u32) as u64) << 32 | (j as u64);
                    fallback.push((fake, f64::INFINITY));
                }
                scratch.clear_touched();
                out = Some(fallback);
                return;
            }

            // Rotation densification
            for j in 0..self.k {
                if buckets[j].is_none() {
                    // offset in {1,..,k-1}
                    let offset =
                        (self.t_rot.hash(j as u32) as usize % (self.k.saturating_sub(1)).max(1)) + 1;
                    let mut idx = (j + offset) % self.k;

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
            let mut result = Vec::with_capacity(self.k);
            for j in 0..self.k {
                let (id, t) = buckets[j].unwrap();
                result.push((id, t as f64));
            }

            scratch.clear_touched();
            out = Some(result);
        });

        out.expect("ERS_SCRATCH closure must set out")
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