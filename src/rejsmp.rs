//! Weighted MinHash via Rejection Sampling (RS) and Efficient Rejection Sampling (ERS).
//!
//! Implements:
//!   - RS (Shrivastava 2016): Algorithm 1/3 with constant-time ISGREEN via
//!     integer-to-component and component-to-M maps.
//!   - ERS (Li & Li 2021): single shared stream r_t, early stopping, and a safe
//!     densification fallback using tabulation hashing.
//!
//! Inputs: sparse weighted vector `&[(u64, f64)]` where id ∈ [0, D) and weight ≥ 0.
//! Randomness: purely via Tab32/Tab64 tabulation hashing (no stateful RNG required).
//! 
//! Important note: the max weight for a vector must be known in advance. Otherwise if w_i > w_max
//! for any vector, the ISGREEN test becomes wrong on that coordinate and your estimate biases high

use std::cmp::Ordering;

use rand_core::RngCore;
use tab_hash::{Tab32Simple, Tab64Simple};

use crate::rng_utils::MtRng;
use crate::hash_utils::{tab32_from_rng, tab64_from_rng, to_unit};

/// A single (id, rank) pair compatible with your DartMinHash plumbing.
pub type Dart = (u64, f64);

/// Integer line partition for the red–green test.
#[derive(Clone)]
pub struct RedGreenIndex {
    /// prefix sums M_i = sum_{j < i} m_j  (length D)
    comp_to_m: Vec<u64>,
    /// length M, entry j -> component i such that j in [M_i, M_{i+1})
    int_to_comp: Vec<u32>,
    /// total length M = sum_i m_i
    m_total: u64,
}

impl RedGreenIndex {
    pub fn from_m(m_per_dim: &[u32]) -> Self {
        let d = m_per_dim.len();
        let mut comp_to_m = Vec::with_capacity(d);
        let mut int_to_comp = Vec::new();
        int_to_comp.reserve(m_per_dim.iter().map(|&mi| mi as usize).sum());

        let mut prefix: u64 = 0;
        for (i, &mi) in m_per_dim.iter().enumerate() {
            comp_to_m.push(prefix);
            for _ in 0..mi {
                int_to_comp.push(i as u32);
            }
            prefix += mi as u64;
        }
        Self { comp_to_m, int_to_comp, m_total: prefix }
    }

    #[inline]
    pub fn m_total(&self) -> u64 { self.m_total }

    /// Map r∈[0,M) to its component, clamped at M-1 to guard FP edge cases.
    #[inline]
    pub fn comp_of(&self, r: f64) -> (u32, u64) {
        let idx = ((r as u64).min(self.m_total.saturating_sub(1))) as usize;
        let i = self.int_to_comp[idx] as usize;
        (i as u32, self.comp_to_m[i])
    }
}

/// Dense array of weights (length D) from a sparse vector.
#[inline]
fn dense_weights(d: usize, x: &[(u64, f64)]) -> Vec<f64> {
    let mut w = vec![0.0f64; d];
    for &(i, xi) in x {
        debug_assert!((i as usize) < d);
        if xi > 0.0 { w[i as usize] = xi; }
    }
    w
}

/// Counter-based U(0,1) from tabulation: (seed, counter) -> 64-bit -> [0,1).
#[inline]
fn u01_from_tab(tab: &Tab64Simple, seed: u64, counter: u64) -> f64 {
    let z = tab.hash(seed ^ counter);
    to_unit(z)
}

/// Rejection Samplng: For each hash, scan a **shared**, x-independent
/// sequence r_t = f(seed, t) over [0,M) and return the first accepted r*’s identity.
/// Two sets collide iff they accept the same r* → unbiased for Jaccard.
pub struct RsWmh {
    d: usize,
    index: RedGreenIndex,
    t_u: Tab64Simple,   // U(0,1) generator
    t_sig: Tab64Simple, // signature from r_bits
    seeds: Vec<u64>,    // per-hash seeds
}

impl RsWmh {
    /// m_per_dim: caps m_i (use 1 if normalized). k: number of hashes.
    pub fn new_mt(rng: &mut MtRng, m_per_dim: &[u32], k: usize) -> Self {
        let d = m_per_dim.len();
        let index = RedGreenIndex::from_m(m_per_dim);
        let t_u  = tab64_from_rng(rng);
        let t_sig = tab64_from_rng(rng);

        let mut seeds = Vec::with_capacity(k);
        for _ in 0..k { seeds.push(rng.next_u64()); }
        Self { d, index, t_u, t_sig, seeds }
    }

    #[inline]
    fn is_green(&self, w_dense: &[f64], r: f64) -> bool {
        let (i, mi) = self.index.comp_of(r);
        let xi = unsafe { *w_dense.get_unchecked(i as usize) };
        r <= (mi as f64) + xi
    }

    /// First accepted r_t identity for one hash.
    #[inline]
    fn one_id(&self, w_dense: &[f64], seed: u64) -> u64 {
        let m = self.index.m_total() as f64;
        let mut t: u64 = 1;
        loop {
            let u = u01_from_tab(&self.t_u, seed, t);
            let r = m * u;
            if self.is_green(w_dense, r) {
                return self.t_sig.hash(r.to_bits());
            }
            t += 1;
        }
    }

    /// RS signature as **k IDs** (collision rate estimates J).
    pub fn sketch_ids(&self, x: &[(u64, f64)]) -> Vec<u64> {
        let w = dense_weights(self.d, x);
        let mut out = Vec::with_capacity(self.seeds.len());
        for &seed in &self.seeds {
            out.push(self.one_id(&w, seed));
        }
        out
    }

    /// Optional diagnostic: geometric “trial counts” (first accepted t).
    pub fn sketch_counts(&self, x: &[(u64, f64)]) -> Vec<u16> {
        let w = dense_weights(self.d, x);
        let m = self.index.m_total() as f64;
        let mut out = Vec::with_capacity(self.seeds.len());
        for &seed in &self.seeds {
            let mut t: u64 = 1;
            loop {
                let u = u01_from_tab(&self.t_u, seed, t);
                let r = m * u;
                if self.is_green(&w, r) {
                    out.push(t as u16);
                    break;
                }
                t += 1;
            }
        }
        out
    }
}

/// Efficient Rejection Sampling: define a **shared global order** r_t = f(0, t) on [0,M).
/// For a set x, accept r_t if green(x). Route each accepted r* to a bucket
/// using a tabulated key from r_bits; rank = t (lower is better). Early-stop
/// once all k buckets have ≥1 candidate; otherwise densify safely.
pub struct ErsWmh {
    index: RedGreenIndex,
    d: usize,
    // tabulation “generators”
    t_u: Tab64Simple,      // for U(0,1) → r
    t_key: Tab64Simple,    // key id from r_bits
    t_bucket: Tab64Simple, // bucket routing from key id
    t_frac: Tab64Simple,   // fractional tie-breaker from r_bits
    t_rot: Tab32Simple,    // rotation seed for densification
    k: u64,
}

#[derive(Clone, Copy)]
struct BucketKey {
    time: u64,   // global attempt index; lower is better
    frac: u32,   // stable fractional tiebreaker
    hash_id: u64 // stable identity derived from r*
}
impl BucketKey { #[inline] fn rank(&self) -> (u64, u32) { (self.time, self.frac) } }
impl PartialEq for BucketKey { fn eq(&self, o: &Self) -> bool { self.rank()==o.rank() } }
impl Eq for BucketKey {}
impl PartialOrd for BucketKey { fn partial_cmp(&self, o: &Self)->Option<Ordering>{Some(self.cmp(o))} }
impl Ord for BucketKey { fn cmp(&self, o: &Self)->Ordering{ self.rank().cmp(&o.rank()) } }

impl ErsWmh {
    pub fn new_mt(rng: &mut MtRng, m_per_dim: &[u32], k: u64) -> Self {
        let index = RedGreenIndex::from_m(m_per_dim);
        let t_u = tab64_from_rng(rng);
        let t_key = tab64_from_rng(rng);
        let t_bucket = tab64_from_rng(rng);
        let t_frac = tab64_from_rng(rng);
        let t_rot = tab32_from_rng(rng);
        Self { index, d: m_per_dim.len(), t_u, t_key, t_bucket, t_frac, t_rot, k }
    }

    #[inline]
    fn is_green(&self, w_dense: &[f64], r: f64) -> bool {
        let (i, mi) = self.index.comp_of(r);
        let xi = unsafe { *w_dense.get_unchecked(i as usize) };
        r <= (mi as f64) + xi
    }

    #[inline] fn route_bucket(&self, key_id: u64) -> usize {
        (self.t_bucket.hash(key_id) % self.k) as usize
    }
    #[inline] fn frac32(&self, r_bits: u64) -> u32 {
        (self.t_frac.hash(r_bits) & 0xFFFF_FFFF) as u32
    }

    /// One ERS sketch: k buckets; minimal (time, frac) per bucket wins.
    /// `max_attempts`: optional cap on attempts before densification.
    pub fn sketch(&self, x: &[(u64, f64)], max_attempts: Option<u64>) -> Vec<Dart> {
        let w = dense_weights(self.d, x);
        let m = self.index.m_total() as f64;
        let mut buckets: Vec<Option<BucketKey>> = vec![None; self.k as usize];

        let wanted = self.k as usize;
        let mut filled = 0usize;
        let mut attempts: u64 = 0;

        // Early-stopping loop
        loop {
            attempts += 1;
            let u = u01_from_tab(&self.t_u, 0, attempts);
            let r = m * u;

            if self.is_green(&w, r) {
                let r_bits = r.to_bits();
                let key_id = self.t_key.hash(r_bits);
                let b = self.route_bucket(key_id);
                let key = BucketKey { time: attempts, frac: self.frac32(r_bits), hash_id: key_id };
                match &mut buckets[b] {
                    None => { buckets[b] = Some(key); filled += 1;
                              if filled == wanted && max_attempts.is_none() { break; } }
                    Some(best) => { if key < *best { *best = key; } }
                }
            }

            if let Some(cap) = max_attempts { if attempts >= cap { break; } }
        }

        // Safe densification:
        if filled == 0 {
            // Deterministic fallback: produce k fake slots
            let mut out = Vec::with_capacity(wanted);
            for j in 0..wanted {
                let fake = (self.t_rot.hash(j as u32) as u64) << 32 | (j as u64);
                out.push((fake, u64::MAX as f64));
            }
            return out;
        }
        if filled < wanted {
            for j in 0..wanted {
                if buckets[j].is_none() {
                    // rotate to next non-empty
                    let mut off = (self.t_rot.hash(j as u32) as usize) % wanted;
                    if off == 0 { off = 1; }
                    let key: BucketKey = {
                        let mut t = 0usize;
                        loop {
                            let jj = (j + off * (t + 1)) % wanted;
                            if let Some(kv) = buckets[jj] { break kv; }
                            t += 1;
                        }
                    };
                    buckets[j] = Some(key);
                }
            }
        }

        // Convert to (id, rank) = (hash_id, time as f64)
        let mut out = Vec::with_capacity(wanted);
        for j in 0..wanted {
            let key = buckets[j].expect("bucket must be filled after densification");
            out.push((key.hash_id, key.time as f64));
        }
        out
    }

    /// Convenience: early stop with no attempt cap.
    pub fn sketch_early_stop(&self, x: &[(u64, f64)]) -> Vec<Dart> {
        self.sketch(x, None)
    }

    /// 1-bit ERS sketch from the key ids (LSB).
    pub fn onebit(&self, x: &[(u64, f64)]) -> Vec<bool> {
        self.sketch_early_stop(&x).into_iter().map(|(id, _)| (id & 1) == 1).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng_utils::{mt_from_seed, MtRng};
    use crate::similarity::jaccard_similarity;

    /// Generate a random weighted set with ids in [0, d)
    fn generate_weighted_set(d: usize, l0: u64, l1: f64, rng: &mut MtRng) -> Vec<(u64, f64)> {
        use std::collections::HashSet;
        let mut elements = HashSet::with_capacity(l0 as usize);
        while elements.len() < l0 as usize {
            let id = (rng.next_u64() as usize) % d;
            elements.insert(id as u64);
        }
        fn uniform01(rng: &mut MtRng) -> f64 { mt19937::gen_res53(rng) }
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
        if excess > 0.0 { y.push((free_id, excess)); }
        y.sort_by_key(|p| p.0);
        y
    }

    /// Build per-dimension caps m_i that dominate all provided sets:
    /// m_i = max(1, ceil(max_s w_i(s))) for any set s in `sets`.
    fn caps_from_sets(d: usize, sets: &[&[(u64, f64)]]) -> Vec<u32> {
        let mut m = vec![0u32; d];
        for s in sets {
            for &(id, w) in *s {
                if w > 0.0 {
                    let idx = id as usize;
                    let cap = (w.ceil() as u32).max(1);
                    if cap > m[idx] { m[idx] = cap; }
                }
            }
        }
        m
    }

    #[test]
    fn rs_counts_are_small() {
        let mut rng = mt_from_seed(7);
        let d = 10_000usize;
        // Build caps from x so that x_i <= m_i
        let x = vec![(1u64, 0.7), (123u64, 0.4), (9999u64, 1.8)];
        let m = caps_from_sets(d, &[&x]);
        let k = 1024;

        let rs = RsWmh::new_mt(&mut rng, &m, k);
        let h = rs.sketch_counts(&x);
        assert_eq!(h.len(), k);
        assert!(h.iter().all(|&v| v > 0));
    }

    #[test]
    fn ers_early_stop_fills_all_buckets() {
        let mut rng = mt_from_seed(1337);
        let d = 200_000usize;
        let k = 4096;

        let x = generate_weighted_set(d, 50_000, 10_000.0, &mut rng);
        let m = caps_from_sets(d, &[&x]);           // <- caps from x
        let ers = ErsWmh::new_mt(&mut rng, &m, k as u64);

        let sk = ers.sketch_early_stop(&x);
        assert_eq!(sk.len(), k);
    }

    #[test]
    fn rs_approximates_weighted_jaccard() {
        let mut rng = mt_from_seed(4242);
        let d = 200_000usize;
        let k = 4096;

        let l0 = 50_000u64;
        let l1 = 10_000.0;
        let x = generate_weighted_set(d, l0, l1, &mut rng);
        let targets = [0.99, 0.96,0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01];

        for &rel in &targets {
            let y = generate_similar_weighted_set(d, &x, rel, &mut rng);
            let j_true = jaccard_similarity(&x, &y);
            println!("true weighted Jaccard: {:?}", j_true);
            // Build caps FROM THIS PAIR (or, in production, dataset-wide).
            let m = caps_from_sets(d, &[&x, &y]);
            let rs = RsWmh::new_mt(&mut rng, &m, k);

            let sig_x = rs.sketch_ids(&x);
            let sig_y = rs.sketch_ids(&y);

            let mut hits = 0usize;
            for i in 0..k { if sig_x[i] == sig_y[i] { hits += 1; } }
            let j_est = hits as f64 / k as f64;
            println!("estimated weighted Jaccard: {:?}", j_est);
            let sd = (j_true * (1.0 - j_true) / (k as f64)).sqrt();
            let tol = (3.0 * sd).max(1.1 / (k as f64).sqrt()); // slightly conservative

            let err = (j_true - j_est).abs();
            assert!(
                err <= tol,
                "rel={rel:.3}, true={j_true:.6}, est={j_est:.6}, err={err:.6}, tol={tol:.6}"
            );
        }
    }

    #[test]
    fn ers_approximates_weighted_jaccard() {
        let mut rng = mt_from_seed(8675309);
        let d = 200_000usize;
        let k = 4096;

        let l0 = 50_000u64;
        let l1 = 10_000.0;
        let x = generate_weighted_set(d, l0, l1, &mut rng);
        let targets = [0.99, 0.96,0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01];

        for &rel in &targets {
            let y = generate_similar_weighted_set(d, &x, rel, &mut rng);
            let j_true = jaccard_similarity(&x, &y);
            println!("true weighted Jaccard: {:?}", j_true);
            // Caps from this pair (or dataset-wide in production)
            let m = caps_from_sets(d, &[&x, &y]);
            let ers = ErsWmh::new_mt(&mut rng, &m, k as u64);

            // Collision rate of key ids across buckets
            let sk_x = ers.sketch_early_stop(&x);
            let sk_y = ers.sketch_early_stop(&y);

            let mut hits = 0usize;
            for i in 0..k as usize { if sk_x[i].0 == sk_y[i].0 { hits += 1; } }
            let j_est = hits as f64 / k as f64;
            println!("estimated weighted Jaccard: {:?}", j_est);
            let sd = (j_true * (1.0 - j_true) / (k as f64)).sqrt();
            let tol = (3.0 * sd).max(1.1 / (k as f64).sqrt());

            let err = (j_true - j_est).abs();
            assert!(
                err <= tol,
                "rel={rel:.3}, true={j_true:.6}, est={j_est:.6}, err={err:.6}, tol={tol:.6}"
            );
        }
    }
}