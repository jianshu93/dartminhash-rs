//! Weighted MinHash via Efficient Rejection Sampling (ERS).
//!
//! Implements:
//!   - ERS (Li & Li 2021): k independent fixed-length sequences
//!     r_{j,1..L} per hash position j; take the first green if any, otherwise mark
//!     empty; then densify empties with a shared, data-independent candidate
//!     sequence as specified by Algorithm 3.
//!
//! Inputs: sparse weighted vector `&[(u64, f64)]` where id ∈ [0, D) and weight ≥ 0.
//! Randomness: purely via Tab32/Tab64 tabulation hashing (no stateful RNG required).
//! By default this uses simple tabulation; with the `mixed_tab` feature it uses
//! mixed tabulation.
//!
//! IMPORTANT: Caps `m_i` are **real-valued** (`f64`) and should be set to the
//! *tight* per-dimension maxima across the dataset: `m_i = max_s x_i(s)`.
//! Using tight caps reduces total M = sum_i m_i, increases acceptance probability,
//! and lets you use much smaller L in ERS.

use crate::hash_utils::*;
use crate::rng_utils::MtRng;

use std::cell::RefCell;

#[cfg(feature = "mixed_tab")]
type Tab64Ers = tab_hash::Tab64Mixed;
#[cfg(not(feature = "mixed_tab"))]
type Tab64Ers = tab_hash::Tab64Simple;

#[cfg(feature = "mixed_tab")]
fn tab64_ers_from_rng(rng: &mut MtRng) -> Tab64Ers {
    mixed_tab64_from_rng(rng)
}

#[cfg(not(feature = "mixed_tab"))]
fn tab64_ers_from_rng(rng: &mut MtRng) -> Tab64Ers {
    tab64_from_rng(rng)
}

/// A single (id, rank) pair compatible with your DartMinHash plumbing.
pub type Dart = (u64, f64);

/// Continuous (real-valued) line partition for the red–green test.
///
/// We store:
/// - base[i] = sum_{h<i} m_h  (left boundary)
/// - cap[i]  = m_i
#[derive(Clone)]
pub struct RedGreenIndex {
    base: Vec<f64>,
    cap: Vec<f64>,
    d: usize,
    m_total: f64,
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
        Self {
            base,
            cap: m_per_dim.to_vec(),
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

    #[inline]
    pub fn base_of(&self, i: usize) -> f64 {
        unsafe { *self.base.get_unchecked(i) }
    }

    #[inline]
    pub fn cap_of(&self, i: usize) -> f64 {
        unsafe { *self.cap.get_unchecked(i) }
    }

    /// Map `r` in `[0, M)` to its cap interval and left boundary.
    #[inline]
    fn comp_of(&self, mut r: f64) -> (usize, f64) {
        if r >= self.m_total {
            r = f64::from_bits(self.m_total.to_bits() - 1);
        }

        // upper_bound: smallest j whose interval base is greater than r.
        // Using upper_bound also skips zero-width intervals at duplicate bases.
        let mut lo = 1usize;
        let mut hi = self.base.len();
        while lo < hi {
            let mid = (lo + hi) >> 1;
            if self.base[mid] > r {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        let i = lo - 1;
        (i, unsafe { *self.base.get_unchecked(i) })
    }

    /// Sample an interval i with P(i)=cap[i]/M, plus offset off ∈ [0,cap[i]).
    /// Uses only tab-hash-derived uniforms (stateless).
    #[cfg(not(feature = "mixed_tab"))]
    #[inline]
    pub fn sample_interval_and_offset(
        &self,
        t_u: &tab_hash::Tab64Simple,
        key: u64,
    ) -> (usize, f64) {
        self.sample_interval_and_offset_impl(t_u, key)
    }

    /// Sample an interval i with P(i)=cap[i]/M, plus offset off ∈ [0,cap[i]).
    /// Uses only tab-hash-derived uniforms (stateless).
    #[cfg(feature = "mixed_tab")]
    #[inline]
    pub fn sample_interval_and_offset(&self, t_u: &tab_hash::Tab64Mixed, key: u64) -> (usize, f64) {
        self.sample_interval_and_offset_impl(t_u, key)
    }

    #[inline]
    fn sample_interval_and_offset_impl(&self, t_u: &Tab64Ers, key: u64) -> (usize, f64) {
        debug_assert!(self.d > 0);
        debug_assert!(self.m_total > 0.0);

        let mut u = to_unit(t_u.hash(key));
        if u >= 1.0 {
            u = f64::from_bits(0x3fefffffffffffff);
        }
        let r = self.m_total * u;
        let (i, base) = self.comp_of(r);
        (i, r - base)
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
/// Then densify empties with Algorithm 3's shared candidate sequence.
pub struct ErsWmh {
    index: RedGreenIndex,
    // tabulation generators
    t_u: Tab64Ers,             // U(0,1) for r_{j,t}
    t_id: Tab64Ers,            // ID from accepted draw r (via r.to_bits())
    t_dense_bucket: Tab64Ers,  // densification hash: output bucket
    t_dense_attempt: Tab64Ers, // densification hash: retry number
    k: usize,
}

impl ErsWmh {
    /// `caps`: real-valued caps (tight upper bounds). `k`: number of hashes.
    pub fn new_mt(rng: &mut MtRng, caps: &[f64], k: u64) -> Self {
        let index = RedGreenIndex::from_caps(caps);
        // Keep the original accepted-draw seed order for reproducibility.
        let t_u = tab64_ers_from_rng(rng);
        let t_id = tab64_ers_from_rng(rng);
        let t_dense_bucket = tab64_ers_from_rng(rng);
        let t_dense_attempt = tab64_ers_from_rng(rng);
        Self {
            index,
            t_u,
            t_id,
            t_dense_bucket,
            t_dense_attempt,
            k: k as usize,
        }
    }

    #[inline]
    fn is_green(&self, w_dense: &[f64], r: f64) -> bool {
        let (i, base) = self.index.comp_of(r);
        let xi = unsafe { *w_dense.get_unchecked(i) };
        r <= base + xi
    }

    #[inline]
    fn index_from_hash(hash: u64, len: usize) -> usize {
        // Lemire's multiply-high range reduction avoids the modulo bias that is
        // otherwise most visible when k is small.
        ((hash as u128 * len as u128) >> 64) as usize
    }

    /// Algorithm 3's H(j, attempt), represented by two independently seeded
    /// tabulation hashes. The same candidate sequence is therefore used for a
    /// given output bucket in every vector sketched by this ERS instance.
    #[inline]
    fn densification_candidate(&self, j: usize, attempt: usize) -> usize {
        let hash = self.t_dense_bucket.hash(j as u64) ^ self.t_dense_attempt.hash(attempt as u64);
        Self::index_from_hash(hash, self.k)
    }

    fn densification_source(&self, original_nonempty: &[bool], j: usize) -> usize {
        debug_assert_eq!(original_nonempty.len(), self.k);
        debug_assert!(original_nonempty.iter().any(|&present| present));

        // Independent uniform candidates need k/n_nonempty probes on average.
        // 32k probes makes failure less likely than exp(-32), even with only
        // one donor, while still bounding runtime for a defective hash family.
        let max_probes = self.k.saturating_mul(32).max(32);
        for attempt in 1..=max_probes {
            let candidate = self.densification_candidate(j, attempt);
            if original_nonempty[candidate] {
                return candidate;
            }
        }

        // Deterministic, consistent fallback: give every original donor a
        // candidate-specific random priority and take the minimum. This path is
        // practically unreachable, but unlike a linear scan it remains
        // independent of the donor layout.
        let bucket_hash = self.t_dense_bucket.hash(j as u64);
        original_nonempty
            .iter()
            .enumerate()
            .filter(|(_, present)| **present)
            .min_by_key(|(candidate, _)| bucket_hash ^ self.t_dense_attempt.hash(*candidate as u64))
            .map(|(candidate, _)| candidate)
            .expect("at least one original ERS bucket is non-empty")
    }

    fn densify(&self, buckets: &mut [Option<(u64, u32)>]) {
        let original_nonempty: Vec<bool> = buckets.iter().map(Option::is_some).collect();

        for j in 0..self.k {
            if !original_nonempty[j] {
                let source = self.densification_source(&original_nonempty, j);
                buckets[j] = buckets[source];
            }
        }
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
                    let fake = self.t_dense_bucket.hash(j as u64) ^ j as u64;
                    fallback.push((fake, f64::INFINITY));
                }
                scratch.clear_touched();
                out = Some(fallback);
                return;
            }

            let w = &scratch.w;

            // Fixed-length sequences r_{j,t}; accept first green per j.
            for j in 0..self.k {
                for t in 1..=l_per_hash {
                    // One uniform draw over the complete red-green line [0, M).
                    let key = ((j as u64) << 32) ^ (t as u64);
                    let mut u = to_unit(self.t_u.hash(key));
                    if u >= 1.0 {
                        u = f64::from_bits(0x3fefffffffffffff);
                    }
                    let r = m * u;

                    if self.is_green(w, r) {
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
                    let fake = self.t_dense_bucket.hash(j as u64) ^ j as u64;
                    fallback.push((fake, f64::INFINITY));
                }
                scratch.clear_touched();
                out = Some(fallback);
                return;
            }

            self.densify(&mut buckets);

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
    use rand_core::RngCore;

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
        let mut ids: Vec<u64> = elements.into_iter().collect();
        ids.sort_unstable();
        for idx in ids {
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
    fn ers_sparse_entry_order_is_not_part_of_the_vector() {
        let caps = vec![1.0, 0.8, 0.6, 0.4];
        let ordered = vec![(0, 0.75), (1, 0.5), (2, 0.25), (3, 0.125)];
        let reversed = ordered.iter().copied().rev().collect::<Vec<_>>();
        let mut hash_rng = mt_from_seed(0xe255_2001);
        let ers = ErsWmh::new_mt(&mut hash_rng, &caps, 4096);

        assert_eq!(
            ers.sketch(&ordered, Some(16)),
            ers.sketch(&reversed, Some(16))
        );
    }

    #[test]
    fn algorithm3_selects_original_donors_uniformly() {
        let k = 4096usize;
        let mut hash_rng = mt_from_seed(0xe255_2002);
        let ers = ErsWmh::new_mt(&mut hash_rng, &[1.0], k as u64);
        let mut original_nonempty = vec![false; k];
        original_nonempty[0] = true;
        original_nonempty[1] = true;

        let mut counts = [0usize; 2];
        for j in 2..k {
            counts[ers.densification_source(&original_nonempty, j)] += 1;
        }

        let expected = (k - 2) as f64 / 2.0;
        for (source, observed) in counts.into_iter().enumerate() {
            let relative_error = (observed as f64 - expected).abs() / expected;
            assert!(
                relative_error < 0.06,
                "source={source}, observed={observed}, expected={expected}, relative_error={relative_error}"
            );
        }
    }

    #[test]
    fn algorithm3_never_uses_a_densified_bucket_as_a_donor() {
        let mut hash_rng = mt_from_seed(0xe255_2003);
        let ers = ErsWmh::new_mt(&mut hash_rng, &[1.0], 128);
        let donor_a = (0xa11c_e001, 3);
        let donor_b = (0xb22d_f002, 7);
        let mut buckets = vec![None; 128];
        buckets[11] = Some(donor_a);
        buckets[97] = Some(donor_b);
        let original_nonempty: Vec<bool> = buckets.iter().map(Option::is_some).collect();

        for j in 0..buckets.len() {
            if !original_nonempty[j] {
                let source = ers.densification_source(&original_nonempty, j);
                assert!(
                    original_nonempty[source],
                    "bucket {j} selected source {source}"
                );
            }
        }

        ers.densify(&mut buckets);

        assert!(buckets
            .iter()
            .all(|bucket| matches!(bucket, Some(v) if *v == donor_a || *v == donor_b)));
    }

    #[test]
    fn ers_tracks_exact_weighted_jaccard_across_hash_seeds() {
        use crate::similarity::jaccard_similarity;

        let x = vec![(0, 0.5)];
        let y = vec![(0, 0.25), (1, 0.25)];
        let k = 4096usize;
        let l = 16u64;
        let seeds = 32u64;
        let expected = jaccard_similarity(&x, &y);
        let standard_error = (expected * (1.0 - expected) / (k as f64 * seeds as f64)).sqrt();
        let tolerance = 6.0 * standard_error + 0.006;

        // Loose caps leave roughly 1% of buckets empty at L=16, exercising
        // densification and verifying that correctness depends on drawing from
        // the full M=sum(m_i), not on every cap being tight. Average independent
        // hash seeds because simple tabulation does not concentrate on every
        // fixed low-dimensional input.
        for (label, caps) in [("tight", vec![0.5, 0.25]), ("loose", vec![1.0, 1.0])] {
            let mut hits = 0usize;
            for seed in 0..seeds {
                let mut hash_rng = mt_from_seed(0xe255_2004 + seed);
                let ers = ErsWmh::new_mt(&mut hash_rng, &caps, k as u64);
                let sk_x = ers.sketch(&x, Some(l));
                let sk_y = ers.sketch(&y, Some(l));
                hits += sk_x.iter().zip(&sk_y).filter(|(a, b)| a.0 == b.0).count();
            }
            let observed = hits as f64 / (k as f64 * seeds as f64);

            assert!(
                (observed - expected).abs() <= tolerance,
                "{label} caps: exact Jaccard={expected:.6}, observed={observed:.6}, tolerance={tolerance:.6}"
            );
        }
    }

    #[test]
    fn ers_early_stop_fills_all_buckets() {
        let mut data_rng = mt_from_seed(1337);
        let mut hash_rng = mt_from_seed(0xe255_0001);
        let d = 200_000usize;
        let k = 4096;

        // Base set
        let x = generate_weighted_set(d, 50_000, 10_000.0, &mut data_rng);

        // Build tight caps **from the data actually being sketched**
        let m = caps_from_sets(d, &[&x]);

        // ERS with data-consistent caps
        let ers = ErsWmh::new_mt(&mut hash_rng, &m, k as u64);

        // Algorithm 2: per-hash sequence length L
        let l: u64 = 512;
        let sk = ers.sketch(&x, Some(l));

        assert_eq!(sk.len(), k);
    }

    #[test]
    fn ers_approximates_weighted_jaccard() {
        use crate::similarity::jaccard_similarity;

        let mut data_rng = mt_from_seed(8675309);
        let d = 200_000usize;
        let k = 4096;

        // Fixed per-hash sequence length (Algorithm 2)
        let l: u64 = 1024; // with tight caps, even smaller L often works

        let l0 = 50_000u64;
        let l1 = 10_000.0;
        let x = generate_weighted_set(d, l0, l1, &mut data_rng);
        let targets = [
            0.99, 0.96, 0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1,
            0.05, 0.01,
        ];

        for (target_idx, &rel) in targets.iter().enumerate() {
            let y = generate_similar_weighted_set(d, &x, rel, &mut data_rng);
            let j_true = jaccard_similarity(&x, &y);
            println!("true weighted Jaccard: {:?}", j_true);

            // Tight caps must dominate BOTH vectors in the comparison
            let m_per_dim = caps_from_sets(d, &[&x, &y]);

            // Rebuild ERS for this pair with valid caps
            let mut hash_rng = mt_from_seed(0xe255_1000 ^ (target_idx as u64));
            let ers = ErsWmh::new_mt(&mut hash_rng, &m_per_dim, k as u64);

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
