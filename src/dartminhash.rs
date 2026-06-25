//! DartMinHash: DartHash + repeatedly throws darts until all buckets filled.

use crate::darthash::{Dart, DartHash};
use crate::hash_utils::*;
use crate::rng_utils::MtRng;
use std::f64::INFINITY;

#[cfg(feature = "mixed_tab")]
type Tab64Bucket = tab_hash::Tab64Mixed;
#[cfg(not(feature = "mixed_tab"))]
type Tab64Bucket = tab_hash::Tab64Simple;

#[cfg(feature = "mixed_tab")]
fn tab64_bucket_from_rng(rng: &mut MtRng) -> Tab64Bucket {
    mixed_tab64_from_rng(rng)
}

#[cfg(not(feature = "mixed_tab"))]
fn tab64_bucket_from_rng(rng: &mut MtRng) -> Tab64Bucket {
    tab64_from_rng(rng)
}

// Sketch = k slots of (id, rank)
pub type MinHashSketch = Vec<Dart>;

pub struct DartMinHash {
    k: u64,
    bucket_hasher: Tab64Bucket,
    dart_hash: DartHash,
}

impl DartMinHash {
    // t = k*ln(k) + 2k
    pub fn new_mt(rng: &mut MtRng, k: u64) -> Self {
        let t = ((k as f64) * (k as f64).ln() + 2.0 * (k as f64)).ceil() as u64;
        let bucket_hasher = tab64_bucket_from_rng(rng);
        let dart_hash = DartHash::new_mt(rng, t);
        Self {
            k,
            bucket_hasher,
            dart_hash,
        }
    }

    // Returns k minhash darts. Ensures every bucket got something by increasing theta if needed.
    pub fn sketch(&self, x: &[(u64, f64)]) -> MinHashSketch {
        let mut minhashes = vec![(0u64, INFINITY); self.k as usize];
        let mut theta = 1.0;
        loop {
            let mut filled = vec![false; self.k as usize];
            let darts = self.dart_hash.darts(x, theta);
            for &(id, rank) in &darts {
                let j = (self.bucket_hasher.hash(id) % self.k) as usize;
                filled[j] = true;
                if rank < minhashes[j].1 {
                    minhashes[j] = (id, rank);
                }
            }
            if filled.iter().all(|&b| b) {
                break;
            }
            theta += 0.5;
        }
        minhashes
    }
}

#[cfg(test)]
mod tests {
    use rand_core::RngCore;
    use std::collections::HashSet;

    use crate::{
        dartminhash::DartMinHash,
        rng_utils::{MtRng, mt_from_seed},
        similarity::{jaccard_estimate_from_minhashes, jaccard_similarity},
        treeminhash::TreeMinHash,
    };
    /// Uniform(0,1) using the same MT19937 rng.
    fn uniform01(rng: &mut MtRng) -> f64 {
        mt19937::gen_res53(rng)
    }

    /// Generate a random weighted set:
    /// Pick L0 distinct random indices (u64)
    ///  Draw L0-1 uniform(0,1), sort, use the gaps * L1 as weights
    /// Returns sorted by id.
    fn generate_weighted_set(l0: u64, l1: f64, rng: &mut MtRng) -> Vec<(u64, f64)> {
        let mut elements = HashSet::with_capacity(l0 as usize);
        while elements.len() < l0 as usize {
            elements.insert(rng.next_u64());
        }

        // Uniform splitters
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
            out.push((idx, w));
            prev = z[j];
            j += 1;
        }
        out.sort_by_key(|p| p.0);
        out
    }

    /// Generate Y from X with a target relative overlap:
    /// y = relative_overlap * x  (element-wise scaling)
    /// plus the remaining mass as a new element not in x.
    fn generate_similar_weighted_set(
        x: &[(u64, f64)],
        relative_overlap: f64,
        rng: &mut MtRng,
    ) -> Vec<(u64, f64)> {
        // Pick a free id not in x
        let mut free_id;
        'find_id: loop {
            free_id = rng.next_u64();
            if x.binary_search_by_key(&free_id, |p| p.0).is_err() {
                break 'find_id;
            }
        }

        let mut excess = 0.0;
        let mut y = Vec::with_capacity(x.len() + 1);
        for &(id, w) in x {
            let w_scaled = w * relative_overlap;
            excess += w - w_scaled;
            y.push((id, w_scaled));
        }
        y.push((free_id, excess));
        y.sort_by_key(|p| p.0);
        y
    }

    #[test]
    fn dartminhash_approximates_weighted_jaccard() {
        let mut data_rng = mt_from_seed(1337);
        let mut hash_rng = mt_from_seed(0xd417_0001);

        // data size and mass
        let l0 = 50_000; // number of nonzeros
        let l1 = 10_000.0; // total weight (approximately)
        let k = 4096; // sketch size

        let dm = DartMinHash::new_mt(&mut hash_rng, k);

        // Generate a base set
        let x = generate_weighted_set(l0, l1, &mut data_rng);
        assert_eq!(x.len(), l0 as usize);

        // Try a few overlaps
        let targets = [
            0.99, 0.96, 0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1,
            0.05, 0.01,
        ];
        // tolerance on Jaccard estimate (tweak with k)

        for &rel in &targets {
            let y = generate_similar_weighted_set(&x, rel, &mut data_rng);

            // True weighted Jaccard
            let j_true = jaccard_similarity(&x, &y);
            println!("true weighted Jaccard: {:?}", j_true);
            // Sketch
            let sk_x = dm.sketch(&x);
            let sk_y = dm.sketch(&y);
            assert_eq!(sk_x.len(), k as usize);
            assert_eq!(sk_y.len(), k as usize);

            let j_est = jaccard_estimate_from_minhashes(&sk_x, &sk_y);
            println!("estimated weighted Jaccard: {:?}", j_est);
            // Check accuracy
            let sd = (j_true * (1.0 - j_true) / (k as f64)).sqrt();
            let tol = (3.2 * sd).max(1.25 / (k as f64).sqrt());
            let err = (j_true - j_est).abs();
            assert!(
                err <= tol,
                "rel_overlap={rel}, true={j_true:.4}, est={j_est:.4}, err={err:.4} > tol={tol}"
            );
        }
    }
    #[test]
    fn dartminhash2_approximates_weighted_jaccard_sparse() {
        let mut data_rng = mt_from_seed(2025);
        let mut hash_rng = mt_from_seed(0xd417_0002);

        // "5% sparse": pick a much smaller l0 than the (implicit) universe size.
        // (Your helpers don't take D; sparsity here is "few nonzeros" relative to a huge ID space.)
        let l0 = 5_000; // number of nonzeros (~5% of a conceptual D=1,000,000)
        let l1 = 3_000.0; // total weight (kept moderate)
        let k = 4096; // sketch size (smaller than 4096 so this test is fast)

        let dm = DartMinHash::new_mt(&mut hash_rng, k);

        // Base set
        let x = generate_weighted_set(l0, l1, &mut data_rng);
        assert_eq!(x.len(), l0 as usize);

        // A range of true Jaccards, as in your other tests
        let targets = [
            0.99, 0.96, 0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1,
            0.05, 0.01,
        ];

        for &rel in &targets {
            let y = generate_similar_weighted_set(&x, rel, &mut data_rng);

            // Ground truth
            let j_true = jaccard_similarity(&x, &y);
            println!("true weighted Jaccard (DMH sparse5): {:?}", j_true);

            // Sketch & estimate
            let sk_x = dm.sketch(&x);
            let sk_y = dm.sketch(&y);
            assert_eq!(sk_x.len(), k as usize);
            assert_eq!(sk_y.len(), k as usize);

            let j_est = jaccard_estimate_from_minhashes(&sk_x, &sk_y);
            println!("estimated weighted Jaccard (DMH sparse5): {:?}", j_est);

            // σ-aware tolerance (matches the RS/ERS tests you have)
            let sd = (j_true * (1.0 - j_true) / (k as f64)).sqrt();
            let tol = (3.2 * sd).max(1.25 / (k as f64).sqrt());

            let err = (j_true - j_est).abs();
            assert!(
                err <= tol,
                "DMH sparse5: true={j_true:.6}, est={j_est:.6}, err={err:.6}, tol={tol:.6}"
            );
        }
    }

    /// Empirical DartMinHash accuracy sweep across several independent data and
    /// hasher seeds. Run twice to compare hash families:
    ///
    ///     cargo test --release --no-default-features dartminhash_multi_seed_accuracy_sweep -- --ignored --nocapture
    ///     cargo test --release --features mixed_tab dartminhash_multi_seed_accuracy_sweep -- --ignored --nocapture
    #[test]
    #[ignore]
    fn dartminhash_multi_seed_accuracy_sweep() {
        let seeds = [7, 19, 42, 1_337, 2_025, 86_753_09];
        let targets = [
            0.99, 0.96, 0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1,
            0.05, 0.01,
        ];

        let l0 = 10_000;
        let l1 = 10_000.0;
        let k = 2048;

        let mut cases = 0usize;
        let mut sum_abs = 0.0;
        let mut sum_sq = 0.0;
        let mut max_abs = 0.0;
        let mut worst_seed = 0;
        let mut worst_rel = 0.0;
        let mut worst_true = 0.0;
        let mut worst_est = 0.0;

        for &seed in &seeds {
            let mut data_rng = mt_from_seed(seed);
            let x = generate_weighted_set(l0, l1, &mut data_rng);

            // Keep the input data stream independent from hasher construction so
            // simple and mixed builds compare on exactly the same weighted sets.
            let mut hash_rng = mt_from_seed(seed ^ 0x9e37_79b9_7f4a_7c15);
            let dm = DartMinHash::new_mt(&mut hash_rng, k);
            let sk_x = dm.sketch(&x);

            for &rel in &targets {
                let y = generate_similar_weighted_set(&x, rel, &mut data_rng);
                let j_true = jaccard_similarity(&x, &y);
                let sk_y = dm.sketch(&y);
                let j_est = jaccard_estimate_from_minhashes(&sk_x, &sk_y);
                let err = (j_true - j_est).abs();

                cases += 1;
                sum_abs += err;
                sum_sq += err * err;
                if err > max_abs {
                    max_abs = err;
                    worst_seed = seed;
                    worst_rel = rel;
                    worst_true = j_true;
                    worst_est = j_est;
                }
            }
        }

        let mean_abs = sum_abs / cases as f64;
        let rmse = (sum_sq / cases as f64).sqrt();
        let mode = if cfg!(feature = "mixed_tab") {
            "mixed_tab"
        } else {
            "simple_tab"
        };

        println!(
            "DMH_MULTI_SEED mode={mode} seeds={} targets={} cases={cases} k={k} l0={l0} mean_abs={mean_abs:.8} rmse={rmse:.8} max_abs={max_abs:.8} worst_seed={worst_seed} worst_rel={worst_rel:.3} worst_true={worst_true:.8} worst_est={worst_est:.8}",
            seeds.len(),
            targets.len()
        );

        assert!(mean_abs.is_finite());
        assert!(rmse.is_finite());
        assert!(max_abs.is_finite());
    }

    /// Large raw-count / absolute-weight simulation.
    ///
    /// This is intentionally ignored because it is a timing-style stress test,
    /// not a small unit test. It is meant to mimic large-weight scenarios,
    /// where the sum of set masses are very large,
    /// which is the case where TreeMinHash should be tested against DartMinHash.
    ///
    /// Run with:
    ///
    ///     cargo test dartminhash_large_weight_sum_vs_treeminhash --release -- --ignored --nocapture
    #[test]
    #[ignore]
    fn dartminhash_large_weight_sum_vs_treeminhash() {
        use std::hint::black_box;
        use std::time::Instant;

        let mut data_rng = mt_from_seed(98_001);

        // Slightly larger sparse set than the ordinary correctness tests.
        // l1 is deliberately huge to simulate raw-count absolute weighted data.
        let l0 = 120_000u64;
        let l1 = 1.0e12f64;
        let k = 4096u64;

        let x = generate_weighted_set(l0, l1, &mut data_rng);
        assert_eq!(x.len(), l0 as usize);

        let actual_sum: f64 = x.iter().map(|&(_, w)| w).sum();
        println!(
            "large raw-count simulation: nonzeros={}, requested_sum={:.3e}, actual_sum={:.3e}, k={}",
            l0, l1, actual_sum, k
        );

        let mut dmh_rng = mt_from_seed(12_345);
        let dmh = DartMinHash::new_mt(&mut dmh_rng, k);
        let dmh_start = Instant::now();
        let sk_dmh = black_box(dmh.sketch(black_box(&x)));
        let dmh_elapsed = dmh_start.elapsed();
        assert_eq!(sk_dmh.len(), k as usize);
        assert!(sk_dmh.iter().all(|&(_, r)| r.is_finite()));

        let mut tmh_rng = mt_from_seed(12_345);
        let tmh = TreeMinHash::new_mt(&mut tmh_rng, k);
        let tmh_start = Instant::now();
        let sk_tmh = black_box(tmh.sketch(black_box(&x)));
        let tmh_elapsed = tmh_start.elapsed();
        assert_eq!(sk_tmh.len(), k as usize);
        assert!(sk_tmh.iter().all(|&(_, r)| r.is_finite()));

        println!("DartMinHash elapsed: {:?}", dmh_elapsed);
        println!("TreeMinHash elapsed: {:?}", tmh_elapsed);
        if tmh_elapsed.as_nanos() > 0 {
            println!(
                "DMH / TMH elapsed ratio: {:.3}x",
                dmh_elapsed.as_secs_f64() / tmh_elapsed.as_secs_f64()
            );
        }

        // Do not assert timing. Runtime ratios are machine/build dependent.
        // The purpose of this ignored test is to expose the large-weight regime
        // clearly under identical input size, total mass, sketch size, and seed.
    }

    #[test]
    fn conversions_match() {
        let x_w = 10.0;
        let y_w = 8.0;
        let j = 0.4;
        let l1 = crate::similarity::l1_from_jaccard(x_w, y_w, j);
        let j2 = crate::similarity::jaccard_from_l1(x_w, y_w, l1);
        assert!((j - j2).abs() < 1e-12);
    }
}
