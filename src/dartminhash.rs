//! DartMinHash: DartHash + repeatedly throws darts until all buckets filled.

use std::f64::INFINITY;
use crate::darthash::{Dart, DartHash};
use crate::hash_utils::tab64_from_rng;
use crate::rng_utils::MtRng;
use tab_hash::Tab64Simple;

// Sketch = k slots of (id, rank)
pub type MinHashSketch = Vec<Dart>;

pub struct DartMinHash {
    k: u64,
    bucket_hasher: Tab64Simple,
    dart_hash: DartHash,
}

impl DartMinHash {
    // t = k*ln(k) + 2k
    pub fn new_mt(rng: &mut MtRng, k: u64) -> Self {
        let t = ((k as f64) * (k as f64).ln() + 2.0 * (k as f64)).ceil() as u64;
        let bucket_hasher = tab64_from_rng(rng);
        let dart_hash = DartHash::new_mt(rng, t);
        Self { k, bucket_hasher, dart_hash }
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
            if filled.iter().all(|&b| b) { break; }
            theta += 0.5;
        }
        minhashes
    }

    // 1-bit sketch: take LSB of hash id of each bucket winner
    pub fn onebit(&self, x: &[(u64, f64)]) -> Vec<bool> {
        self.sketch(x).into_iter().map(|(id, _)| (id & 1) == 1).collect()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use rand_core::RngCore;

    use crate::{
        dartminhash::DartMinHash,
        rng_utils::{mt_from_seed, MtRng},
        similarity::{
            jaccard_estimate_from_minhashes, jaccard_similarity},
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
        for idx in elements {
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
        let mut rng = mt_from_seed(1337);

        // data size and mass
        let l0 = 50_000;      // number of nonzeros
        let l1 = 10_000.0;    // total weight (approximately)
        let k  = 4096;      // sketch size

        let dm = DartMinHash::new_mt(&mut rng, k);

        // Generate a base set
        let x = generate_weighted_set(l0, l1, &mut rng);
        assert_eq!(x.len(), l0 as usize);

        // Try a few overlaps
        let targets = [0.99, 0.96,0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01];
         // tolerance on Jaccard estimate (tweak with k)

        for &rel in &targets {
            let y = generate_similar_weighted_set(&x, rel, &mut rng);

            // True weighted Jaccard
            let j_true = jaccard_similarity(&x, &y);
            let tol = 1.0 / (k as f64).sqrt();
            println!("true weighted Jaccard: {:?}", j_true);
            // Sketch
            let sk_x = dm.sketch(&x);
            let sk_y = dm.sketch(&y);
            assert_eq!(sk_x.len(), k as usize);
            assert_eq!(sk_y.len(), k as usize);

            let j_est = jaccard_estimate_from_minhashes(&sk_x, &sk_y);
            println!("estimated weighted Jaccard: {:?}", j_est);
            // Check accuracy
            let err = (j_true - j_est).abs();
            assert!(
                err <= tol as f64,
                "rel_overlap={rel}, true={j_true:.4}, est={j_est:.4}, err={err:.4} > tol={tol}"
            );
        }
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
