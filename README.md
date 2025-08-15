[![Latest Version](https://img.shields.io/crates/v/dartminhash?style=for-the-badge&color=mediumpurple&logo=rust)](https://crates.io/crates/dartminhash)
[![docs.rs](https://img.shields.io/docsrs/dartminhash?style=for-the-badge&logo=docs.rs&color=mediumseagreen)](https://docs.rs/dartminhash/latest/dartminhash/)



<div align="center">
  <img width="50%" src ="DartMinHash_logo.png">
</div>

# DartMinHash & Rejection Sampling: Fast Sketching for Weighted Sets
This crate provides the implementation of [DartMinHash](https://arxiv.org/abs/2005.11547) (1), [Rejection Sampling](https://proceedings.neurips.cc/paper/2016/hash/c2626d850c80ea07e7511bbae4c76f4b-Abstract.html) (2) and [Efficient Rejection Sampling](https://ojs.aaai.org/index.php/AAAI/article/view/16543) (3) algorithm for estimation of weighted Jaccard similarity. To reproduce the algorithm in the paper, we use the same tabulation hashing idea (4). Mersenne Twister PRNG was used as seed.  Other high quality 64-bit hash functions such as xxhash-rust or whyash-rs should also work as well. 

Note: DartMinHash is only significantly faster than (Efficient) Rejection Sampling (2,3) for very sparse vectors, that is the number of nonezero elements (d) is less than ~2% of vector dimension (D) on average for all vectors. This is especially true for large-scale datasets. However, For RS and ERS, a maxmimum value of weight for input vector must be known. Otherwise, the estimation is significantly biased (6). Therefore, general applicability is limited by the required priori knowledge of sharp upper bounds for $w_{max}(d)$. Also, ERS is not unbiased (3).

# Install & test
Add below lines to your Cargo.toml dependencies. Official release in crates.io is [here](https://crates.io/crates/dartminhash).
```bash
dartminhash = "0.1.0"
```

Test case to evaulate the accuracy of the DartMinHash algorithm.
```bash
cargo test --release dartminhash_approximates_weighted_jaccard -- --nocapture
```
Test for true weighted Jaccard from 0.005 to 0.98. Output for DartMinHash:

```bash
true weighted Jaccard: 0.9801980198019941
estimated weighted Jaccard: 0.98046875
true weighted Jaccard: 0.9230769230769282
estimated weighted Jaccard: 0.921142578125
true weighted Jaccard: 0.8691588785046781
estimated weighted Jaccard: 0.872802734375
true weighted Jaccard: 0.8181818181818195
estimated weighted Jaccard: 0.8193359375
true weighted Jaccard: 0.7391304347825998
estimated weighted Jaccard: 0.73876953125
true weighted Jaccard: 0.6666666666666655
estimated weighted Jaccard: 0.67138671875
true weighted Jaccard: 0.5999999999999979
estimated weighted Jaccard: 0.6025390625
true weighted Jaccard: 0.5384615384615447
estimated weighted Jaccard: 0.53662109375
true weighted Jaccard: 0.48148148148148356
estimated weighted Jaccard: 0.484619140625
true weighted Jaccard: 0.42857142857142844
estimated weighted Jaccard: 0.4267578125
true weighted Jaccard: 0.37931034482758635
estimated weighted Jaccard: 0.37060546875
true weighted Jaccard: 0.3333333333333311
estimated weighted Jaccard: 0.336181640625
true weighted Jaccard: 0.25000000000000155
estimated weighted Jaccard: 0.245361328125
true weighted Jaccard: 0.17647058823529352
estimated weighted Jaccard: 0.17626953125
true weighted Jaccard: 0.11111111111111141
estimated weighted Jaccard: 0.109130859375
true weighted Jaccard: 0.05263157894736854
estimated weighted Jaccard: 0.05224609375
true weighted Jaccard: 0.02564102564102588
estimated weighted Jaccard: 0.026123046875
true weighted Jaccard: 0.005025125628140735
estimated weighted Jaccard: 0.004638671875

```

Test case to evaulate the accuracy of RS and ERS algorithms.
```bash
cargo test --release  rs_approximates_weighted_jaccard_only  -- --nocapture
cargo test --release ers_approximates_weighted_jaccard -- --nocapture
```

Test for true weighted Jaccard from 0.005 to 0.98. Output for RS:

```bash
true weighted Jaccard: 0.9801980198019689
estimated weighted Jaccard: 0.979736328125
true weighted Jaccard: 0.9230769230769127
estimated weighted Jaccard: 0.921875
true weighted Jaccard: 0.8691588785046708
estimated weighted Jaccard: 0.866455078125
true weighted Jaccard: 0.8181818181817918
estimated weighted Jaccard: 0.82373046875
true weighted Jaccard: 0.7391304347825958
estimated weighted Jaccard: 0.740966796875
true weighted Jaccard: 0.6666666666666508
estimated weighted Jaccard: 0.67041015625
true weighted Jaccard: 0.5999999999999872
estimated weighted Jaccard: 0.6005859375
true weighted Jaccard: 0.5384615384615266
estimated weighted Jaccard: 0.530029296875
true weighted Jaccard: 0.48148148148147907
estimated weighted Jaccard: 0.482421875
true weighted Jaccard: 0.42857142857142794
estimated weighted Jaccard: 0.435546875
true weighted Jaccard: 0.37931034482758264
estimated weighted Jaccard: 0.38818359375
true weighted Jaccard: 0.3333333333333351
estimated weighted Jaccard: 0.32373046875
true weighted Jaccard: 0.24999999999999478
estimated weighted Jaccard: 0.243896484375
true weighted Jaccard: 0.17647058823529485
estimated weighted Jaccard: 0.1796875
true weighted Jaccard: 0.11111111111111124
estimated weighted Jaccard: 0.112548828125
true weighted Jaccard: 0.052631578947367655
estimated weighted Jaccard: 0.04833984375
true weighted Jaccard: 0.025641025641025373
estimated weighted Jaccard: 0.025390625
true weighted Jaccard: 0.005025125628140681
estimated weighted Jaccard: 0.004638671875
```


Test for true weighted Jaccard from 0.005 to 0.98. Output for ERS:
```bash
true weighted Jaccard: 0.9801980198019894
estimated weighted Jaccard: 0.979248046875
true weighted Jaccard: 0.923076923076927
estimated weighted Jaccard: 0.916259765625
true weighted Jaccard: 0.8691588785046818
estimated weighted Jaccard: 0.870361328125
true weighted Jaccard: 0.8181818181818254
estimated weighted Jaccard: 0.815185546875
true weighted Jaccard: 0.7391304347826214
estimated weighted Jaccard: 0.728271484375
true weighted Jaccard: 0.6666666666666645
estimated weighted Jaccard: 0.661865234375
true weighted Jaccard: 0.6000000000000002
estimated weighted Jaccard: 0.602783203125
true weighted Jaccard: 0.538461538461545
estimated weighted Jaccard: 0.54296875
true weighted Jaccard: 0.48148148148149517
estimated weighted Jaccard: 0.4736328125
true weighted Jaccard: 0.4285714285714275
estimated weighted Jaccard: 0.42724609375
true weighted Jaccard: 0.3793103448275884
estimated weighted Jaccard: 0.373291015625
true weighted Jaccard: 0.3333333333333343
estimated weighted Jaccard: 0.320068359375
true weighted Jaccard: 0.2500000000000003
estimated weighted Jaccard: 0.241943359375
true weighted Jaccard: 0.17647058823529319
estimated weighted Jaccard: 0.187255859375
true weighted Jaccard: 0.11111111111111192
estimated weighted Jaccard: 0.110595703125
true weighted Jaccard: 0.05263157894736851
estimated weighted Jaccard: 0.04638671875
true weighted Jaccard: 0.025641025641025394
estimated weighted Jaccard: 0.026611328125
true weighted Jaccard: 0.00502512562814068
estimated weighted Jaccard: 0.00537109375

```

# Usage
DartMinhash: 
```rust
use dartminhash::dartminhash::DartMinHash;
use dartminhash::rng_utils::mt_from_seed;
use dartminhash::similarity::jaccard_estimate_from_minhashes;

fn main() {
    let mut rng = mt_from_seed(42);
    let k = 128;

    let dartminhash = DartMinHash::new_mt(&mut rng, k);

    // Weighted inputs: overlap in IDs, but weights differ a bit
    let sample_a = vec![
        (5, 1.2),
        (17, 0.9),
        (23, 1.1),
        (42, 0.95),
        (100, 1.0),
    ];
    let sample_b = vec![
        (5, 1.0),
        (17, 1.0),
        (44, 1.1),
        (100, 1.05),
    ];

    let sketch_a = dartminhash.sketch(&sample_a);
    let sketch_b = dartminhash.sketch(&sample_b);

    let est_jaccard = jaccard_estimate_from_minhashes(&sketch_a, &sketch_b);

    println!("Estimated weighted Jaccard similarity: {:.4}", est_jaccard);
}

```

Rejection Sampling:

```rust
use dartminhash::{RsWmh};
use dartminhash::rng_utils::mt_from_seed;

/// Build caps that dominate the given vectors: m_i = max(1, ceil(max(x_i,y_i))).
fn caps_from_pair(d: usize, a: &[(u64, f64)], b: &[(u64, f64)]) -> Vec<u32> {
    let mut m = vec![1u32; d];
    for &(i, w) in a.iter().chain(b.iter()) {
        if w > 0.0 {
            let cap = (w.ceil() as u32).max(1);
            let idx = i as usize;
            if cap > m[idx] { m[idx] = cap; }
        }
    }
    m
}

fn main() {
    // Suppose your ids are < D
    let d: usize = 1_000;       // feature universe size (just an example)
    let k: usize = 128;         // number of hashes
    let mut rng = mt_from_seed(42);

    // Weighted inputs
    let sample_a = vec![
        (5, 1.2),
        (17, 0.9),
        (23, 1.1),
        (42, 0.95),
        (100, 1.0),
    ];
    let sample_b = vec![
        (5, 1.0),
        (17, 1.0),
        (44, 1.1),
        (100, 1.05),
    ];

    // Build per-dimension caps (for production, precompute once across your dataset)
    let m_per_dim = caps_from_pair(d, &sample_a, &sample_b);

    // RS: k independent hashes → k 64-bit ids (collisions estimate J)
    let rs = RsWmh::new_mt(&mut rng, &m_per_dim, k);
    let sig_a = rs.sketch_ids(&sample_a);
    let sig_b = rs.sketch_ids(&sample_b);

    // Collision-rate estimator
    let hits = sig_a.iter().zip(sig_b.iter()).filter(|(x, y)| x == y).count();
    let est_jaccard = hits as f64 / (k as f64);

    println!("RS estimated weighted Jaccard: {:.4}", est_jaccard);

    // (Optional) You can also inspect geometric trial counts:
    // let counts_a = rs.sketch_counts(&sample_a);
    // println!("example RS counts head: {:?}", &counts_a[..8]);
}

```

Efficient Rejection Sampling:

```rust
use dartminhash::{ErsWmh};
use dartminhash::rng_utils::mt_from_seed;

fn caps_from_sets(d: usize, sets: &[&[(u64, f64)]]) -> Vec<u32> {
    let mut m = vec![1u32; d];
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

fn main() {
    let mut rng = mt_from_seed(1337);

    let d: usize = 200_000;
    let k: u64   = 1024;
    let L: u64   = 512;   // try 256–1024; larger L → fewer pre-densify empties

    // Two weighted vectors
    let a = vec![(5, 1.2), (17, 0.9), (23, 1.1), (42, 0.95), (100, 1.0)];
    let b = vec![(5, 1.0), (17, 1.0), (44, 1.1), (100, 1.05)];

    // Caps must dominate both vectors
    let m_per_dim = caps_from_sets(d, &[&a, &b]);

    let ers = ErsWmh::new_mt(&mut rng, &m_per_dim, k);

    // ERS returns k (id, rank) pairs; collisions on id estimate Jaccard
    let sk_a = ers.sketch(&a, Some(L));
    let sk_b = ers.sketch(&b, Some(L));

    let hits = sk_a.iter().zip(&sk_b).filter(|(x, y)| x.0 == y.0).count();
    let j_est = hits as f64 / k as f64;

    println!("ERS (L={}) estimated weighted Jaccard: {:.4}", L, j_est);
}

```



# References
1.Christiani, T., 2020. Dartminhash: Fast sketching for weighted sets. arXiv preprint arXiv:2005.11547.

2.Shrivastava, A., 2016. Simple and efficient weighted minwise hashing. Advances in Neural Information Processing Systems, 29.

3.Li, X. and Li, P., 2021, May. Rejection sampling for weighted jaccard similarity revisited. In Proceedings of the AAAI Conference on Artificial Intelligence (Vol. 35, No. 5, pp. 4197-4205).

4.Pǎtraşcu, M. and Thorup, M., 2012. The power of simple tabulation hashing. Journal of the ACM (JACM), 59(3), pp.1-50.

5.Ertl, O. (2025) “TreeMinHash: Fast Sketching for Weighted Jaccard Similarity Estimation”. Zenodo. doi: 10.5281/zenodo.16730965.

6.Ertl, O., 2018, July. Bagminhash-minwise hashing algorithm for weighted sets. In Proceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining (pp. 1368-1377).
