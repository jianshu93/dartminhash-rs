<div align="center">
  <img width="50%" src ="DartMinHash_logo.png">
</div>

# DartMinHash: Fast Sketching for Weighted Sets
This crate provides the implementation of [DartMinHash](https://arxiv.org/abs/2005.11547) (1) algorithm for estimation of weighted Jaccard similarity. To reproduce the algorithm in the paper, we use the same tabulation hashing idea (2). Other high quality 64-bit hash functions such as xxhash-rust or whyash-rs should also work as well. 

# Install & test
Add to your Cargo.toml. Official release in crates.io will come soon.
```bash
dartminhash = { git = "https://github.com/jianshu93/dartminhash-rs.git" }
```

Test case to evaulate the accuracy of the algorithm
```bash
cargo test --release dartminhash_approximates_weighted_jaccard -- --nocapture
```
Test for true weighted Jaccard from 0.005 to 0.8. Output:

```bash
true weighted Jaccard: 0.81818181818183
estimated weighted Jaccard: 0.8203125
true weighted Jaccard: 0.5384615384615511
estimated weighted Jaccard: 0.53515625
true weighted Jaccard: 0.3333333333333311
estimated weighted Jaccard: 0.334228515625
true weighted Jaccard: 0.17647058823529124
estimated weighted Jaccard: 0.1806640625
true weighted Jaccard: 0.05263157894736854
estimated weighted Jaccard: 0.052734375
true weighted Jaccard: 0.02564102564102568
estimated weighted Jaccard: 0.025390625
true weighted Jaccard: 0.005025125628140737
estimated weighted Jaccard: 0.0048828125


```

# Usage

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


# References
Christiani, T., 2020. Dartminhash: Fast sketching for weighted sets. arXiv preprint arXiv:2005.11547.

Pǎtraşcu, M. and Thorup, M., 2012. The power of simple tabulation hashing. Journal of the ACM (JACM), 59(3), pp.1-50.
