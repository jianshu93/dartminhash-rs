<div align="center">
  <img width="50%" src ="DartMinHash_logo.png">
</div>

# DartMinHash: Fast Sketching for Weighted Sets
This crate provides the implementation of [DartMinHash](https://arxiv.org/abs/2005.11547) (1) algorithm for estimation of weighted Jaccard similarity. To reproduce the algorithm in the paper, we use the same tabulation hashing idea (2). Mersenne Twister PRNG was used as seed.  Other high quality 64-bit hash functions such as xxhash-rust or whyash-rs should also work as well. 

# Install & test
Add to your Cargo.toml. Official release in crates.io will come soon.
```bash
dartminhash = { git = "https://github.com/jianshu93/dartminhash-rs.git" }
```

Test case to evaulate the accuracy of the algorithm
```bash
cargo test --release dartminhash_approximates_weighted_jaccard -- --nocapture
```
Test for true weighted Jaccard from 0.005 to 0.98. Output:

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

Otmar Ertl. 2023. TreeMinHash: https://github.com/oertl/treeminhash
