//! dartminhash_rs: A Rust implementation of weighted MinHash, including DartMinHash, TreeMinHash and Efficient Rejection Sampling
//! inspired by https://github.com/tobc/dartminhash (C++).
//!
//! Main items:
//! - [`dart_hash::DartHash`] : produces darts from a weighted feature vector
//! - [`dart_minhash::DartMinHash`] : turns darts into a k-sized MinHash sketch
//!
//! Feature universe element = `(u64 id, f64 weight)`

pub mod darthash;
pub mod dartminhash;
pub mod hash_utils;
pub mod rejsmp;
pub mod rng_utils;
pub mod similarity;

pub use crate::darthash::DartHash;
pub use crate::dartminhash::DartMinHash;
pub use crate::rejsmp::ErsWmh;
pub use crate::similarity::{
    count_collisions, hamming_distance, intersection, jaccard_estimate_from_minhashes,
    jaccard_from_l1, jaccard_similarity, l1_from_jaccard, l1_similarity,
    onebit_minhash_jaccard_estimate, weight,
};

pub mod treeminhash;

pub use crate::treeminhash::TreeMinHash;
