//! dartminhash_rs: A Rust implementation of DartMinHash (weighted MinHash)
//! inspired by https://github.com/tobc/dartminhash (C++).
//!
//! Main items:
//! - [`dart_hash::DartHash`] : produces darts from a weighted feature vector
//! - [`dart_minhash::DartMinHash`] : turns darts into a k-sized MinHash sketch
//!
//! Feature universe element = `(u64 id, f64 weight)`

pub mod hash_utils;
pub mod rng_utils;
pub mod darthash;
pub mod dartminhash;
pub mod similarity;
pub mod rejsmp;

pub use crate::darthash::DartHash;
pub use crate::dartminhash::DartMinHash;
pub use crate::similarity::{
    weight,
    intersection,
    jaccard_similarity,
    l1_similarity,
    hamming_distance,
    onebit_minhash_jaccard_estimate,
    jaccard_from_l1,
    l1_from_jaccard,
    count_collisions,
    jaccard_estimate_from_minhashes
};
pub use crate::rejsmp::ErsWmh;