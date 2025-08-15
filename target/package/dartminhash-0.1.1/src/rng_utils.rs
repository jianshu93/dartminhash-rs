//! Thin wrapper to build an MT19937 rng that implements `RngCore`.

use mt19937::{MT19937, Seed};
use rand_core::SeedableRng;

// Convenience function: build MT19937 from a single u64 seed.
pub fn mt_from_seed(seed64: u64) -> MT19937 {
    // Expand u64 into 624 u32s deterministically.
    // Simpler: use MT19937::new_with_slice_seed(&[seed_low, seed_high, ...]).
    let mut seed_arr = [0u32; mt19937::N];
    seed_arr[0] = (seed64 & 0xFFFF_FFFF) as u32;
    seed_arr[1] = (seed64 >> 32) as u32;
    let seed = Seed(seed_arr);
    MT19937::from_seed(seed)
}

// Allow MT19937 to be used as a mutable RNGCore directly
pub type MtRng = MT19937;