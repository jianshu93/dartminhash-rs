use tab_hash::{Tab32Simple, Tab64Simple};
use rand_core::RngCore;

// Convert a u64 to a uniform double in [0,1)
#[inline]
pub fn to_unit(x: u64) -> f64 {
    // divide by 2^64 - 1
    (x as f64) / 18446744073709551615.0
}

// Split a u64 into two 32-bit halves → two uniforms in [0,1)
#[inline]
pub fn to_units(x: u64) -> (f64, f64) {
    let hi = (x >> 32) as u32;
    let lo = (x & 0xFFFF_FFFF) as u32;
    (
        (hi as f64) / 4294967295.0,
        (lo as f64) / 4294967295.0,
    )
}

// Simple wrapper to create a Tab32Simple table seeded by an RNG
pub fn tab32_from_rng<R: RngCore>(rng: &mut R) -> Tab32Simple {
    let mut table = vec![vec![0u32; 256]; 4];
    for i in 0..4 {
        for j in 0..256 {
            table[i][j] = rng.next_u32();
        }
    }
    Tab32Simple::from_vec(table)
}

// Wrapper to create a Tab64Simple table seeded by an RNG
pub fn tab64_from_rng<R: RngCore>(rng: &mut R) -> Tab64Simple {
    let mut table = vec![vec![0u64; 256]; 8];
    for i in 0..8 {
        for j in 0..256 {
            table[i][j] = rng.next_u64();
        }
    }
    Tab64Simple::from_vec(table)
}

// Weighted sum helper

pub fn total_weight(x: &[(u64, f64)]) -> f64 {
    x.iter().map(|(_, w)| w).sum()
}