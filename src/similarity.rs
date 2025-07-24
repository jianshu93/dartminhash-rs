//! src/similarity.rs
//! All vectors of `(id, weight)` are assumed sorted by id ascending.

use std::cmp::Ordering;

/// Sum of weights.
#[inline]
pub fn weight(x: &[(u64, f64)]) -> f64 {
    x.iter().map(|&(_, w)| w).sum()
}

// Intersection sum:  Σ min(w_x, w_y) over shared ids.
// Requires both x and y sorted by id.
pub fn intersection(x: &[(u64, f64)], y: &[(u64, f64)]) -> f64 {
    let (mut i, mut j, mut s) = (0usize, 0usize, 0.0);
    while i < x.len() && j < y.len() {
        match x[i].0.cmp(&y[j].0) {
            Ordering::Equal => {
                s += x[i].1.min(y[j].1);
                i += 1;
                j += 1;
            }
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
        }
    }
    s
}

// Jaccard similarity:  |x ∩ y| / |x ∪ y|  (weights).
#[inline]
pub fn jaccard_similarity(x: &[(u64, f64)], y: &[(u64, f64)]) -> f64 {
    let s = intersection(x, y);
    let wx = weight(x);
    let wy = weight(y);
    s / (wx + wy - s)
}

// L1 similarity (a.k.a. normalized intersection):  |x ∩ y| / min(|x|, |y|).
#[inline]
pub fn l1_similarity(x: &[(u64, f64)], y: &[(u64, f64)]) -> f64 {
    let s = intersection(x, y);
    let wx = weight(x);
    let wy = weight(y);
    s / wx.min(wy)
}

/// Hamming distance between two 1-bit sketches.
#[inline]
pub fn hamming_distance(x: &[bool], y: &[bool]) -> f64 {
    assert_eq!(x.len(), y.len(), "bit vectors must be same length");
    let mut h = 0.0;
    for (a, b) in x.iter().zip(y.iter()) {
        if a != b {
            h += 1.0;
        }
    }
    h
}

// One-bit MinHash Jaccard estimate (see original code):
// max(0, 2*(1 - H/T) - 1), where H is Hamming distance, T = length.
#[inline]
pub fn onebit_minhash_jaccard_estimate(x: &[bool], y: &[bool]) -> f64 {
    let h = hamming_distance(x, y);
    let t = x.len() as f64;
    (2.0 * (1.0 - h / t) - 1.0).max(0.0)
}

// Convert L1 similarity → Jaccard similarity.
#[inline]
pub fn jaccard_from_l1(x_weight: f64, y_weight: f64, l1_sim: f64) -> f64 {
    let inter = x_weight.min(y_weight) * l1_sim;
    let uni = x_weight + y_weight - inter;
    inter / uni
}

// Convert Jaccard similarity → L1 similarity.
#[inline]
pub fn l1_from_jaccard(x_weight: f64, y_weight: f64, j_sim: f64) -> f64 {
    let inter = j_sim * (x_weight + y_weight) / (1.0 + j_sim);
    inter / x_weight.min(y_weight)
}

// Count collisions (same id) between two MinHash sketches (id, rank) pairs.
// Only id is checked (matches original C++).
pub fn count_collisions(x: &[(u64, f64)], y: &[(u64, f64)]) -> u64 {
    assert_eq!(x.len(), y.len(), "sketches must be same length");
    let mut c = 0u64;
    for i in 0..x.len() {
        if x[i].0 == y[i].0 {
            c += 1;
        }
    }
    c
}

// Jaccard estimate from MinHash sketches: collisions / k
#[inline]
pub fn jaccard_estimate_from_minhashes(x: &[(u64, f64)], y: &[(u64, f64)]) -> f64 {
    count_collisions(x, y) as f64 / x.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weight_and_intersection() {
        let mut a = vec![(1, 0.4), (2, 0.1)];
        let mut b = vec![(1, 0.4), (3, 0.3)];
        a.sort_by_key(|p| p.0);
        b.sort_by_key(|p| p.0);
        assert!((weight(&a) - 0.5).abs() < 1e-12);
        assert!((weight(&b) - 0.7).abs() < 1e-12);
        assert!((intersection(&a, &b) - 0.4).abs() < 1e-12);
        let j = jaccard_similarity(&a, &b);
        assert!((j - 0.4 / (0.5 + 0.7 - 0.4)).abs() < 1e-12);
    }

    #[test]
    fn test_onebit_estimate() {
        let x = vec![true, false, true, true];
        let y = vec![true, true, false, true];
        let h = hamming_distance(&x, &y);
        assert_eq!(h, 2.0);
        let est = onebit_minhash_jaccard_estimate(&x, &y);
        // Just sanity: within [0,1]
        assert!(est >= 0.0 && est <= 1.0);
    }

    #[test]
    fn test_conversions() {
        let wx = 10.0;
        let wy = 8.0;
        let j = 0.4;
        let l1 = l1_from_jaccard(wx, wy, j);
        let j_back = jaccard_from_l1(wx, wy, l1);
        assert!((j - j_back).abs() < 1e-12);
    }
}