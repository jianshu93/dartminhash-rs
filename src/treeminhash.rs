//! TreeMinHash for weighted Jaccard sketching.
//!
//! This is a Rust port of Otmar Ertl's TreeMinHash idea from
//! `weighted_minwise_hashing.hpp`, adapted to the API style used by this crate:
//! sparse weighted vectors are `&[(u64, f64)]`, and sketches are `Vec<(u64, f64)>`.
//!
//! Randomness is provided by simple tabulation hashing (`tab_hash::Tab64Simple`).
//! Each logical random stream is addressed by `(feature_id, stream_id)` and then
//! expanded with a counter.  This keeps sketching stateless and deterministic.

use std::f64::INFINITY;
use tab_hash::Tab64Simple;

use crate::hash_utils::tab64_from_rng;
use crate::rng_utils::MtRng;

/// Same sketch representation as DartMinHash: k slots of `(id, rank)`.
pub type MinHashSketch = Vec<(u64, f64)>;

#[derive(Clone, Copy, Debug)]
struct Node {
    lower_bound: f64,
    inv_rate: f64,
    ratio: f64,
}

impl Node {
    #[inline]
    fn new(lower_bound: f64, mid_bound: f64, upper_bound: f64) -> Self {
        let inv_rate = 1.0 / (upper_bound - lower_bound);
        let ratio = inv_rate * (mid_bound - lower_bound);
        Self {
            lower_bound,
            inv_rate,
            ratio,
        }
    }
}

/// Build the same implicit binary tree
///
/// Leaves correspond to exponentially spaced weight intervals. Internal nodes
/// store the lower bound of their interval and the rate parameters needed for
/// exponential thinning.
fn pre_calculate_tree(factor: f64, max: f64) -> Vec<Node> {
    assert!(
        max > 0.0 && max.is_finite() || max == f64::MAX,
        "max must be positive"
    );
    assert!(factor > 0.0 && factor < 1.0, "factor must be in (0, 1)");

    #[derive(Clone, Copy, Debug, Default)]
    struct TmpNode {
        min_idx: u32,
        mid_idx: u32,
        max_idx: u32,
    }

    let mut boundaries = Vec::<f64>::new();
    boundaries.push(max);
    loop {
        let prev = *boundaries.last().unwrap();
        if prev == 0.0 {
            break;
        }
        let next = (factor * prev).min(next_toward_zero(prev));
        boundaries.push(next);
        if next == 0.0 {
            break;
        }
    }
    boundaries.reverse();

    let num_nodes = boundaries.len() - 1;
    assert!(num_nodes > 0);
    assert!(num_nodes <= (u32::MAX as usize));

    let mut counts = vec![0u32; 2 * num_nodes - 1];
    for i in 0..num_nodes {
        counts[num_nodes - 1 + i] = 1;
    }
    for idx in (0..=(num_nodes - 2)).rev() {
        let left = 2 * idx + 1;
        let right = 2 * idx + 2;
        counts[idx] = counts[left] + counts[right];
    }

    let mut tmp = vec![TmpNode::default(); 2 * num_nodes - 1];
    tmp[0] = TmpNode {
        min_idx: 0,
        mid_idx: 0,
        max_idx: num_nodes as u32,
    };

    for idx in 0..(num_nodes - 1) {
        let left = 2 * idx + 1;
        let right = 2 * idx + 2;
        let mid = tmp[idx].min_idx + counts[left];
        tmp[idx].mid_idx = mid;
        tmp[left] = TmpNode {
            min_idx: tmp[idx].min_idx,
            mid_idx: 0,
            max_idx: mid,
        };
        tmp[right] = TmpNode {
            min_idx: mid,
            mid_idx: 0,
            max_idx: tmp[idx].max_idx,
        };
    }

    let mut tree = Vec::with_capacity(2 * num_nodes - 1);
    for t in tmp {
        let mid_bound = if t.mid_idx == 0 {
            -1.0
        } else {
            boundaries[t.mid_idx as usize]
        };
        tree.push(Node::new(
            boundaries[t.min_idx as usize],
            mid_bound,
            boundaries[t.max_idx as usize],
        ));
    }
    tree
}

/// `std::nexttoward(x, 0)` equivalent for positive finite f64 values.
#[inline]
fn next_toward_zero(x: f64) -> f64 {
    debug_assert!(x >= 0.0);
    if x == 0.0 {
        0.0
    } else {
        f64::from_bits(x.to_bits() - 1)
    }
}

/// A deterministic random stream based on simple tabulation hashing.
///
/// The C++ reference uses a bit-stream RNG seeded by `(id, stream_id)`.  Here we
/// generate each 64-bit word by tab-hashing a mixed key derived from
/// `(id, stream_id, counter)`.  All hashers are `Tab64Simple` tables seeded in
/// `TreeMinHash::new_mt`.
#[derive(Clone)]
struct TabStream<'a> {
    h0: &'a Tab64Simple,
    h1: &'a Tab64Simple,
    id: u64,
    stream_id: u64,
    counter: u64,
}

impl<'a> TabStream<'a> {
    #[inline]
    fn new(h0: &'a Tab64Simple, h1: &'a Tab64Simple, id: u64, stream_id: u64) -> Self {
        Self {
            h0,
            h1,
            id,
            stream_id,
            counter: 0,
        }
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let c = self.counter;
        self.counter = self.counter.wrapping_add(1);

        // Keep the construction simple-tabulation based, but domain-separated.
        let a = self.id ^ c.wrapping_mul(0x9e37_79b9_7f4a_7c15);
        let b = self.stream_id.rotate_left(17) ^ c.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        self.h0.hash(a) ^ self.h1.hash(b)
    }

    #[inline]
    fn uniform_open01(&mut self) -> f64 {
        // Use the top 53 bits. Add 0.5 so the result is strictly inside (0, 1).
        const DEN: f64 = (1u64 << 53) as f64;
        let v = self.next_u64() >> 11;
        ((v as f64) + 0.5) / DEN
    }

    #[inline]
    fn exponential1(&mut self) -> f64 {
        -self.uniform_open01().ln()
    }

    #[inline]
    fn bernoulli(&mut self, p: f64) -> bool {
        if p <= 0.0 {
            false
        } else if p >= 1.0 {
            true
        } else {
            self.uniform_open01() < p
        }
    }

    #[inline]
    fn uniform_index(&mut self, upper_exclusive: u32) -> u32 {
        debug_assert!(upper_exclusive > 0);
        // Rejection sampling avoids modulo bias. For unit tests and sketches this
        // branch almost always exits on the first iteration.
        let bound = u64::MAX - (u64::MAX % (upper_exclusive as u64));
        loop {
            let x = self.next_u64();
            if x < bound {
                return (x % (upper_exclusive as u64)) as u32;
            }
        }
    }
}

/// Lazy random permutation stream over `0..m`.
///
/// This is the standard partial Fisher-Yates shuffle: the kth call samples one
/// still-unused index uniformly.  It matches the role of `PermutationStream` in
/// the C++ reference implementation.
#[derive(Clone, Debug)]
struct PermutationStream {
    permutation: Vec<u32>,
    pos: u32,
}

impl PermutationStream {
    fn new(m: u32) -> Self {
        let permutation = (0..m).collect();
        Self {
            permutation,
            pos: 0,
        }
    }

    #[inline]
    fn reset(&mut self) {
        self.pos = 0;
        // Restore only the prefix that may have been modified. Simpler and less
        // bug-prone than keeping a touched list; O(m) is acceptable because reset
        // happens only when a leaf interval produces candidates.
        for (i, v) in self.permutation.iter_mut().enumerate() {
            *v = i as u32;
        }
    }

    #[inline]
    fn next(&mut self, rng: &mut TabStream<'_>) -> u32 {
        let n = self.permutation.len() as u32;
        debug_assert!(self.pos < n);
        let j = self.pos + rng.uniform_index(n - self.pos);
        self.permutation.swap(self.pos as usize, j as usize);
        let out = self.permutation[self.pos as usize];
        self.pos += 1;
        out
    }
}

/// TreeMinHash sketcher for weighted Jaccard similarity.
pub struct TreeMinHash {
    k: u32,
    h0: Tab64Simple,
    h1: Tab64Simple,
    tree: Vec<Node>,
    num_non_leaf_nodes: u32,
    initial_limit_factor: f64,
    factors: Vec<f64>,
}

impl TreeMinHash {
    /// Build with MT19937-seeded simple tabulation hash tables.
    ///
    /// Defaults match the C++ constructor: `max = f64::MAX`, `factor = 0.5`,
    /// and first-run success probability `0.9`.
    pub fn new_mt(rng: &mut MtRng, k: u64) -> Self {
        Self::with_params_mt(rng, k, f64::MAX, 0.5, 0.9)
    }

    pub fn with_params_mt(
        rng: &mut MtRng,
        k: u64,
        max: f64,
        factor: f64,
        success_probability_first_run: f64,
    ) -> Self {
        assert!(k > 0 && k <= (u32::MAX as u64), "k must fit into u32");
        assert!(success_probability_first_run > 0.0 && success_probability_first_run < 1.0);

        let tree = pre_calculate_tree(factor, max);
        let num_non_leaf_nodes = (tree.len() - (tree.len() + 1) / 2) as u32;
        let k_f = k as f64;
        // Equivalent to the C++ expression:
        // -log(-expm1(log(p) / m)) * m
        let initial_limit_factor =
            -(-(success_probability_first_run.ln() / k_f).exp_m1()).ln() * k_f;

        let mut factors = Vec::with_capacity(k.saturating_sub(1) as usize);
        for i in 0..(k as u32).saturating_sub(1) {
            factors.push(k_f / (k_f - (i as f64) - 1.0));
        }

        let h0 = tab64_from_rng(rng);
        let h1 = tab64_from_rng(rng);

        Self {
            k: k as u32,
            h0,
            h1,
            tree,
            num_non_leaf_nodes,
            initial_limit_factor,
            factors,
        }
    }

    /// Return k weighted MinHash slots.  The rank component is useful for
    /// debugging and compatibility; Jaccard estimation only needs id collisions.
    pub fn sketch(&self, x: &[(u64, f64)]) -> MinHashSketch {
        let weight_sum: f64 = x.iter().filter(|(_, w)| *w > 0.0).map(|(_, w)| *w).sum();
        if !(weight_sum > 0.0) || !weight_sum.is_finite() {
            return vec![(0, INFINITY); self.k as usize];
        }

        let limit_increment = self.initial_limit_factor / weight_sum;
        let mut limit = limit_increment;
        let mut result = vec![(0u64, limit); self.k as usize];
        let mut buffer: Vec<(f64, u32)> = Vec::with_capacity(self.num_non_leaf_nodes as usize);
        let mut permutation_stream = PermutationStream::new(self.k);

        loop {
            for &(id, w) in x {
                if !(w > 0.0) || !w.is_finite() {
                    continue;
                }

                buffer.clear();
                let mut node_idx = 0u32;
                let mut rng = self.rng(id, node_idx as u64);
                let mut point = rng.exponential1() * self.tree[node_idx as usize].inv_rate;
                if !(point < limit) {
                    continue;
                }

                loop {
                    while node_idx < self.num_non_leaf_nodes
                        && self.tree[node_idx as usize].lower_bound < w
                    {
                        let node = self.tree[node_idx as usize];
                        let inherit_to_left = rng.bernoulli(node.ratio);

                        node_idx <<= 1;
                        let sibling_idx = node_idx + 1 + (inherit_to_left as u32);
                        node_idx += 2 - (inherit_to_left as u32);

                        let sibling_node = self.tree[sibling_idx as usize];
                        let sibling_point = point + rng.exponential1() * sibling_node.inv_rate;
                        if sibling_point < limit && sibling_node.lower_bound < w {
                            buffer.push((sibling_point, sibling_idx));
                        }
                    }

                    let node = self.tree[node_idx as usize];
                    if node.lower_bound < w {
                        let inv_rate = node.inv_rate;
                        let acceptance_probability = (w - node.lower_bound) * inv_rate;

                        permutation_stream.reset();
                        for kk in 0..self.k {
                            let next_point = if (kk as usize) < self.factors.len() {
                                point + rng.exponential1() * inv_rate * self.factors[kk as usize]
                            } else {
                                INFINITY
                            };

                            let idx = permutation_stream.next(&mut rng) as usize;
                            if point < result[idx].1 {
                                if !(acceptance_probability < 1.0) {
                                    result[idx] = (id, point);
                                } else {
                                    let stream_id = ((node_idx as u64) << 32) | (idx as u64);
                                    let mut rng2 = self.rng(id, stream_id);
                                    let mut p = point;
                                    while p < result[idx].1 {
                                        if rng2.uniform_open01() < acceptance_probability {
                                            result[idx] = (id, p);
                                            break;
                                        }
                                        p += rng2.exponential1() * inv_rate * (self.k as f64);
                                    }
                                }
                            }

                            if !(next_point < limit) {
                                break;
                            }
                            point = next_point;
                        }
                    }

                    if let Some((p, ni)) = buffer.pop() {
                        point = p;
                        node_idx = ni;
                        rng = self.rng(id, node_idx as u64);
                    } else {
                        break;
                    }
                }
            }

            if result.iter().all(|&(_, r)| r != limit) {
                return result;
            }

            let old_limit = limit;
            limit += limit_increment;
            for slot in &mut result {
                if slot.1 == old_limit {
                    slot.1 = limit;
                }
            }
        }
    }

    #[inline]
    fn rng(&self, id: u64, stream_id: u64) -> TabStream<'_> {
        TabStream::new(&self.h0, &self.h1, id, stream_id)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use rand_core::RngCore;

    use crate::{
        dartminhash::DartMinHash,
        rng_utils::{MtRng, mt_from_seed},
        similarity::{jaccard_estimate_from_minhashes, jaccard_similarity},
        treeminhash::TreeMinHash,
    };

    /// Uniform(0,1) using the same MT19937 rng as the existing DartMinHash tests.
    fn uniform01(rng: &mut MtRng) -> f64 {
        mt19937::gen_res53(rng)
    }

    /// Generate a random weighted set:
    /// Pick L0 distinct random indices, draw L0-1 splitters, sort them, and use
    /// the gaps times L1 as weights. Returns sorted by id.
    fn generate_weighted_set(l0: u64, l1: f64, rng: &mut MtRng) -> Vec<(u64, f64)> {
        let mut elements = HashSet::with_capacity(l0 as usize);
        while elements.len() < l0 as usize {
            elements.insert(rng.next_u64());
        }

        let mut z: Vec<f64> = (0..(l0 - 1)).map(|_| uniform01(rng)).collect();
        z.push(1.0);
        z.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut prev = 0.0;
        let mut j = 0usize;
        let mut out = Vec::with_capacity(l0 as usize);
        for idx in elements {
            let w = l1 * (z[j] - prev);
            out.push((idx, w));
            prev = z[j];
            j += 1;
        }
        out.sort_by_key(|p| p.0);
        out
    }

    /// Generate Y from X with target relative overlap:
    /// y = relative_overlap * x plus the remaining mass as a new element.
    fn generate_similar_weighted_set(
        x: &[(u64, f64)],
        relative_overlap: f64,
        rng: &mut MtRng,
    ) -> Vec<(u64, f64)> {
        let free_id;
        loop {
            let candidate = rng.next_u64();
            if x.binary_search_by_key(&candidate, |p| p.0).is_err() {
                free_id = candidate;
                break;
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
    fn treeminhash_approximates_weighted_jaccard() {
        let mut rng = mt_from_seed(1337);

        // Same structure as the DartMinHash test, with slightly reduced size so
        // `cargo test` remains practical in debug mode.
        let l0 = 50_000;
        let l1 = 10_000.0;
        let k = 4096;

        let tmh = TreeMinHash::new_mt(&mut rng, k);
        let x = generate_weighted_set(l0, l1, &mut rng);
        assert_eq!(x.len(), l0 as usize);

        let targets = [
            0.99, 0.96, 0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1,
            0.05, 0.01,
        ];

        for &rel in &targets {
            let y = generate_similar_weighted_set(&x, rel, &mut rng);
            let j_true = jaccard_similarity(&x, &y);
            println!("true weighted Jaccard: {:?}", j_true);

            let sk_x = tmh.sketch(&x);
            let sk_y = tmh.sketch(&y);
            assert_eq!(sk_x.len(), k as usize);
            assert_eq!(sk_y.len(), k as usize);

            let j_est = jaccard_estimate_from_minhashes(&sk_x, &sk_y);
            println!("estimated weighted Jaccard: {:?}", j_est);
            let sd = (j_true * (1.0 - j_true) / (k as f64)).sqrt();
            let tol = (3.2 * sd).max(1.25 / (k as f64).sqrt());
            let err = (j_true - j_est).abs();
            assert!(
                err <= tol,
                "TMH: rel_overlap={rel}, true={j_true:.6}, est={j_est:.6}, err={err:.6}, tol={tol:.6}"
            );
        }
    }

    #[test]
    fn treeminhash_approximates_weighted_jaccard_sparse() {
        let mut rng = mt_from_seed(2025);

        let l0 = 5_000;
        let l1 = 3_000.0;
        let k = 4096;

        let tmh = TreeMinHash::new_mt(&mut rng, k);
        let x = generate_weighted_set(l0, l1, &mut rng);
        assert_eq!(x.len(), l0 as usize);

        let targets = [
            0.99, 0.96, 0.93, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2, 0.1,
            0.05, 0.01,
        ];

        for &rel in &targets {
            let y = generate_similar_weighted_set(&x, rel, &mut rng);
            let j_true = jaccard_similarity(&x, &y);

            let sk_x = tmh.sketch(&x);
            let sk_y = tmh.sketch(&y);
            assert_eq!(sk_x.len(), k as usize);
            assert_eq!(sk_y.len(), k as usize);

            let j_est = jaccard_estimate_from_minhashes(&sk_x, &sk_y);
            let sd = (j_true * (1.0 - j_true) / (k as f64)).sqrt();
            let tol = (3.2 * sd).max(1.25 / (k as f64).sqrt());
            let err = (j_true - j_est).abs();
            assert!(
                err <= tol,
                "TMH sparse: true={j_true:.6}, est={j_est:.6}, err={err:.6}, tol={tol:.6}"
            );
        }
    }


    /// Large raw-count / absolute-weight simulation.
    ///
    /// This is intentionally ignored because it is a timing-style stress test,
    /// not a small unit test. It mirrors the companion test in dartminhash.rs
    /// and is useful for testing the regime where raw counts create very large
    /// total set mass.
    ///
    /// Run with:
    ///
    ///     cargo test treeminhash_large_weight_sum_vs_dartminhash --release -- --ignored --nocapture
    #[test]
    #[ignore]
    fn treeminhash_large_weight_sum_vs_dartminhash() {
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

        let mut tmh_rng = mt_from_seed(12_345);
        let tmh = TreeMinHash::new_mt(&mut tmh_rng, k);
        let tmh_start = Instant::now();
        let sk_tmh = black_box(tmh.sketch(black_box(&x)));
        let tmh_elapsed = tmh_start.elapsed();
        assert_eq!(sk_tmh.len(), k as usize);
        assert!(sk_tmh.iter().all(|&(_, r)| r.is_finite()));

        let mut dmh_rng = mt_from_seed(12_345);
        let dmh = DartMinHash::new_mt(&mut dmh_rng, k);
        let dmh_start = Instant::now();
        let sk_dmh = black_box(dmh.sketch(black_box(&x)));
        let dmh_elapsed = dmh_start.elapsed();
        assert_eq!(sk_dmh.len(), k as usize);
        assert!(sk_dmh.iter().all(|&(_, r)| r.is_finite()));

        println!("TreeMinHash elapsed: {:?}", tmh_elapsed);
        println!("DartMinHash elapsed: {:?}", dmh_elapsed);
        if tmh_elapsed.as_nanos() > 0 {
            println!(
                "DMH / TMH elapsed ratio: {:.3}x",
                dmh_elapsed.as_secs_f64() / tmh_elapsed.as_secs_f64()
            );
        }

        // Do not assert timing. Runtime ratios are machine/build dependent.
        // This test is an explicit stress benchmark for the high-weight-sum
        // regime, not a correctness proof.
    }

    #[test]
    fn treeminhash_is_deterministic_for_fixed_seed() {
        let x = vec![(1, 0.25), (3, 0.75), (10, 1.5), (42, 0.125)];
        let mut rng1 = mt_from_seed(7);
        let mut rng2 = mt_from_seed(7);
        let tmh1 = TreeMinHash::new_mt(&mut rng1, 256);
        let tmh2 = TreeMinHash::new_mt(&mut rng2, 256);
        assert_eq!(tmh1.sketch(&x), tmh2.sketch(&x));
    }

    #[test]
    fn treeminhash_empty_or_zero_input_returns_empty_slots() {
        let mut rng = mt_from_seed(11);
        let tmh = TreeMinHash::new_mt(&mut rng, 16);
        let sk = tmh.sketch(&[(1, 0.0), (2, -1.0)]);
        assert_eq!(sk.len(), 16);
        assert!(sk.iter().all(|&(id, rank)| id == 0 && rank.is_infinite()));
    }
}
