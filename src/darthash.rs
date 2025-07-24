//! DartHash: produces darts (id, rank) from a weighted vector.

use std::f64::INFINITY;

use crate::hash_utils::*;
use crate::rng_utils::MtRng;
use tab_hash::{Tab32Simple, Tab64Simple};

// A single dart = (hashed_id, rank)
pub type Dart = (u64, f64);

pub struct DartHash {
    t: u64,
    // 32-bit tabulation hashers
    t_nu: Tab32Simple,
    t_rho: Tab32Simple,
    t_w: Tab32Simple,
    t_r: Tab32Simple,
    // 64-bit tabulation hashers
    t_i: Tab64Simple,
    t_p: Tab64Simple,
    t_q: Tab64Simple,
    f_h: Tab64Simple,
    m_h: Tab64Simple,
    // Precomputed powers of 2 and poisson CDF
    powers_of_two: Vec<f64>,
    neg_powers_of_two: Vec<f64>,
    poisson_cdf: Vec<f64>,
}

impl DartHash {
    // t: expected number of darts (usually k ln k + 2k)
    pub fn new_mt(rng: &mut MtRng, t: u64) -> Self {
        let t_nu = tab32_from_rng(rng);
        let t_rho = tab32_from_rng(rng);
        let t_w = tab32_from_rng(rng);
        let t_r = tab32_from_rng(rng);

        let t_i = tab64_from_rng(rng);
        let t_p = tab64_from_rng(rng);
        let t_q = tab64_from_rng(rng);
        let f_h = tab64_from_rng(rng);
        let m_h = tab64_from_rng(rng);

        // Precompute 2^k up to ~1000
        let mut pow2 = Vec::with_capacity(1000);
        let mut p = 1.0;
        for _ in 0..1000 {
            pow2.push(p);
            p *= 2.0;
        }
        let mut neg_pow2 = Vec::with_capacity(1000);
        let mut q = 1.0;
        for _ in 0..1000 {
            neg_pow2.push(q);
            q *= 0.5;
        }

        // Poisson(1) CDF up to 100
        let mut poisson_cdf = Vec::with_capacity(100);
        let mut pdf = (-1.0f64).exp();
        let mut cdf = pdf;
        for i in 0..100 {
            poisson_cdf.push(cdf);
            pdf = pdf / ((i + 1) as f64);
            cdf += pdf;
        }

        Self {
            t,
            t_nu,
            t_rho,
            t_w,
            t_r,
            t_i,
            t_p,
            t_q,
            f_h,
            m_h,
            powers_of_two: pow2,
            neg_powers_of_two: neg_pow2,
            poisson_cdf,
        }
    }

    // Generate darts for a weighted vector x.
    // x: vector of (feature_id, weight)
    // theta: search parameter (default 1.0)
    pub fn darts(&self, x: &[(u64, f64)], theta: f64) -> Vec<Dart> {
        let mut darts = Vec::with_capacity((2 * self.t) as usize);
        let total_w = total_weight(x);
        if total_w == 0.0 {
            return darts;
        }
        let max_rank = theta / total_w;
        let t_inv = 1.0 / (self.t as f64);
        let rho_upper = ((1.0 + max_rank).log2().floor()).max(0.0) as u32;

        for &(i, xi) in x {
            if xi <= 0.0 { continue; }

            let i_hash = self.t_i.hash(i);
            let nu_upper = ((1.0 + (self.t as f64) * xi).log2().floor()).max(0.0) as u32;

            for nu in 0..=nu_upper {
                let nu_hash = self.t_nu.hash(nu as u32);
                for rho in 0..=rho_upper {
                    let region_hash = (nu_hash as u64) ^ (self.t_rho.hash(rho as u32) as u64);

                    let two_nu = self.powers_of_two[nu as usize];
                    let two_rho = self.powers_of_two[rho as usize];
                    let w_base = (two_nu - 1.0) * t_inv;
                    let r_base = two_rho - 1.0;

                    let delta_nu = two_nu * t_inv * self.neg_powers_of_two[rho as usize];
                    let delta_rho = two_rho * self.neg_powers_of_two[nu as usize];

                    let mut w0 = w_base;
                    let w_max = if rho < 32 { 1u32 << rho } else { 1u32 << 31 };
                    for w in 0..w_max {
                        if xi < w0 { break; }
                        let w_hash = self.t_w.hash(w);
                        let mut r0 = r_base;
                        let r_max = if nu < 32 { 1u32 << nu } else { 1u32 << 31 };

                        for r in 0..r_max {
                            if max_rank < r0 { break; }

                            let area_hash = (w_hash as u64) ^ (self.t_r.hash(r) as u64);
                            let z = i_hash ^ region_hash ^ area_hash;

                            // Poisson draw via CDF table
                            let p_z = to_unit(self.t_p.hash(z));
                            let mut p_count: usize = 0;
                            while p_count < self.poisson_cdf.len() && p_z > self.poisson_cdf[p_count] {
                                p_count += 1;
                            }

                            let mut q_idx: usize = 0;
                            while q_idx < p_count {
                                // combine z & q to make unique
                                let z_q = z ^ ((q_idx as u64) << 56)
                                    ^ ((q_idx as u64) << 48)
                                    ^ ((q_idx as u64) << 40)
                                    ^ ((q_idx as u64) << 32)
                                    ^ ((q_idx as u64) << 24)
                                    ^ ((q_idx as u64) << 16)
                                    ^ ((q_idx as u64) << 8)
                                    ^ (q_idx as u64);

                                let (u_w, u_r) = to_units(self.t_q.hash(z_q));
                                let weight = w0 + delta_nu * u_w;
                                let rank   = r0 + delta_rho * u_r;

                                if weight < xi && rank < max_rank {
                                    darts.push((self.f_h.hash(z_q), rank));
                                }
                                q_idx += 1;
                            }

                            r0 += delta_rho;
                        }
                        w0 += delta_nu;
                    }
                }
            }
        }
        darts
    }

    // Convert darts to k buckets, keep min rank in each
    pub fn minhash(&self, x: &[(u64, f64)], k: u64) -> Vec<Dart> {
        let darts = self.darts(x, 1.0);
        let mut mh = vec![(0u64, INFINITY); k as usize];
        for &(id, rank) in &darts {
            let j = (self.m_h.hash(id) % k) as usize;
            if rank < mh[j].1 {
                mh[j] = (id, rank);
            }
        }
        mh
    }

    pub fn onebit_minhash(&self, x: &[(u64, f64)], k: u64) -> Vec<bool> {
        let mh = self.minhash(x, k);
        mh.into_iter().map(|(id, _)| (id & 1) == 1).collect()
    }
}