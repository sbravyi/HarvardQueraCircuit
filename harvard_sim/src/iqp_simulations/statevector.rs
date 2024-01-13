use anyhow::{Context, Result};
use std::ops::Mul;

use sparse_bin_mat::{SparseBinMat, SparseBinVec};

use crate::phase_polynomial::PolynomialGraph;

use super::simulation_params::SimulationParams;
use rand::Rng;

pub fn generate_random_statevector(params: &SimulationParams) -> SparseBinVec {
    let mut rng = rand::thread_rng();
    let rand: u32 = rng.gen_range(0..params.n_qubits);
    if rand == 0 {
        return SparseBinVec::zeros(0);
    }
    let n_bits = rand.ilog2() as usize;
    let mut nontrivial_positions = Vec::with_capacity(n_bits);
    for i in 0..n_bits {
        let bit_index = 1 << i;
        if rand & bit_index != 0 {
            nontrivial_positions.push(i);
        }
    }
    SparseBinVec::new(params.n_qubits as usize, nontrivial_positions)
}

// Sergey's matrices Gamma, delta^B, and delta^G

pub struct LinearSystems {
    pub gamma: SparseBinMat,
    pub delta_b: SparseBinVec,
    pub delta_g: SparseBinVec,
    pub x_r: SparseBinVec,
}

impl LinearSystems {
    pub fn new(params: &SimulationParams, phase_graph: &PolynomialGraph) -> Self {
        let nodes = params.nodes as usize;
        let n_qubits = params.n_qubits as usize;
        let mut gamma = SparseBinMat::zeros(n_qubits, n_qubits);
        for (b, g) in &phase_graph.bg_monomials {
            gamma = gamma.emplace_at(1, *b as usize, *g as usize);
        }
        let delta_b = SparseBinVec::zeros(nodes);
        let delta_g = SparseBinVec::zeros(nodes);
        let x_r = SparseBinVec::zeros(nodes);
        Self {
            gamma,
            delta_b,
            delta_g,
            x_r,
        }
    }

    pub fn is_gamma_null_space_quick_check(
        &self,
        s_b: &SparseBinVec,
        s_g: &SparseBinVec,
    ) -> Result<bool> {
        let s_b_difference = s_b.bitwise_xor_with(&self.delta_b)?;
        let s_b_xr_overlap = self.x_r.mul(&s_b_difference);
        let s_b_xr_overlap_even_parity = s_b_xr_overlap.is_zero();
        let s_g_difference = s_g.bitwise_xor_with(&self.delta_g)?;
        let s_g_xr_overlap = self.x_r.mul(&s_g_difference);
        let s_g_xr_overlap_even_parity = s_g_xr_overlap.is_zero();
        Ok(s_b_xr_overlap_even_parity && s_g_xr_overlap_even_parity)
    }

    pub fn update_with_flip_bit(
        mut self,
        flip_bit: u32,
        phase_graph: &PolynomialGraph,
    ) -> Result<Self> {
        let flip_bit_vec = SparseBinVec::new(self.x_r.len(), vec![flip_bit as usize]);
        self.x_r = self
            .x_r
            .bitwise_xor_with(&flip_bit_vec)
            .context("x_r[flip_bit]^= 1")?;
        println!("flip bit: {flip_bit}");
        for (b, g) in phase_graph.rbg_monomials.get(&flip_bit).unwrap() {
            let row = *b as usize;
            let column = *g as usize;
            println!("row: {row}, column: {column}");
            let old_val = self.gamma.get(row, column).unwrap();
            self.gamma =
                self.gamma
                    .emplace_at(if old_val.is_zero() { 1 } else { 0 }, row, column);
        }
        println!("flip bit rb: {flip_bit}");
        if phase_graph.rb_monomials.contains_key(&flip_bit) {
            for h in phase_graph.rb_monomials.get(&flip_bit).unwrap() {
                println!("h: {h}, nodes: {}", self.delta_b.len());
                let flip_h_vec = SparseBinVec::new(self.delta_b.len(), vec![*h as usize]);
                self.delta_b = self
                    .delta_b
                    .bitwise_xor_with(&flip_h_vec)
                    .context("deltaB[h]^=1")?;
            }
        }
        println!("flip bit rg: {flip_bit}");
        if phase_graph.rg_monomials.contains_key(&flip_bit) {
            for h in phase_graph.rg_monomials.get(&flip_bit).unwrap() {
                let flip_h_vec = SparseBinVec::new(self.delta_g.len(), vec![*h as usize]);
                self.delta_g = self
                    .delta_g
                    .bitwise_xor_with(&flip_h_vec)
                    .context("deltaB[h]^=1")?;
            }
        }
        Ok(self)
    }
}
