// Sergey's matrices Gamma, delta^B, and delta^G

use anyhow::Result;
use bitvec::vec::BitVec;
use nalgebra::SimdBool;

use crate::{phase_polynomial::PolynomialGraph, bit_matrix::BitMatrix};

use super::simulation_params::SimulationParams;
use bitvec::prelude::*;


pub struct LinearSystems {
    pub gamma: BitMatrix,
    pub delta_b: BitVec,
    pub delta_g: BitVec,
    pub x_r: BitVec,
}

impl LinearSystems {
    pub fn new(params: &SimulationParams, phase_graph: &PolynomialGraph) -> Self {
        let nodes = params.nodes as usize;
        let n_qubits = params.n_qubits as usize;
        let mut gamma = BitMatrix{};
        // for (b, g) in &phase_graph.bg_monomials {
        //     gamma = gamma.emplace_at(1, *b as usize, *g as usize);
        // }
        let delta_b = bitvec![usize, Lsb0; 0; params.nodes as usize ];
        let delta_g = bitvec![usize, Lsb0; 0; params.nodes as usize ];
        let x_r = bitvec![usize, Lsb0; 0; params.nodes as usize ];
        Self {
            gamma,
            delta_b,
            delta_g,
            x_r,
        }
    }

    // pub fn solve_if_gamma_null_space_quick_check(
    //     &self,
    //     s_b: &SparseBinMat,
    //     s_g: &SparseBinMat,
    // ) -> Result<Option<()>> {
    //     let s_b_difference = s_b.bitwise_xor_with(&self.delta_b)?.transposed();
    //     let s_b_xr_overlap = self.x_r.mul(&s_b_difference);
    //     let s_b_xr_overlap_even_parity = s_b_xr_overlap.number_of_ones() % 2 == 0;
    //     let s_g_difference = s_g.bitwise_xor_with(&self.delta_g)?.transposed();
    //     let s_g_xr_overlap = self.x_r.mul(&s_g_difference);
    //     let s_g_xr_overlap_even_parity = s_g_xr_overlap.number_of_ones() % 2 == 0;
    //     let not_in_nullspace = s_b_xr_overlap_even_parity && s_g_xr_overlap_even_parity;
    //     if not_in_nullspace {
    //         let solution = self
    //             .gamma
    //             .solve(&s_b_difference)
    //             .context("solving Gamma * X = sB^deltaB")?;
    //     }
    //     Ok(Some(()))
    // }

    pub fn update_with_flip_bit(
        &mut self,
        flip_bit: u32,
        phase_graph: &PolynomialGraph,
    ) {
        let flip_index = flip_bit as usize;
        let to_flip = !self.x_r[flip_index];
        self.x_r.set(flip_index, to_flip);
    }
}
