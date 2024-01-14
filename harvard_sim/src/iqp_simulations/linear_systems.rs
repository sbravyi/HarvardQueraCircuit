// Sergey's matrices Gamma, delta^B, and delta^G

use anyhow::Result;
use bitvec::vec::BitVec;
use itertools::{min, Itertools};

use crate::{
    bit_matrix::{gauss_jordan::GaussJordan, matrix::BitMatrix},
    phase_polynomial::PolynomialGraph,
};

use super::simulation_params::SimulationParams;
use bitvec::prelude::*;

pub struct LinearSystems {
    pub gamma: BitMatrix,
    pub delta_b: BitVec,
    pub delta_g: BitVec,
    pub x_r: BitVec,
    pub gj: GaussJordan,
}

impl LinearSystems {
    pub fn new(params: &SimulationParams, phase_graph: &PolynomialGraph) -> Self {
        let nodes = params.nodes as usize;
        let mut gamma = BitMatrix::zeroes(nodes, nodes);
        for (b, g) in &phase_graph.bg_monomials {
            gamma.set(*b as usize, *g as usize, true);
        }
        let delta_b = bitvec![usize, Lsb0; 0; nodes];
        let delta_g = bitvec![usize, Lsb0; 0; nodes];
        let x_r = bitvec![usize, Lsb0; 0; nodes];
        let mut gj = GaussJordan::zero(nodes, nodes);
        Self {
            gamma,
            delta_b,
            delta_g,
            x_r,
            gj,
        }
    }

    pub fn solve_if_gamma_null_space_quick_check(
        &mut self,
        s_b: &BitVec,
        s_g: &BitVec,
    ) -> Result<Option<usize>> {
        let mut s_b_xr_overlap_bits = 0;
        let min_length = min([s_b.len(), self.delta_b.len(), self.x_r.len()]).unwrap();
        for (idx, bit) in s_b[0..min_length].iter().by_vals().enumerate() {
            let delta = self.delta_b[idx];
            let x = self.x_r[idx];
            let overlap = (bit ^ delta) & x;
            if overlap {
                s_b_xr_overlap_bits += 1;
            }
        }
        let s_b_xr_overlap_even_parity = s_b_xr_overlap_bits % 2 == 0;
        let mut s_g_xr_overlap_bits = 0;
        for (idx, bit) in s_g[0..min_length].iter().by_vals().enumerate() {
            let delta = self.delta_b[idx];
            let x = self.x_r[idx];
            let overlap = (bit ^ delta) & x;
            if overlap {
                s_g_xr_overlap_bits += 1;
            }
        }
        let s_g_xr_overlap_even_parity = s_g_xr_overlap_bits % 2 == 0;
        let not_in_nullspace = s_b_xr_overlap_even_parity && s_g_xr_overlap_even_parity;
        if not_in_nullspace {
            self.gj.copy_from_matrix(&self.gamma);
            self.gj.go_to_echelon_form();
            let rank = self.gj.rank();
            return Ok(Some(rank));
        }
        Ok(None)
    }

    pub fn update_with_flip_bit(&mut self, flip_bit: u32, phase_graph: &PolynomialGraph) {
        let flip_index = flip_bit as usize;
        let to_flip = !self.x_r[flip_index];
        self.x_r.set(flip_index, to_flip);
        for (h0, h1) in &phase_graph.rbg_monomials[&flip_bit] {
            self.gamma.flip(*h0 as usize, *h1 as usize);
        }
        let rb_monomials = &phase_graph.rb_monomials[&flip_bit];
        for h in rb_monomials {
            let h_index = *h as usize;
            let to_flip = !self.delta_b[h_index];
            unsafe {
                self.delta_b.set_unchecked(flip_index, to_flip);
            }
        }
        let rg_monomials = &phase_graph.rg_monomials[&flip_bit];
        for h in rg_monomials {
            let h_index = *h as usize;
            let to_flip = !self.delta_g[h_index];
            unsafe {
                self.delta_g.set_unchecked(flip_index, to_flip);
            }
        }
    }
}
