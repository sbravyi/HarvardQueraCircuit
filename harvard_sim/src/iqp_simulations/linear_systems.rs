// Sergey's matrices Gamma, delta^B, and delta^G

use bitvec::vec::BitVec;
use num_traits::Pow;

use crate::{
    bit_matrix::{
        backwards_substitution::BackwardsSubstitution, gauss_jordan::GaussJordan, matrix::BitMatrix,
    },
    phase_polynomial::PolynomialGraph,
};

use super::simulation_params::SimulationParams;
use bitvec::prelude::*;

pub struct LinearSystems {
    gamma: BitMatrix,
    delta_b: BitVec,
    delta_g: BitVec,
    x_r: BitVec,
    gj: GaussJordan,
    solver: BackwardsSubstitution,
    // caching allocations for
    // often used intermediary vectors
    sb_delta_b: BitVec,
    sg_delta_g: BitVec,
}

impl LinearSystems {
    // handles fast execution of the main linear system solving. Ensures that there are no allocations
    // in the main loop by holding onto all the allocations that are necessary as a one-time initialization
    // cost.
    pub fn new(params: &SimulationParams, phase_graph: &PolynomialGraph) -> Self {
        let nodes = params.nodes as usize;
        let mut gamma = BitMatrix::zeroes(nodes, nodes);
        for (b, g) in &phase_graph.bg_monomials {
            gamma.set(*b as usize, *g as usize, true);
        }
        let delta_b = bitvec![usize, Lsb0; 0; nodes];
        let delta_g = bitvec![usize, Lsb0; 0; nodes];
        let sb_delta_b = bitvec![usize, Lsb0; 0; nodes];
        let sg_delta_g = bitvec![usize, Lsb0; 0; nodes];
        let x_r = bitvec![usize, Lsb0; 0; nodes];
        let gj = GaussJordan::zero(nodes, nodes);
        let solver = BackwardsSubstitution::zero(nodes);
        Self {
            gamma,
            delta_b,
            delta_g,
            sb_delta_b,
            sg_delta_g,
            x_r,
            gj,
            solver,
        }
    }

    // extremely performance sensitive function - this is called an exponential amount of times
    pub fn solve_if_gamma_null_space_quick_check(
        &mut self,
        s_b: &BitVec,
        s_g: &BitVec,
        s_r: &BitVec,
    ) -> Option<f64> {
        let mut s_b_xr_overlap_bits = 0;
        // let min_length: usize = min([s_b.len(), self.delta_b.len(), self.x_r.len()]).unwrap();
        for (idx, bit) in s_b.iter().by_vals().enumerate() {
            let delta = self.delta_b[idx];
            let xord = bit ^ delta;
            unsafe {
                self.sb_delta_b.set_unchecked(idx, xord);
            }
            let x = self.x_r[idx];
            let overlap = xord & x;
            if overlap {
                s_b_xr_overlap_bits += 1;
            }
        }
        let s_b_xr_overlap_even_parity = s_b_xr_overlap_bits % 2 == 0;
        let mut s_g_xr_overlap_bits = 0;
        for (idx, bit) in s_g.iter().by_vals().enumerate() {
            let delta = self.delta_g[idx];
            let xord = bit ^ delta;
            unsafe {
                self.sg_delta_g.set_unchecked(idx, xord);
            }
            let x = self.x_r[idx];
            let overlap = xord & x;
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
            self.solver.solve(&self.gj, &self.sb_delta_b)?;
            let xg = &self.solver.solution;
            let has_amplitude_contributions = if rank == self.gj.number_of_columns {
                true
            } else {
                // TODO
                false
            };
            if !has_amplitude_contributions {
                return None;
            }
            let mut sg_delta_g_xg_overlap = 0;
            for (idx, bit) in xg.iter().by_vals().enumerate() {
                let sg_delta_g_i = self.sg_delta_g[idx];
                if sg_delta_g_i & bit {
                    sg_delta_g_xg_overlap += 1;
                }
            }
            let mut sr_xr_overlap = 0;
            for (idx, bit) in s_r.iter().by_vals().enumerate() {
                if bit & self.x_r[idx] {
                    sr_xr_overlap += 1;
                }
            }
            let phase_exponent = sg_delta_g_xg_overlap + sr_xr_overlap % 2;
            let phase: f64 = (-1.0).pow(phase_exponent);
            let amplitude_increment: f64 = phase / ((1 << rank) as f64);
            return Some(amplitude_increment);
        }
        None
    }

    // extremely performance sensitive function - this is called an exponential amount of times
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
