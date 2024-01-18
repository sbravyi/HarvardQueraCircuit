// Sergey's matrices Gamma, delta^B, and delta^G

use std::collections::{BTreeSet, HashMap};

use bitvec::vec::BitVec;
use num_traits::Pow;

use crate::{
    bit_matrix::{matrix::BitMatrix, sergey_solver::SergeySolver},
    phase_polynomial::PolynomialGraph,
};

use super::{simulation_params::SimulationParams, swap_check::SwapSymmetries};
use bitvec::prelude::*;

pub struct LinearSystems {
    pub gamma: BitMatrix,
    pub delta_b: BitVec,
    pub delta_g: BitVec,
    pub x_r: BitVec,
    pub solver: SergeySolver,
    // caching allocations for
    // often used intermediary vectors
    pub sb_delta_b: BitVec,
    pub sg_delta_g: BitVec,
    pub symmetry_checker: Option<SwapSymmetries>,
    pub symmetry_groups: HashMap<u16, BTreeSet<u16>>,
    pub dmitri_skips: BTreeSet<u16>,
    pub increments: HashMap<u16, f64>
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
        let solver = SergeySolver::zero(nodes);
        let symmetry_checker = if params.nodes <= 16 {
            Some(SwapSymmetries::new())
        } else {
            None
        };
        Self {
            gamma,
            delta_b,
            delta_g,
            sb_delta_b,
            sg_delta_g,
            x_r,
            solver,
            symmetry_checker,
            symmetry_groups: HashMap::new(),
            dmitri_skips: BTreeSet::new(),
            increments: HashMap::new()
        }
    }

    #[inline(never)]
    fn sergey_check(&mut self, s_b: &BitVec, s_g: &BitVec) -> Option<()> {
        let mut s_b_xr_overlap_bits = 0;
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
        if s_b_xr_overlap_even_parity && s_g_xr_overlap_even_parity {
            Some(())
        } else {
            None
        }
    }

    #[inline(never)]
    pub fn dmitri_check(&mut self) -> Option<BTreeSet<u16>> {
        if let Some(checker) = &mut self.symmetry_checker {
            if let Some(symmetries) = checker.check_for_symmetries() {
                // #[cfg(debug_assertions)]
                // println!("symmetries detected for {:?} - {:?}", checker.current_bit_pattern, symmetries);
                let mut clone = symmetries.clone();
                clone.insert(checker.current_bit_pattern);
                self.symmetry_groups.insert(checker.current_bit_pattern, clone);
                for symm in symmetries.clone() {
                    let mut clone = symmetries.clone();
                    clone.insert(checker.current_bit_pattern);
                    self.symmetry_groups.insert(symm, clone);
                }
                return Some(symmetries)
            } else {
                #[cfg(debug_assertions)]
                self.dmitri_skips.insert(checker.current_bit_pattern);
                // println!("dmitri check woulda skipped bit [{}]: group: [{:?}]", checker.current_bit_pattern, self.symmetry_groups.get(&checker.current_bit_pattern).unwrap());
                return None
            }
        }
        unreachable!()
    }

    // extremely performance sensitive function - this is called an exponential amount of times
    pub fn solve_if_gamma_null_space_quick_check(
        &mut self,
        s_b: &BitVec,
        s_g: &BitVec,
        s_r: &BitVec,
    ) -> Option<f64> {
        self.dmitri_check();
        self.sergey_check(s_b, s_g)?;
        let amplitude_multiplier = 1.0;
        self.solver.solve(&self.gamma, &self.sb_delta_b)?;
        if !self.solver.is_full_rank() && !self.solver.is_nullspace_codeword(&self.sg_delta_g) {
            return None;
        }
        let xg = &mut self.solver.solution;
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
        let rank = self.solver.rank();
        let amplitude_increment: f64 = phase / ((1 << rank) as f64);
        self.increments.insert(self.symmetry_checker.as_ref().unwrap().current_bit_pattern, amplitude_increment);
        Some(amplitude_increment * amplitude_multiplier)
    }

    // extremely performance sensitive function - this is called an exponential amount of times
    pub fn update_with_flip_bit(&mut self, flip_bit: u32, phase_graph: &PolynomialGraph) {
        let flip_index = flip_bit as usize;
        if let Some(ss) = &mut self.symmetry_checker {
            ss.increment_bit(Some(flip_bit));
            // ss.increment_bit(None);
        }
        let to_flip = !self.x_r[flip_index];
        unsafe {
            self.x_r.set_unchecked(flip_index, to_flip);
        }
        for (h0, h1) in &phase_graph.rbg_monomials[&flip_bit] {
            self.gamma.flip(*h0 as usize, *h1 as usize);
        }
        let rb_monomials = &phase_graph.rb_monomials[&flip_bit];
        for h in rb_monomials {
            let h_index = *h as usize;
            let to_flip = !self.delta_b[h_index];
            unsafe {
                self.delta_b.set_unchecked(h_index, to_flip);
            }
        }
        let rg_monomials = &phase_graph.rg_monomials[&flip_bit];
        for h in rg_monomials {
            let h_index = *h as usize;
            let to_flip = !self.delta_g[h_index];
            unsafe {
                self.delta_g.set_unchecked(h_index, to_flip);
            }
        }
    }
}
