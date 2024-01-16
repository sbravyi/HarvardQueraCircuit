use super::{
    iqp_circuit::build_iqp_circuit, simulation::Simulation, simulation_params::SimulationParams,
    statevector::generate_random_statevector,
};
use crate::{
    gray_code_flip_bit::GrayCodeFlipBit,
    iqp_simulations::{iqp_circuit::QubitColoringIndexes, linear_systems::LinearSystems}
};
use anyhow::{Context, Result};
use bitvec::vec::BitVec;
use itertools::Itertools;
use std::time::Instant;

#[cfg(debug_assertions)]
use crate::debug_bv::debug_linear_system;

pub struct CPUSmallIntSimulation {
    params: SimulationParams,
}

impl CPUSmallIntSimulation {
    pub fn new(boolean_cube_dimension: u32) -> Self {
        Self {
            params: SimulationParams::new(boolean_cube_dimension),
        }
    }
}

impl Simulation for CPUSmallIntSimulation {
    fn run(&mut self) -> Result<()> {
        let start = Instant::now();
        let s: BitVec = generate_random_statevector(&self.params);
        let (phase_polynomial, coloring) =
            build_iqp_circuit(&self.params).context("Building IQP circuit")?;
        let QubitColoringIndexes { red, green, blue } = coloring;
        let phase_graph = phase_polynomial
            .into_polynomial_graph()
            .context("converting phase polynomial into graph")?;
        let gc_flip_bit = GrayCodeFlipBit::new(self.params.nodes);
        let s_r: BitVec = red.iter().map(|idx| s[*idx as usize]).collect();
        let s_b: BitVec = blue.iter().map(|idx| s[*idx as usize]).collect();
        let s_g: BitVec = green.iter().map(|idx| s[*idx as usize]).collect();
        let mut amplitude: f64 = 0.0;
        let mut ls = LinearSystems::new(&self.params, &phase_graph);
        for flip_bit in gc_flip_bit {
            #[cfg(debug_assertions)]
            println!("BEFORE SOLVING");
            #[cfg(debug_assertions)]
            println!("flip bit: {flip_bit}");
            #[cfg(debug_assertions)]
            debug_linear_system(&ls);
            if let Some(amplitude_increment) =
                ls.solve_if_gamma_null_space_quick_check(&s_b, &s_g, &s_r)
            {
            amplitude += amplitude_increment;
            #[cfg(debug_assertions)]
            println!("AMPLITUDE UPDATE: {:?}", (flip_bit, amplitude, amplitude_increment));
            amplitude += 0.0;
            }
            #[cfg(debug_assertions)]
            println!("AFTER SOLVING");
            #[cfg(debug_assertions)]
            debug_linear_system(&ls);
            ls.update_with_flip_bit(flip_bit, &phase_graph);
        }
        amplitude /= (1 << self.params.nodes) as f64;
        let end = Instant::now();
        log::debug!("Time to execute: {:#?}", end - start);
        println!(
            "Amplitude <s|U|00...0> (S = <{}|) ::= {:#?}",
            s.iter()
                .by_vals()
                .map(|x| {
                    if x {
                        "1"
                    } else {
                        "0"
                    }
                })
                .collect_vec()
                .join(""),
            amplitude
        );
        Ok(())
    }
}
