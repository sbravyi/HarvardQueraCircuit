use super::{
    iqp_circuit::{build_iqp_circuit, QubitColoring},
    simulation::Simulation,
    simulation_params::SimulationParams,
    statevector::{generate_random_statevector, LinearSystems},
};
use crate::gray_code_flip_bit::GrayCodeFlipBit;
use anyhow::{Context, Result};
use itertools::Itertools;
use std::time::Instant;

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
        let s = generate_random_statevector(&self.params);
        log::debug!(
            "Statevector: <{}|",
            s.iter_dense()
                .map(|c| if c.is_one() { "1" } else { "0" })
                .collect_vec()
                .join("")
        );
        let start = Instant::now();
        let (phase_polynomial, coloring) =
            build_iqp_circuit(&self.params).context("Building IQP circuit")?;
        let QubitColoring { red, green, blue } = coloring;
        let phase_graph = phase_polynomial
            .into_polynomial_graph()
            .context("converting phase polynomial into graph")?;
        let gc_flip_bit = GrayCodeFlipBit::new(self.params.nodes);
        let s_r = s
            .keep_only_positions(&red.into_iter().map(|q| q.index as usize).collect_vec())
            .context("collecting red components of statevector")?;
        let s_b = s
            .keep_only_positions(&blue.into_iter().map(|q| q.index as usize).collect_vec())
            .context("collecting blue components of statevector")?;
        let s_g = s
            .keep_only_positions(&green.into_iter().map(|q| q.index as usize).collect_vec())
            .context("collecting green components of statevector")?;
        let mut ls = LinearSystems::new(&self.params, &phase_graph);
        let mut amplitude: f64 = 0.0;
        for flip_bit in gc_flip_bit {
            if ls.is_gamma_null_space_quick_check(&s_b, &s_g)? {}
            ls = ls
                .update_with_flip_bit(flip_bit, &phase_graph)
                .context("updating flip bit")?;
        }
        let end = Instant::now();
        log::debug!("Time to execute: {:#?}", end - start);
        Ok(())
    }
}
