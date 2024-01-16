use super::{
    iqp_circuit::build_iqp_circuit, simulation::Simulation, simulation_params::SimulationParams,
    statevector::generate_random_statevector,
};
use crate::{
    gray_code_flip_bit::GrayCodeFlipBit,
    iqp_simulations::{iqp_circuit::QubitColoringIndexes, linear_systems::LinearSystems},
};
use anyhow::{Context, Result};
use bitvec::vec::BitVec;
use itertools::Itertools;
use std::time::Instant;

pub struct CPUSmallIntSimulation {
    params: SimulationParams,
    statevector: Option<BitVec>
}

impl CPUSmallIntSimulation {
    pub fn new(boolean_cube_dimension: u32) -> Self {
        Self {
            params: SimulationParams::new(boolean_cube_dimension),
            statevector: None
        }
    }

    pub fn with_sv(boolean_cube_dimension: u32, sv: BitVec) -> Self {
        Self {
            params: SimulationParams::new(boolean_cube_dimension),
            statevector: Some(sv)
        }
    }
}

impl Simulation for CPUSmallIntSimulation {
    fn run(&mut self) -> Result<f64> {
        let start = Instant::now();
        let s: BitVec = self.statevector.take().unwrap_or_else(|| {
            generate_random_statevector(&self.params)
        });
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
            if let Some(amplitude_increment) =
                ls.solve_if_gamma_null_space_quick_check(&s_b, &s_g, &s_r)
            {
                amplitude += amplitude_increment;
            }
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
        Ok(amplitude)
    }
}

#[cfg(test)]
mod test {
    use bitvec::vec::BitVec;
    use bitvec::prelude::*;
    use super::*;


    fn run_amplitude_baseline(statevector: BitVec, expected_amplitude: f64) {
        const TOLERANCE: f64 = 0.000001;
        let mut sim = CPUSmallIntSimulation::with_sv(2, statevector);
        let res = sim.run().expect("simulation ran with no error");
        let difference = res - expected_amplitude;
        assert!(difference.abs() < TOLERANCE, "GOT: {res} | EXPECTED: {expected_amplitude}");
    }

    #[test]
    fn test_amplitudes() {
        let testcases = vec![
            (bitvec![1,0,0,1,0,0,0,0,0,1,1,0], -0.019531249999999986),
            (bitvec![0,1,1,0,0,1,0,1,1,1,0,0], -0.003906249999999994),
            (bitvec![0,1,0,0,1,1,0,1,0,1,0,0], 0.011718749999999991),
            (bitvec![1,1,0,1,0,1,0,1,0,0,0,1], 0.027343749999999983),
            (bitvec![0,0,1,1,0,1,1,1,1,0,0,0], 0.0039062499999999952),
            (bitvec![0,1,1,0,0,1,0,0,1,0,0,1], 0.027343749999999983),
        ];

        for (sv, expected_amplutde) in testcases {
            run_amplitude_baseline(sv, expected_amplutde)
        }
    }

}