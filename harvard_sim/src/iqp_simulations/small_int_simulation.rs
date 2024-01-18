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
use std::{time::Instant, collections::{HashMap, BTreeMap, BTreeSet, HashSet}};

pub struct CPUSmallIntSimulation {
    params: SimulationParams,
    statevector: Option<BitVec>,
}

impl CPUSmallIntSimulation {
    pub fn new(boolean_cube_dimension: u32) -> Self {
        Self {
            params: SimulationParams::new(boolean_cube_dimension),
            statevector: None,
        }
    }

    pub fn with_sv(boolean_cube_dimension: u32, sv: BitVec) -> Self {
        Self {
            params: SimulationParams::new(boolean_cube_dimension),
            statevector: Some(sv),
        }
    }
}

impl Simulation for CPUSmallIntSimulation {
    fn run(&mut self) -> Result<f64> {
        let start = Instant::now();
        let s: BitVec = self
            .statevector
            .take()
            .unwrap_or_else(|| generate_random_statevector(&self.params));
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
                // println!("Updated amplitude: [{}] {} + {} = {}",ls.symmetry_checker.as_ref().unwrap().current_bit_pattern,amplitude - amplitude_increment, amplitude_increment, amplitude);
            }
            ls.update_with_flip_bit(flip_bit, &phase_graph);
        }
        amplitude /= (1 << self.params.nodes) as f64;
        let end = Instant::now();
        log::debug!("Time to execute: {:#?}", end - start);
        let mut increments = ls.increments.clone();
        let missed_increments: BTreeMap<u16, f64> = ls.dmitri_skips.iter().filter_map(|skip| {
            increments.remove_entry(skip)
        }).collect();
        let missed_symmetry_groups: BTreeMap<BTreeSet<u16>, f64> = missed_increments.iter().map(|(skipped_key, _), | {
            (ls.symmetry_groups.get(&skipped_key).unwrap().clone(), 0.0_f64)
        }).collect();
        let mut actual_group_accumulation: BTreeMap<BTreeSet<u16>, Vec<f64>> = missed_symmetry_groups.keys().map(|group| {
            let mut acc = vec![];
            for member in group.iter() {
                 acc.push(ls.increments.get(member).copied().unwrap_or(0.0));
            }
            (group.clone(), acc)
        }).collect();
        let dmitri_check_group_accumulation: BTreeMap<BTreeSet<u16>, Vec<f64>> = missed_symmetry_groups.keys().map(|group| {
            let first_val = ls.increments.get(group.first().unwrap()).copied().unwrap_or(0.0);
            let mut acc = vec![];
            for _ in group.iter() {
                 acc.push(first_val);
            }
            (group.clone(), acc)
        }).collect();
        let error_terms: BTreeMap<BTreeSet<u16>, (Vec<f64>, Vec<f64>)> = actual_group_accumulation.iter().filter_map(|(group, acc)| {
            let mut acc_summed = 0.0;
            for val in acc {
                acc_summed += val;
            }
            let mut dmitri_sum = 0.0;
            let dmitri_val = dmitri_check_group_accumulation.get(group).unwrap();
            for val in dmitri_val {
                dmitri_sum += val;
            }
            if (dmitri_sum - acc_summed).abs() < 0.000001 {
                return None
            }
            Some((group.clone(), (acc.clone(), dmitri_val.clone())))
        }).collect();
        let good_terms: BTreeMap<BTreeSet<u16>, (Vec<f64>, Vec<f64>)> = actual_group_accumulation.iter().filter_map(|(group, acc)| {
            let mut acc_summed = 0.0;
            for val in acc {
                acc_summed += val;
            }
            let mut dmitri_sum = 0.0;
            let dmitri_val = dmitri_check_group_accumulation.get(group).unwrap();
            for val in dmitri_val {
                dmitri_sum += val;
            }
            for (idx, val) in acc.iter().enumerate() {
                let dmitri_val = dmitri_val[idx];
                if (val - dmitri_val).abs() > 0.000001 {
                    return None
                }
            }
            Some((group.clone(), (acc.clone(), dmitri_val.clone())))
        }).collect();
        println!("total missed terms: {}", actual_group_accumulation.len());
        println!("erroneously dropped terms: {}", error_terms.len());
        println!("good terms: {}", good_terms.len());
        for error_term in error_terms {
            println!("{:?}", error_term);
        }
        let all_symmetry_groups: BTreeSet<BTreeSet<u16>> = ls.symmetry_groups.values().cloned().collect();
        let symmetry_group_results: BTreeSet<BTreeMap<u16, String>> = all_symmetry_groups.into_iter().map(|group| {
            group.into_iter().map(|member| {
                (member, ls.increments.get(&member).copied().unwrap_or(0.0).to_string())
            }).collect()
        }).collect();
        println!("symmetry results");
        for symmetry_result in symmetry_group_results {
            println!("{:?}", symmetry_result);
        }
        let had_results: BTreeSet<(u16, String)> = ls.increments.iter().filter_map(|(key, res)| {
            if (res).abs() > 0.0 {Some(( *key, res.to_string() ))} else { None }
        }).collect();
        println!("actually had {} results!", had_results.len());
        for (key, res) in had_results {
            println!("0b{:b}: {:?}", key, res);
        }
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
    use super::*;
    use bitvec::prelude::*;
    use bitvec::vec::BitVec;

    fn run_amplitude_baseline(statevector: BitVec, expected_amplitude: f64) {
        const TOLERANCE: f64 = 1e-8;
        let dimension = (statevector.len() / 3).ilog2();
        let mut sim = CPUSmallIntSimulation::with_sv(dimension, statevector);
        let res = sim.run().expect("simulation ran with no error");
        let difference = res - expected_amplitude;
        assert!(
            difference.abs() < TOLERANCE,
            "GOT: {res} | EXPECTED: {expected_amplitude}"
        );
    }

    #[test]
    fn test_amplitudes() {
        let testcases = vec![
            (
                bitvec![1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0],
                -0.019531249999999986,
            ),
            (
                bitvec![0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0],
                -0.003906249999999994,
            ),
            (
                bitvec![0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0],
                0.011718749999999991,
            ),
            (
                bitvec![1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1],
                0.027343749999999983,
            ),
            (
                bitvec![0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0],
                0.0039062499999999952,
            ),
            (
                bitvec![0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1],
                0.027343749999999983,
            ),
        ];

        for (sv, expected_amplutde) in testcases {
            run_amplitude_baseline(sv, expected_amplutde)
        }
    }

    #[test]
    fn test_amplitudes_b_3 () {
        let testcases = vec![
            (
                bitvec![1,1,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,1,1],
                -0.0001220703125
            )
        ];
        for (sv, expected_amplutde) in testcases {
            run_amplitude_baseline(sv, expected_amplutde)
        }
    }

    #[test]
    fn test_amplitudes_b_4() { 
        let testcases = vec![
            (
                bitvec![1,0,0,1,0,0,1,1,1,0,1,0,0,1,0,0,0,1,1,0,0,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,0,1,0,0,1,1,1],
                3.003_515_303_134_918e-8,
            ),
            (
                bitvec![1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0],
                1.606_531_441_211_700_4e-8
            )
        ];
        for (sv, expected_amplutde) in testcases {
            run_amplitude_baseline(sv, expected_amplutde)
        }
    }
}
