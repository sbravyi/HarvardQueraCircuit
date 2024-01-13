use crate::{
    phase_polynomial::PhasePolynomial,
    qubit::{Color, Qubit},
};
use anyhow::Result;

use super::simulation_params::SimulationParams;

pub struct QubitColoring {
    pub red: Vec<Qubit>,
    pub blue: Vec<Qubit>,
    pub green: Vec<Qubit>,
}

impl QubitColoring {
    pub fn new_for_n_qubits(n_qubits: u32) -> Self {
        let mut coloring = Self {
            red: vec![],
            blue: vec![],
            green: vec![],
        };
        for i in 0..n_qubits {
            let qubit = Qubit::new(i);
            match qubit.color {
                Color::Red => coloring.red.push(qubit),
                Color::Blue => coloring.blue.push(qubit),
                Color::Green => coloring.green.push(qubit),
            }
        }
        coloring
    }
}

pub fn build_iqp_circuit(params: &SimulationParams) -> Result<(PhasePolynomial, QubitColoring)> {
    let c = QubitColoring::new_for_n_qubits(params.n_qubits);
    let mut pp = PhasePolynomial::new();
    // apply the initial layer of "A-rectangles", see page 29 in
    // https://arxiv.org/pdf/2312.03982.pdf
    for i in 0..params.nodes as usize {
        pp.ccz(c.red[i], c.blue[i], c.green[i])?;
        pp.cz(c.red[i], c.blue[i])?;
        pp.cz(c.blue[i], c.green[i])?;
        pp.cz(c.red[i], c.green[i])?;
        // we ignore pauli Z gates since they can be absorbed into a Pauli frame
    }
    for direction in 0..params.boolean_cube_dimension {
        // apply CNOTs oriented along this direction on the cube
        // cube nodes with even pariry = control qubits
        // cube nodes with odd parity = target qubits
        for x in 0..params.nodes as usize {
            if x.count_ones() % 2 == 0 {
                continue;
            }
            let y = x ^ (1 << direction) as usize;
            pp.cnot(c.red[x], c.red[y])?;
            pp.cnot(c.blue[x], c.blue[y])?;
            pp.cnot(c.green[x], c.green[y])?;
        }
        // alternate between layers of A or B rectangles, see page 29 in
        // https://arxiv.org/pdf/2312.03982.pdf
        // some A/B rectangles acting on nodes with even parity cancel each other
        for i in 0..params.nodes as usize {
            pp.ccz(c.red[i], c.blue[i], c.green[i])?;
            pp.cz(c.red[i], c.blue[i])?;
            pp.cz(c.blue[i], c.green[i])?;
            if direction % 2 != 0 {
                pp.cz(c.red[i], c.green[i])?;
            }
        }
    }
    Ok((pp, c))
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use itertools::Itertools;

    use super::*;

    fn indexes(qbs: Vec<Qubit>) -> Vec<u32> {
        qbs.into_iter().map(|q| q.index).collect_vec()
    }

    #[test]
    fn test_iqp_construction() {
        let params = SimulationParams::new(2);
        let (pp, c) = build_iqp_circuit(&params).unwrap();
        let pg = pp.into_polynomial_graph().unwrap();
        assert_eq!(indexes(c.red), vec![0, 3, 6, 9]);
        assert_eq!(indexes(c.blue), vec![1, 4, 7, 10]);
        assert_eq!(indexes(c.green), vec![2, 5, 8, 11]);
        assert_eq!(
            pg.rbg_monomials,
            HashMap::<u32, Vec<(u32, u32)>>::from_iter(vec![
                (
                    0,
                    vec![(3, 1), (1, 3), (1, 0), (0, 1), (3, 2), (2, 3), (1, 1),]
                ),
                (1, vec![(0, 1), (0, 3), (3, 0), (0, 0), (1, 0), (1, 1)]),
                (2, vec![(3, 2), (2, 3), (2, 2), (0, 3), (3, 0), (3, 3)]),
                (
                    3,
                    vec![(0, 1), (3, 2), (2, 2), (2, 3), (2, 0), (1, 0), (0, 2)]
                )
            ])
        );
        assert_eq!(
            pg.rg_monomials,
            HashMap::<u32, Vec<u32>>::from_iter(vec![
                (0, vec![1, 2]),
                (1, vec![3, 0]),
                (2, vec![3, 0]),
                (3, vec![2, 1]),
            ])
        );
        assert_eq!(
            pg.rb_monomials,
            HashMap::<u32, Vec<u32>>::from_iter(vec![
                (0, vec![1]),
                (1, vec![1, 0]),
                (2, vec![2, 3]),
                (3, vec![2]),
            ])
        );
        assert_eq!(
            pg.bg_monomials,
            vec![(1, 1), (1, 0), (0, 1), (3, 2), (2, 3), (2, 2),]
        );
    }
}
