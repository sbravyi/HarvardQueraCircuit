use anyhow::{ensure, Result};
use itertools::Itertools;
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::qubit::{Color, Qubit};

// ### Functions to simulate CCZ/CZ/CNOT gates via the phase polynomial formalism
// # We only consider phase polynomials with degree-2 and degree-3 terms
// # A phase polynomial is described by a list M of tuples of variable indices
// # Example:
// # f(x_0,x_1,x_2,x_3) = x_0*x_1*x_3 + x_1*x_2 (mod 2)
// # is described by M = [(0,1,3),(1,2)]
// # A phase polynomial M can be converted to n-qubit state |psi> as follows:
// #   N = 1<<n
// #	psi = np.ones(N)
// #	for i in range(N):
// #		x = to_binary(n,i)
// #		for h in M:
// #			if all([x[i] for i in h]):
// #				psi[i]*=-1
// #	psi = psi/np.sqrt(N)
#[derive(Default)]
pub struct PhasePolynomial {
    monomials: HashSet<BTreeSet<Qubit>>,
}

impl PhasePolynomial {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn ccz(&mut self, q1: Qubit, q2: Qubit, q3: Qubit) -> Result<()> {
        ensure!(HashSet::from([q1.color, q2.color, q3.color]).len() == 3);
        let qubits = BTreeSet::from([q1, q2, q3]);
        if self.monomials.contains(&qubits) {
            self.monomials.remove(&qubits);
        } else {
            self.monomials.insert(qubits);
        }
        Ok(())
    }

    pub fn cz(&mut self, q1: Qubit, q2: Qubit) -> Result<()> {
        ensure!(q1.color != q2.color);
        let qubits = BTreeSet::from([q1, q2]);
        if self.monomials.contains(&qubits) {
            self.monomials.remove(&qubits);
        } else {
            self.monomials.insert(qubits);
        }
        Ok(())
    }

    pub fn cnot(&mut self, c: Qubit, t: Qubit) -> Result<()> {
        ensure!(c.index != t.index);
        ensure!(c.color == t.color);
        let target_impacted_monomials = self
            .monomials
            .iter()
            .filter(|&m| m.contains(&t))
            .cloned()
            .collect_vec();

        for m in target_impacted_monomials {
            let mut controlled_m = m.clone();
            controlled_m.remove(&t);
            controlled_m.insert(c);
            if self.monomials.contains(&controlled_m) {
                self.monomials.remove(&controlled_m);
            } else {
                self.monomials.insert(controlled_m);
            }
        }
        Ok(())
    }

    pub fn into_polynomial_graph(self) -> Result<PolynomialGraph> {
        let mut rbg_monomials: HashMap<u32, Vec<(u32, u32)>> = HashMap::new();
        let mut rg_monomials: HashMap<u32, Vec<u32>> = HashMap::new();
        let mut rb_monomials: HashMap<u32, Vec<u32>> = HashMap::new();
        let mut bg_monomials: Vec<(u32, u32)> = Vec::new();
        for monomial in self.monomials {
            let monomial_colors: HashMap<Color, Qubit> =
                monomial.into_iter().map(|q| (q.color, q)).collect();
            if !monomial_colors.contains_key(&Color::Blue) {
                rg_monomials
                    .entry(monomial_colors.get(&Color::Red).unwrap().index)
                    .or_default()
                    .push(monomial_colors.get(&Color::Green).unwrap().index);
            } else if !monomial_colors.contains_key(&Color::Green) {
                rb_monomials
                    .entry(monomial_colors.get(&Color::Red).unwrap().index)
                    .or_default()
                    .push(monomial_colors.get(&Color::Blue).unwrap().index);
            } else if !monomial_colors.contains_key(&Color::Red) {
                bg_monomials.push((
                    monomial_colors.get(&Color::Blue).unwrap().index,
                    monomial_colors.get(&Color::Green).unwrap().index,
                ))
            } else {
                // no missing colors - this is an rbg monomial
                let mut targeted_qubits = [
                    monomial_colors.get(&Color::Blue).unwrap(),
                    monomial_colors.get(&Color::Green).unwrap(),
                ];
                targeted_qubits.sort();
                rbg_monomials
                    .entry(monomial_colors.get(&Color::Red).unwrap().index)
                    .or_default()
                    .push((targeted_qubits[0].index, targeted_qubits[1].index));
            }
        }

        Ok(PolynomialGraph {
            rbg_monomials,
            rg_monomials,
            rb_monomials,
            bg_monomials,
        })
    }
}

pub struct PolynomialGraph {
    pub rbg_monomials: HashMap<u32, Vec<(u32, u32)>>,
    pub rb_monomials: HashMap<u32, Vec<u32>>,
    pub rg_monomials: HashMap<u32, Vec<u32>>,
    pub bg_monomials: Vec<(u32, u32)>,
}
