use anyhow::{ensure, Result};
use indexmap::{IndexMap, IndexSet};
use itertools::Itertools;
use std::collections::{HashMap, HashSet};

use crate::qubit::{Color, Qubit};

#[derive(Hash, PartialEq, Eq, Debug, Clone)]
enum Monomial {
    Two([Qubit; 2]),
    Three([Qubit; 3]),
}

impl Monomial {
    fn contains(&self, q: &Qubit) -> bool {
        match self {
            Self::Two([q1, q2]) => q1 == q || q2 == q,
            Self::Three([q1, q2, q3]) => q1 == q || q2 == q || q3 == q,
        }
    }

    fn swap_qubit(&mut self, original_qubit: &Qubit, new_qubit: Qubit) {
        match self {
            Self::Two(tup) => {
                if &tup[0] == original_qubit {
                    tup[0] = new_qubit
                } else {
                    tup[1] = new_qubit
                }
                tup.sort()
            }
            Self::Three(tup) => {
                if &tup[0] == original_qubit {
                    tup[0] = new_qubit
                } else if &tup[1] == original_qubit {
                    tup[1] = new_qubit
                } else {
                    tup[2] = new_qubit
                }
                tup.sort()
            }
        }
    }
}

struct MonomialIterator {
    index: usize,
    monomial: Monomial,
}

impl Iterator for MonomialIterator {
    type Item = Qubit;

    fn next(&mut self) -> Option<Self::Item> {
        let old_index = self.index;
        self.index += 1;
        match self.monomial {
            Monomial::Two(tup) => {
                if old_index > 1 {
                    return None;
                }
                return Some(tup[old_index]);
            }
            Monomial::Three(tup) => {
                if old_index > 2 {
                    return None;
                }
                return Some(tup[old_index]);
            }
        }
    }
}

impl IntoIterator for Monomial {
    type Item = Qubit;
    type IntoIter = MonomialIterator;

    fn into_iter(self) -> MonomialIterator {
        MonomialIterator {
            monomial: self,
            index: 0,
        }
    }
}
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
    monomials: IndexSet<Monomial>,
}

impl PhasePolynomial {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn ccz(&mut self, q1: Qubit, q2: Qubit, q3: Qubit) -> Result<()> {
        ensure!(HashSet::from([q1.color, q2.color, q3.color]).len() == 3);
        let mut qubits = [q1, q2, q3];
        qubits.sort();
        let qubits = Monomial::Three(qubits);
        if self.monomials.contains(&qubits) {
            self.monomials.shift_remove(&qubits);
        } else {
            self.monomials.insert(qubits);
        }
        Ok(())
    }

    pub fn cz(&mut self, q1: Qubit, q2: Qubit) -> Result<()> {
        ensure!(q1.color != q2.color);
        let mut qubits = [q1, q2];
        qubits.sort();
        let qubits = Monomial::Two(qubits);
        if self.monomials.contains(&qubits) {
            self.monomials.shift_remove(&qubits);
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
            controlled_m.swap_qubit(&t, c);
            if self.monomials.contains(&controlled_m) {
                self.monomials.shift_remove(&controlled_m);
            } else {
                self.monomials.insert(controlled_m);
            }
        }
        Ok(())
    }

    pub fn into_polynomial_graph(self) -> Result<PolynomialGraph> {
        let mut rbg_monomials: IndexMap<u32, Vec<(u32, u32)>> = IndexMap::new();
        let mut rg_monomials: IndexMap<u32, Vec<u32>> = IndexMap::new();
        let mut rb_monomials: IndexMap<u32, Vec<u32>> = IndexMap::new();
        let mut bg_monomials: Vec<(u32, u32)> = Vec::new();
        for monomial in self.monomials {
            let monomial_colors: HashMap<Color, u32> = Color::seperate_monomial_colors(monomial);
            if !monomial_colors.contains_key(&Color::Blue) {
                rg_monomials
                    .entry(*monomial_colors.get(&Color::Red).unwrap())
                    .or_default()
                    .push(*monomial_colors.get(&Color::Green).unwrap());
            } else if !monomial_colors.contains_key(&Color::Green) {
                rb_monomials
                    .entry(*monomial_colors.get(&Color::Red).unwrap())
                    .or_default()
                    .push(*monomial_colors.get(&Color::Blue).unwrap());
            } else if !monomial_colors.contains_key(&Color::Red) {
                bg_monomials.push((
                    *monomial_colors.get(&Color::Blue).unwrap(),
                    *monomial_colors.get(&Color::Green).unwrap(),
                ))
            } else {
                // no missing colors - this is an rbg monomial
                let targeted_qubits = [
                    *monomial_colors.get(&Color::Blue).unwrap(),
                    *monomial_colors.get(&Color::Green).unwrap(),
                ];
                rbg_monomials
                    .entry(*monomial_colors.get(&Color::Red).unwrap())
                    .or_default()
                    .push((targeted_qubits[0], targeted_qubits[1]));
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
    pub rbg_monomials: IndexMap<u32, Vec<(u32, u32)>>,
    pub rb_monomials: IndexMap<u32, Vec<u32>>,
    pub rg_monomials: IndexMap<u32, Vec<u32>>,
    pub bg_monomials: Vec<(u32, u32)>,
}

#[cfg(test)]
mod tests {
    use indexmap::IndexMap;

    use super::*;

    #[test]
    fn test_phase_graph() {
        let mut pp = PhasePolynomial::new();
        pp.ccz(Qubit::new(0), Qubit::new(5), Qubit::new(10)).unwrap();
        let pg = pp.into_polynomial_graph().unwrap();
        assert_eq!(
            pg.rbg_monomials,
            IndexMap::<u32, Vec<(u32, u32)>>::from_iter(vec![(0, vec![(3, 1)])])
        )
    }
}
