use bitvec::prelude::*;

use super::gauss_jordan::GaussJordan;

pub struct BackwardsSubstitution {
    pub solution: BitVec,
}

impl BackwardsSubstitution {
    pub fn zero(n: usize) -> Self {
        Self {
            solution: bitvec![usize, Lsb0; 0; n],
        }
    }

    pub fn solve(&mut self, u: &GaussJordan, b: &BitVec) -> Option<()> {
        self.solution.set_elements(0x0);
        let n = self.solution.len();
        let mut i = n;
        while i > 0 {
            i -= 1;
            let mut tmp = b[i];
            let row = &u.rows[i];
            if i >= u.rank {
                if tmp {
                    // zero row and b[j] != 0 means there's no
                    // solution to Ux = b
                    return None;
                }
                unsafe {
                    self.solution.set_unchecked(i, false);
                }
            }
            let diag = row[i];
            if !diag {
                // arbitrarily clear the solution at index i if it is a free
                // variable
                unsafe {
                    self.solution.set_unchecked(i, false);
                }
            }
            for j in i..n {
                tmp ^= row[j] & self.solution[j]
            }
            unsafe {
                self.solution.set_unchecked(i, tmp);
            }
        }
        Some(())
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use crate::bit_matrix::matrix::BitMatrix;

    use super::*;

    #[test]
    fn test_b2_input() {
        let mut solver = BackwardsSubstitution::zero(4);
        // [0, 1, 0, 1]
        // [1, 1, 0, 0]
        // [0, 0, 0, 0]
        // [1, 0, 0, 1]
        let rows: Vec<Vec<usize>> = vec![vec![1, 3], vec![0, 1], vec![], vec![0, 3]];
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let mut gj = GaussJordan::zero(4, 4);
        gj.copy_from_matrix(&matrix);
        gj.go_to_echelon_form();
        solver.solve(&gj, &bits![1, 1, 0, 0].to_bitvec());
        assert_eq!(&solver.solution, &bits![0, 1, 0, 0].to_bitvec());
    }
}
