use bitvec::prelude::*;

use super::{gauss_jordan::GaussJordan, matrix::BitMatrix};

pub struct BackwardsSubstitution {
    pub solution: BitVec,
    pub gj: GaussJordan,
}

impl BackwardsSubstitution {
    pub fn zero(n: usize) -> Self {
        Self {
            solution: bitvec![usize, Lsb0; 0; n],
            gj: GaussJordan::zero(n, n + 1)
        }
    }

    pub fn solve(&mut self, u: &BitMatrix, b: &BitVec) -> Option<()> {
        // TODO: make augmented system for gj and extract
        self.solution.set_elements(0x0);
        let n = self.solution.len();
        let mut i = n;
        self.gj.copy_from_augmented_system(u, b);
        self.gj.go_to_echelon_form();
        while i > 0 {
            i -= 1;
            let mut tmp = self.gj.rows[i][self.gj.last_col_idx];
            let row = &self.gj.rows[i];
            if i >= self.gj.rank {
                if tmp {
                    // zero row and b[j] != 0 means there's no
                    // solution to Ux = b
                    return None;
                }
                unsafe {
                    self.solution.set_unchecked(i, false);
                }
                continue;
            }
            let diag = row[i];
            if !diag {
                // arbitrarily clear the solution at index i if it is a free
                // variable
                unsafe {
                    self.solution.set_unchecked(i, false);
                }
                continue;
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
        solver.solve(&matrix, &bits![1, 1, 0, 0].to_bitvec());
        assert_eq!(&solver.solution, &bits![0, 1, 0, 0].to_bitvec());
    }

    #[test]
    fn test_b2_input_2() {
        let mut solver = BackwardsSubstitution::zero(4);
        // [1, 0, 1, 1]
        // [0, 1, 0, 1]
        // [1, 0, 0, 1]
        // [1, 1, 1, 0]
        let rows: Vec<Vec<usize>> = vec![vec![0, 2, 3], vec![1, 3], vec![0, 3], vec![0, 1, 2]];
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let mut gj = GaussJordan::zero(4, 4);
        gj.copy_from_matrix(&matrix);
        gj.go_to_echelon_form();
        assert_eq!(
            gj.rows.iter().map(|bv| bv.to_string()).collect_vec(),
            vec![
                "[1, 0, 1, 1]",
                "[0, 1, 0, 1]",
                "[0, 0, 1, 0]",
                "[0, 0, 0, 0]",
            ]
        );
        let b = bits![0, 1, 0, 1].to_bitvec();
        let x = bits![0, 1, 0, 0].to_bitvec();
        solver.solve(&matrix, &b);
        assert_eq!(&solver.solution, &x);
    }
}
