use bitvec::prelude::*;
use bitvec::{field::BitField, vec::BitVec};
use itertools::Itertools;

use super::column_matrix::ColumnMatrix;
use super::matrix::BitMatrix;

pub struct SergeySolver {
    augmented_system: BitMatrix,
    n: usize,
    x: ColumnMatrix,
    pub solution: BitVec,
    syndrome: BitVec,
    zero_b: bool,
}

fn debug_bitmatrix(m: &BitMatrix) {
    println!("{}", m.rows.iter().map(|r| { r.to_string() }).join("\n"));
}

fn debug_col_matrix(m: &ColumnMatrix) {
    let mut transposed = BitMatrix::zeroes(m.number_of_rows, m.cols.len());
    for (col_idx, col) in m.cols.iter().enumerate() {
        for (row_idx, bit) in col.iter().enumerate() {
            transposed.rows[row_idx].set(col_idx, *bit);
        }
    }
    debug_bitmatrix(&transposed);
}

fn debug_bitvec(b: &BitVec) {
    println!("{}", b);
}

impl SergeySolver {
    pub fn zero(n: usize) -> Self {
        Self {
            augmented_system: BitMatrix::zeroes(n, n + 1),
            n,
            x: ColumnMatrix::zeroes(n + 1, n + 1),
            solution: bitvec![usize, Lsb0; 0; n],
            syndrome: bitvec![usize, Lsb0; 0; n + 1],
            zero_b: false,
        }
    }

    pub fn debug(&self) {
        println!("a");
        debug_bitmatrix(&self.augmented_system);
        println!("x");
        debug_col_matrix(&self.x);
        println!("sol");
        debug_bitvec(&self.solution);
        println!("syndrome");
        debug_bitvec(&self.syndrome);
    }

    pub fn copy_from_augmented_system(&mut self, m: &BitMatrix, b: &BitVec) {
        debug_assert_eq!(m.number_of_columns, self.n);
        debug_assert_eq!(b.len(), self.n);
        debug_assert_eq!(m.rows.len(), self.augmented_system.rows.len());
        self.zero_b = b.first_one().is_none();
        let x_rank = if self.zero_b { self.n } else { self.n + 1 };
        for (idx, row) in m.rows.iter().enumerate() {
            self.augmented_system.rows[idx][..self.n].copy_from_bitslice(&row[..]);
        }
        for (idx, b) in b.iter().enumerate() {
            unsafe {
                self.augmented_system.rows[idx].set_unchecked(self.n, *b);
            }
        }
        self.syndrome.store(0x0);
        self.x.reset();
        if self.zero_b {
            self.x.remove_col(self.n);
        }
        // quick identity matrix creation
        for idx in 0..x_rank {
            self.x.cols[idx].store(1 << idx);
        }
    }

    fn reformulate_x_from_augmented_system(&mut self, current_equation: usize) -> Option<()> {
        let b_row = &self.augmented_system.rows[current_equation];
        for (col_idx, col) in self.x.cols.iter().enumerate() {
            let mut syndrome_val = false;
            // syndrome val is the inner product of each x column with the current b row
            for (b_one_idx, bit) in b_row.iter().enumerate() {
                syndrome_val ^= *bit & col[b_one_idx];
            }
            unsafe {
                self.syndrome.set_unchecked(col_idx, syndrome_val);
            }
        }
        let n_variables = self.x.cols.len();
        let mut bad_cols = self
            .syndrome
            .iter_ones()
            .filter(|col_idx| *col_idx < n_variables);
        let first_bad_col_idx = bad_cols.next();
        // TODO: maybe try iterating in reverse, take em all out at once, and then
        // add them while they're in removed
        if let Some(first_bad_col_idx) = first_bad_col_idx {
            self.x.remove_col(first_bad_col_idx);
            let first_bad_col = self.x.pop_from_removed();
            // add the first bad column to all other bad columns
            // then take the remaining bad columns and project them
            // onto x
            let mut bad_cols_removed = 1;
            for bad_col_idx in bad_cols {
                self.x.remove_col(bad_col_idx - bad_cols_removed);
                let mut bad_col = self.x.pop_from_removed();
                bad_col[..] ^= &first_bad_col;
                self.x.push_back_into_usage(bad_col);
                bad_cols_removed += 1;
            }
            self.x.put_back_in_removed(first_bad_col);
        }
        Some(())
    }

    fn take_arbitrary_solution_from_x(&mut self) -> Option<()> {
        let last_row = self.x.number_of_rows - 1;
        if self.zero_b {
            return self.x.cols.first().map(|first_solution| {
                self.solution[..self.n].copy_from_bitslice(&first_solution[..self.n]);
            });
        }
        let sol_col_idx = if let Some((col_idx, solution)) = self
            .x
            .cols
            .iter()
            .enumerate()
            .find(|(_, col)| col[last_row])
        {
            self.solution[..self.n].copy_from_bitslice(&solution[..self.n]);
            Some(col_idx)
        } else {
            None
        };
        sol_col_idx.map(|sol_col_idx| self.calculate_nullspace_from_solution_column(sol_col_idx))
    }

    fn calculate_nullspace_from_solution_column(&mut self, col_idx: usize) {
        if self.zero_b {
            return;
        }
        self.x.remove_col(col_idx);
        let sol_col = self.x.pop_from_removed();
        if sol_col[self.n] {
            for col in self.x.cols.iter_mut().filter(|col| col[self.n]) {
                col[..self.n] ^= &sol_col;
            }
        }
        self.x.put_back_in_removed(sol_col);
    }

    pub fn rank(&self) -> usize {
        self.n.saturating_sub(self.x.cols.len())
    }

    pub fn is_full_rank(&self) -> bool {
        self.x.cols.is_empty()
    }

    fn full_rank_result(&mut self) -> Option<()> {
        if self.zero_b {
            None
        } else {
            self.solution.store(0x0);
            Some(())
        }
    }

    pub fn solve(&mut self, u: &BitMatrix, b: &BitVec) -> Option<()> {
        self.copy_from_augmented_system(u, b);
        for i in 0..self.n {
            self.reformulate_x_from_augmented_system(i)?;
            if self.is_full_rank() {
                return self.full_rank_result();
            }
        }
        self.take_arbitrary_solution_from_x()?;
        Some(())
    }

    pub fn is_nullspace_codeword(&self, codeword: &BitVec) -> bool {
        for col in self.x.cols.iter() {
            let mut inner_product = false;
            for (bit_idx, bit) in codeword
                .iter()
                .enumerate()
                .filter(|(bit_idx, _)| *bit_idx < self.n)
            {
                let col_val = *bit & col[bit_idx];
                inner_product ^= col_val
            }
            if inner_product {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::bit_matrix::sergey_solver::SergeySolver;

    fn four_by_four_system(mat_rows: Vec<Vec<usize>>, b: BitVec, x: BitVec) {
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in mat_rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let mut solver = SergeySolver::zero(4);
        solver.solve(&matrix, &b).unwrap();
        assert_eq!(solver.solution.to_string(), x.to_string());
    }

    fn four_by_four_system_failure(mat_rows: Vec<Vec<usize>>, b: BitVec) {
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in mat_rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let mut solver = SergeySolver::zero(4);
        assert!(solver.solve(&matrix, &b).is_none());
    }

    fn nullspace_check(
        mat_rows: Vec<Vec<usize>>,
        b: BitVec,
        x: BitVec,
        codeword: BitVec,
        expect_syndrome: bool,
    ) {
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in mat_rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let mut solver = SergeySolver::zero(4);
        solver.solve(&matrix, &b).unwrap();
        assert_eq!(
            solver.solution.to_string(),
            x.to_string(),
            "solution mismatch"
        );
        assert_eq!(
            solver.is_nullspace_codeword(&codeword),
            expect_syndrome,
            "codeword expectation mismatch"
        );
    }

    #[test]
    fn test_flip_code_0() {
        // gamma
        // [0, 1, 0, 0]
        // [1, 1, 0, 0]
        // [0, 0, 1, 1]
        // [0, 0, 1, 0]
        // b
        // [1, 1, 1, 1]
        // x
        // [1, 0, 1, 0]
        four_by_four_system(
            vec![vec![1], vec![0, 1], vec![2, 3], vec![2]],
            bitvec![1, 1, 1, 1].to_bitvec(),
            bitvec![0, 1, 1, 0].to_bitvec(),
        )
    }

    #[test]
    fn test_input_3() {
        // gamma
        // [1, 1, 0, 0]
        // [1, 1, 0, 1]
        // [0, 0, 0, 1]
        // [0, 1, 1, 1]
        // gj
        // [1, 1, 0, 0, 1]
        // [0, 1, 1, 1, 1]
        // [0, 0, 0, 1, 0]
        // [0, 0, 0, 0, 0]
        // b
        // [1, 1, 0, 1]
        // x
        // [1, 0, 1, 0]
        four_by_four_system(
            vec![vec![0, 1], vec![0, 1, 3], vec![3], vec![1, 2, 3]],
            bitvec![1, 1, 0, 1].to_bitvec(),
            bitvec![1, 0, 1, 0].to_bitvec(),
        )
    }

    #[test]
    fn test_nullspace_zerovec() {
        // gamma
        // ([[1, 0, 0, 0],
        // [0, 0, 0, 0],
        // [0, 0, 0, 0],
        // [0, 0, 0, 1]])
        // b
        //[0, 0, 0, 0]
        // x
        // [0, 1, 1, 0]
        // sG^deltaG
        // [0, 1, 1, 0]
        // is codeword
        // False
        nullspace_check(
            vec![vec![0], vec![], vec![], vec![3]],
            bitvec![0, 0, 0, 0].to_bitvec(),
            bitvec![0, 1, 0, 0].to_bitvec(),
            bitvec![0, 1, 1, 0].to_bitvec(),
            false,
        )
    }

    #[test]
    fn test_nullspace_nonzerovec() {
        // Gamma
        // array([[0, 1, 1, 0],
        //        [1, 0, 0, 1],
        //        [1, 0, 0, 1],
        //        [0, 1, 1, 0]])
        // sB^deltaB
        // array([1, 0, 0, 1])
        // sG^deltaG
        // array([0, 1, 1, 0])
        // xG
        // array([0, 1, 0, 0])
        // check_syndome
        // True
        nullspace_check(
            vec![vec![1, 2], vec![0, 3], vec![0, 3], vec![1, 2]],
            bitvec![1, 0, 0, 1].to_bitvec(),
            bitvec![0, 1, 0, 0].to_bitvec(),
            bitvec![0, 1, 1, 0].to_bitvec(),
            true,
        )
    }

    #[test]
    fn test_negative_case() {
        four_by_four_system_failure(
            vec![vec![0], vec![], vec![], vec![3]],
            bitvec![1, 1, 1, 0].to_bitvec(),
        )
    }
}
