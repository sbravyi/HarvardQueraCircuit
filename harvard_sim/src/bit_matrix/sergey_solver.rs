use bitvec::prelude::*;
use bitvec::{vec::BitVec, field::BitField};
use itertools::Itertools;

use super::column_matrix::ColumnMatrix;
use super::matrix::BitMatrix;

pub struct SergeySolver {
    augmented_system: BitMatrix,
    n: usize,
    x: ColumnMatrix,
    solution: BitVec,
}

fn debug_bitmatrix(m: &BitMatrix) {
    println!("{}", m.rows.iter().map(|r|{r.to_string()}).join("\n"));
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
        }
    }

    pub fn debug(&self) {
        println!("a");
        debug_bitmatrix(&self.augmented_system);
        println!("x");
        debug_col_matrix(&self.x);
        println!("sol");
        debug_bitvec(&self.solution);
    }

    pub fn copy_from_augmented_system(&mut self, m: &BitMatrix, b: &BitVec) {
        debug_assert_eq!(m.number_of_columns, self.n);
        debug_assert_eq!(b.len(), self.n);
        debug_assert_eq!(m.rows.len(), self.augmented_system.rows.len());
        let x_rank = if b.first_one().is_some() { self.n + 1 } else { self.n };
        for (idx, row) in m.rows.iter().enumerate() {
            self.augmented_system.rows[idx][..self.n].copy_from_bitslice(&row[..]);
        }
        for (idx, b) in b.iter().enumerate() {
            unsafe {
                self.augmented_system.rows[idx].set_unchecked(self.n, *b);
            }
        }
        self.x.reset();
        // quick identity matrix creation
        for idx in 0..x_rank {
            self.x.cols[idx].store(1 << idx);
        }
    }

    fn reformulate_x_from_augmented_syndrome(&mut self, current_equation: usize) -> Option<()> {
        let b_row = &self.augmented_system.rows[current_equation];
        let mut syndrome = bitvec![usize, Lsb0; 0; self.x.cols.len()];
        for (col_idx, col) in self.x.cols.iter().enumerate() {
            let mut syndrome_val = false;
            // syndrome val is the inner product of each x column with the current b row
            for b_one_idx in b_row.iter_ones() {
                syndrome_val ^= b_row[b_one_idx] & col[b_one_idx];
            }
            unsafe {
                syndrome.set_unchecked(col_idx, syndrome_val);
            }
        }
        let mut bad_cols = syndrome.iter_ones();
        let first_bad_col_idx = bad_cols.next();
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

    pub fn take_arbitrary_solution_from_x(&mut self) -> Option<()> {
        let last_row = self.x.number_of_rows - 1;
        if let Some(solution) = self.x.cols.iter().find_or_first(|col| col[last_row]) {
            self.solution[..self.n].copy_from_bitslice(&solution[..self.n]);
            Some(())
        } else {
            None
        }
    }

    pub fn solve(&mut self, u: &BitMatrix, b: &BitVec) -> Option<()> {
        self.copy_from_augmented_system(u, b);
        for i in 0..self.n {
            println!("equation {i}");
            self.debug();
            self.reformulate_x_from_augmented_syndrome(i)?;
        }
        self.take_arbitrary_solution_from_x()?;
        println!("closing remarks");
        self.debug();
        Some(())
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use crate::bit_matrix::sergey_solver::SergeySolver;


    fn four_by_four_system(
        mat_rows: Vec<Vec<usize>>,
        b: BitVec,
        x: BitVec,
    ) {
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
            bitvec![0, 1, 1, 0].to_bitvec()
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
            bitvec![1, 0, 1, 0].to_bitvec()
        )
    }
}



