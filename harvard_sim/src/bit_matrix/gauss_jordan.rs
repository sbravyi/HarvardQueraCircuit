use bitvec::prelude::*;

use super::matrix::BitMatrix;

pub struct GaussJordan {
    pub number_of_columns: usize,
    pub last_col_idx: usize,
    active_column: usize,
    pub rows: Vec<BitVec>,
    pub rank: usize,
}

impl GaussJordan {
    pub fn zero(rows: usize, cols: usize) -> Self {
        Self {
            number_of_columns: cols,
            last_col_idx: cols - 1,
            active_column: 0,
            rows: (0..rows).map(|_| bitvec![usize, Lsb0; 0; cols]).collect(),
            rank: 0,
        }
    }

    // allows one to run gauss jordan on matrices of the same shape many times
    // without extra allocations, provided the matrix is the same shape.
    pub fn copy_from_matrix(&mut self, m: &BitMatrix) {
        debug_assert_eq!(m.number_of_columns, self.number_of_columns);
        debug_assert_eq!(m.rows.len(), self.rows.len());
        self.active_column = 0;
        self.rank = 0;
        for (idx, row) in m.rows.iter().enumerate() {
            self.rows[idx].copy_from_bitslice(&row[..]);
        }
    }

    // allows one to run gauss jordan on matrices of the same shape many times
    // without extra allocations, provided the matrix is the same shape.
    pub fn copy_from_augmented_system(&mut self, m: &BitMatrix, b: &BitVec) {
        debug_assert_eq!(m.number_of_columns, self.last_col_idx);
        debug_assert_eq!(m.rows.len(), self.rows.len());
        self.active_column = 0;
        self.rank = 0;
        for (idx, row) in m.rows.iter().enumerate() {
            self.rows[idx][..self.last_col_idx].copy_from_bitslice(&row[..]);
        }
        for (idx, b) in b.iter().enumerate() {
            unsafe {
                self.rows[idx].set_unchecked(self.last_col_idx, *b);
            }
        }
    }

    pub fn go_to_echelon_form(&mut self) {
        while self.is_not_in_echelon_form() {
            self.pivot_active_column();
            self.go_to_next_column();
        }
        self.rows.sort_by_cached_key(|row| {
            if let Some(pivot) = row.first_one() {
                // kind of a hack to make rank equal
                // to U in Ux = b, since we expect
                // this gauss-jordan to only
                // be used on augmented systems
                if pivot != self.last_col_idx {
                    self.rank += 1;
                }
                pivot
            } else {
                usize::MAX
            }
        });
    }

    fn is_not_in_echelon_form(&self) -> bool {
        self.active_column < self.number_of_columns
    }

    fn pivot_active_column(&mut self) {
        if let Some(pivot) = self.find_and_remove_pivot() {
            self.pivot_rows_that_start_in_active_column_with(&pivot);
            self.rows.push(pivot);
        }
    }

    fn find_and_remove_pivot(&mut self) -> Option<BitVec> {
        let mut row_index = 0;
        while row_index < self.rows.len() {
            if self.row_at_index_start_at_active_column(row_index) {
                let row = self.get_and_remove_row_at_index(row_index);
                return Some(row);
            }
            row_index += 1;
        }
        None
    }

    fn get_and_remove_row_at_index(&mut self, index: usize) -> BitVec {
        self.rows.swap_remove(index)
    }

    fn pivot_rows_that_start_in_active_column_with(&mut self, pivot: &BitVec) {
        let mut row_index = 0;
        while row_index < self.rows.len() {
            if self.row_at_index_start_at_active_column(row_index) {
                for (idx, mut bit) in self.rows[row_index].iter_mut().enumerate() {
                    unsafe {   
                        *bit ^= *pivot.get_unchecked(idx);
                    }
                }
            }
            row_index += 1;
        }
    }

    fn row_at_index_start_at_active_column(&self, index: usize) -> bool {
        self.rows[index]
            .first_one()
            .map(|column| column == self.active_column)
            .unwrap_or(false)
    }

    fn go_to_next_column(&mut self) {
        self.active_column += 1;
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn do_nothing_if_already_in_echelon_form() {
        let rows: Vec<Vec<usize>> = vec![vec![0, 1, 2], vec![1, 2, 3], vec![3, 4, 5], vec![5, 6]];
        let mut matrix = BitMatrix::zeroes(4, 7);
        for (ridx, r) in rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let mut gj = GaussJordan::zero(4, 7);
        gj.copy_from_matrix(&matrix);
        gj.go_to_echelon_form();
        assert_eq!(
            gj.rows.iter().map(|bv| bv.to_string()).collect_vec(),
            matrix.rows.iter().map(|bv| bv.to_string()).collect_vec()
        );
    }

    #[test]
    fn compute_the_good_echelon_form() {
        let rows: Vec<Vec<usize>> = vec![
            vec![0, 1, 2],
            vec![1, 2, 3],
            vec![0, 3],
            vec![3, 4, 5],
            vec![0, 4, 6],
            vec![5, 6],
        ];
        let mut matrix = BitMatrix::zeroes(rows.len(), 7);
        for (ridx, r) in rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let mut gj = GaussJordan::zero(rows.len(), 7);
        gj.copy_from_matrix(&matrix);
        gj.go_to_echelon_form();

        assert_eq!(
            gj.rows.iter().map(|bv| bv.to_string()).collect_vec(),
            vec![
                "[1, 1, 1, 0, 0, 0, 0]",
                "[0, 1, 1, 1, 0, 0, 0]",
                "[0, 0, 0, 1, 1, 1, 0]",
                "[0, 0, 0, 0, 0, 1, 1]",
                "[0, 0, 0, 0, 0, 0, 0]",
                "[0, 0, 0, 0, 0, 0, 0]"
            ]
        );
    }

    #[test]
    fn test_b2_input() {
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
        assert_eq!(
            gj.rows.iter().map(|bv| bv.to_string()).collect_vec(),
            vec![
                "[1, 1, 0, 0]",
                "[0, 1, 0, 1]",
                "[0, 0, 0, 0]",
                "[0, 0, 0, 0]",
            ]
        )
    }
}
