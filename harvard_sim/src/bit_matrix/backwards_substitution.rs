use bitvec::prelude::*;
use itertools::Itertools;

use super::{gauss_jordan::GaussJordan, matrix::BitMatrix};

pub struct BackwardsSubstitution {
    pub solution: BitVec,
    pub gj: GaussJordan,
    pub intemediary_vec: BitVec,
}

impl BackwardsSubstitution {
    pub fn zero(n: usize) -> Self {
        Self {
            solution: bitvec![usize, Lsb0; 0; n],
            gj: GaussJordan::zero(n, n + 1),
            intemediary_vec: bitvec![usize, Lsb0; 0; n],
        }
    }

    pub fn solve_non_unique(&mut self, u: &BitMatrix, b: &BitVec) -> Option<()> {
        self.gj.copy_from_augmented_system(u, b);
        self.gj.go_to_echelon_form();
        #[cfg(debug_assertions)]
        println!(
            "gj({}):[\n{}\n]",
            self.gj.rank,
            self.gj
                .rows
                .iter()
                .map(|row| row.to_string())
                .collect_vec()
                .join("\n")
        );
        let vec = &mut self.intemediary_vec;
        let mut last_one = 0;
        for (idx, row) in self.gj.rows.iter().enumerate() {
            let val = row[self.gj.last_col_idx];
            unsafe { vec.set_unchecked(idx, val) }
            if val {
                last_one = idx;
            }
        }
        let rank = self.gj.rank;
        if last_one > rank {
            return None;
        }
        let solution = &mut self.solution;
        solution.set_elements(0x0);
        let mut colrank = rank;
        let mut col_idx: usize = u.number_of_columns;
        let last_row = &self.gj.rows[rank - 1];
        if let Some(last_pivot) = last_row.first_one() {
            let last_original_row = u.rows.len() - 1;
            colrank += last_original_row.saturating_sub(last_pivot);
        }
        #[cfg(debug_assertions)]
        println!("b: {vec}");
        for vec_idx in (0..colrank).rev() {
            if col_idx == 0 {
                return None;
            }
            col_idx -= 1;
            #[cfg(debug_assertions)]
            println!("col: {col_idx}, vec: {vec_idx}");
            if vec[vec_idx] {
                while !self.gj.rows[vec_idx][col_idx] {
                    if col_idx == 0 {
                        return None;
                    }
                    col_idx -= 1;
                    #[cfg(debug_assertions)]
                    println!("ranging over col: {col_idx}");
                }
                println!("ip pre check: {vec}");
                for (idx, mut bit) in vec.iter_mut().enumerate() {
                    let cached_v_idx = *bit;
                    *bit = cached_v_idx ^ self.gj.rows[idx][col_idx];
                }
                println!("ip post check: {vec}");
                println!("sol post check: {solution}");
                unsafe {
                    solution.set_unchecked(col_idx, true);
                }
                println!("sol post check: {solution}");
            } else if !b[vec_idx] && self.gj.rows[vec_idx][col_idx] && col_idx > vec_idx {
                col_idx = col_idx.saturating_sub(1);
            }
        }
        println!("ip on completion: {vec}");
        println!("sol on completion: {solution}");
        if vec.first_one().is_some() {
            None
        } else {
            Some(())
        }
        //Some(())
    }

    pub fn solve(&mut self, u: &BitMatrix, b: &BitVec) -> Option<()> {
        // TODO: make augmented system for gj and extract
        #[cfg(debug_assertions)]
        println!("debugging solver!");
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
            for j in (i + 1)..n {
                let rowval = row[j];
                let solval = self.solution[j];
                tmp ^= rowval & solval;
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
    use bitvec::prelude::*;

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
        solver
            .solve_non_unique(&matrix, &bits![1, 1, 0, 0].to_bitvec())
            .unwrap();
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
        solver.solve_non_unique(&matrix, &b).unwrap();
        assert_eq!(&solver.solution, &x);
    }

    #[test]
    fn test_non_unique() {
        let mut solver = BackwardsSubstitution::zero(4);
        // [1, 1, 0, 1]
        // [1, 1, 0, 1]
        // [0, 0, 1, 0]
        // [1, 1, 0, 0]
        let rows: Vec<Vec<usize>> = vec![vec![0, 1, 3], vec![0, 1, 3], vec![2], vec![0, 1]];
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let b = bits![1, 1, 0, 0].to_bitvec();
        solver.solve_non_unique(&matrix, &b).unwrap();
        assert_eq!(
            solver.gj.rows.iter().map(|bv| bv.to_string()).collect_vec(),
            vec![
                "[1, 1, 0, 1, 1]",
                "[0, 0, 1, 0, 0]",
                "[0, 0, 0, 1, 1]",
                "[0, 0, 0, 0, 0]"
            ]
        );
        let x = bits![0, 0, 0, 1].to_bitvec();
        assert_eq!(solver.solution.to_string(), x.to_string());
    }

    #[test]
    fn test_zeroes() {
        // [0, 1, 0, 0]
        // [1, 1, 0, 0]
        // [0, 0, 1, 1]
        // [0, 0, 1, 0]
        let rows: Vec<Vec<usize>> = vec![vec![1], vec![0, 1], vec![2, 3], vec![2]];
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let b = bits![0, 1, 0, 0].to_bitvec();
        let mut solver = BackwardsSubstitution::zero(4);
        solver.solve_non_unique(&matrix, &b).unwrap();
        let x = bits![1, 0, 0, 0].to_bitvec();
        assert_eq!(solver.solution.to_string(), x.to_string());
    }

    #[test]
    fn test_all_ones() {
        // gamma:  [1, 1, 1, 0]
        //         [1, 0, 0, 0]
        //         [1, 0, 1, 1]
        //         [0, 0, 1, 1]

        // b : [0, 1, 0, 1]
        let rows: Vec<Vec<usize>> = vec![vec![0, 1, 2], vec![0], vec![0, 2, 3], vec![2, 3]];
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let b = bits![0, 1, 0, 1].to_bitvec();
        let mut solver = BackwardsSubstitution::zero(4);
        solver.solve_non_unique(&matrix, &b).unwrap();
        assert_eq!(
            solver.gj.rows.iter().map(|bv| bv.to_string()).collect_vec(),
            vec![
                "[1, 1, 1, 0, 0]",
                "[0, 1, 1, 0, 1]",
                "[0, 0, 1, 1, 1]",
                "[0, 0, 0, 0, 0]"
            ]
        );
        solver.solve_non_unique(&matrix, &b).unwrap();
        let x = bits![1, 0, 1, 0].to_bitvec();
        assert_eq!(solver.solution.to_string(), x.to_string());
    }

    fn four_by_four_system(
        mat_rows: Vec<Vec<usize>>,
        b: BitVec,
        x: BitVec,
        expected_gj: Vec<&str>,
    ) {
        let mut matrix = BitMatrix::zeroes(4, 4);
        for (ridx, r) in mat_rows.iter().enumerate() {
            for cidx in r {
                matrix.set(ridx, *cidx, true)
            }
        }
        let mut solver = BackwardsSubstitution::zero(4);
        solver.solve_non_unique(&matrix, &b).unwrap();
        assert_eq!(
            solver.gj.rows.iter().map(|bv| bv.to_string()).collect_vec(),
            expected_gj
        );
        solver.solve_non_unique(&matrix, &b).unwrap();
        assert_eq!(solver.solution.to_string(), x.to_string());
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
            bitvec![0, 1, 0, 0].to_bitvec(),
            vec![
                "[1, 1, 0, 0, 1]",
                "[0, 1, 1, 1, 1]",
                "[0, 0, 0, 1, 0]",
                "[0, 0, 0, 0, 0]"
            ]
        )
    }
}
