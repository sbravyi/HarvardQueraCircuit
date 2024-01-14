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
        let n = self.solution.len();
        let mut i = n;
        while i > 0 {
            i -= 1;
            let mut tmp = b[i];
            let row = &u.rows[i];
            if row.first_one().is_none() {
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
