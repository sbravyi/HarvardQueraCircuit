use super::matrix::BitMatrix;
use bitvec::prelude::*;

pub fn is_in_nullspace(u: &BitMatrix, b: &BitVec) -> bool {
    // do a matrix multiplication ub = x. if you run into a nonzero value in x, then
    // b is not in the nullspace.
    for row in u.rows.iter() {
        let mut bit_overlap = 0;
        for (idx, c) in row.iter().enumerate() {
            let b_i = b[idx];
            if *c & b_i {
                bit_overlap += 1;
            }
        }
        if bit_overlap % 2 != 0 {
            return false;
        }
    }
    true
}
