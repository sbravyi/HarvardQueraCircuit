use bitvec::prelude::*;

pub struct BitMatrix {
    pub number_of_columns: usize,
    pub rows: Vec<BitVec>,
}

impl BitMatrix {
    pub fn zeroes(rows: usize, cols: usize) -> Self {
        Self {
            number_of_columns: cols,
            rows: (0..rows).map(|_| bitvec![usize, Lsb0; 0; cols ]).collect(),
        }
    }

    pub fn get(&self, r: usize, c: usize) -> bool {
        self.rows[r][c]
    }

    pub fn set(&mut self, r: usize, c: usize, v: bool) {
        unsafe { self.rows[r].set_unchecked(c, v) }
    }

    pub fn flip(&mut self, r: usize, c: usize) {
        unsafe {
            let mut val = self.rows[r].get_unchecked_mut(c);
            *val ^= true;
        }
    }
}
