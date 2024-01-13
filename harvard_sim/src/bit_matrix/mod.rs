use bitvec::prelude::*;

pub struct BitMatrix {
    rows: Vec<BitVec>,
}

impl BitMatrix {
    pub fn zeroes(rows: usize, cols: usize) -> Self {
        Self {
            rows: (0..rows).map(|_| bitvec![usize, Lsb0; 0; cols ]).collect(),
        }
    }
}
