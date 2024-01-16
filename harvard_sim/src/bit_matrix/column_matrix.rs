use bitvec::prelude::*;

// alloc-once matrix that allows efficient column operations
pub struct ColumnMatrix {
    pub cols: Vec<BitVec>,
    pub removed_cols: Vec<BitVec>,
    pub number_of_rows: usize
}

impl ColumnMatrix {
    pub fn zeroes(rows: usize, cols: usize) -> Self {
        Self {
            cols: (0..cols).map(|_| bitvec![usize, Lsb0; 0; rows ]).collect(),
            removed_cols: Vec::with_capacity(cols),
            number_of_rows: rows
        }
    }

    pub fn remove_col(&mut self, col_idx: usize) {
        let removed = self.cols.remove(col_idx);
        self.removed_cols.push(removed);
    }

    pub fn reset(&mut self) {
        self.cols.append(&mut self.removed_cols);
    }

    pub fn pop_from_removed(&mut self) -> BitVec {
        self.removed_cols.pop().unwrap()
    }

    pub fn push_back_into_usage(&mut self, to_add: BitVec) {
        self.cols.push(to_add);
    }

    pub fn put_back_in_removed(&mut self, to_add: BitVec) {
        self.removed_cols.push(to_add);
    }
}
