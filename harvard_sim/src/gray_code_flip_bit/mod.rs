pub struct GrayCodeFlipBit {
    index: u64,
    number_of_bits: u32,
    finished: bool
}

impl GrayCodeFlipBit {
    pub fn new(number_of_bits: u32) -> Self {
        Self {
            index: 1,
            number_of_bits,
            finished: false
        }
    }
}


// Loosely inspired by:
// https://www.sciencedirect.com/science/article/pii/S0195669812001321#br000040
//
// return the bitstring index to flip when
// proceeding from the ith gray code bitstring to the
// i+1th gray code bitstring.
impl Iterator for GrayCodeFlipBit {
    type Item = u32;
    fn next(&mut self) -> Option<u32> {
        if self.finished {
            return None;
        }
        if self.index.ilog2() >= self.number_of_bits {
            // always flip the last bit before you leave
            self.finished = true;
            return Some(self.number_of_bits - 1);
        }
        let n = self.index ^ (self.index - 1);
        let set_bits = n.count_ones();
        self.index += 1;
        Some(set_bits - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flip_bits() {
        let flip_bits = GrayCodeFlipBit::new(4);
        let flip_bits_list: Vec<u32> = flip_bits.collect();
        assert_eq!(
            flip_bits_list,
            vec![0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 3]
        )
    }
}
