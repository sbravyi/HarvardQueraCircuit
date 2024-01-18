use bitvec::prelude::*;
use std::collections::{HashSet, BTreeSet};

#[derive(Default)]
pub struct SwapSymmetries {
    pub current_bit_pattern: u16,
    pub symmetries: HashSet<u16>,
}

impl SwapSymmetries {
    pub fn new() -> Self {
        Default::default()
    }

    fn inverted_nybbles(bitstring: u16) -> u16 {
        let mut inverted = 0;
        let to_invert = &BitSlice::<_, Lsb0>::from_element(&bitstring);
        for (bit_idx, bit) in to_invert.iter().by_vals().enumerate() {
            let idx = bit_idx as u16;
            let target_index = (3 - ((idx) % 4)) + (idx / 4) * 4;
            if bit {
                inverted |= 1 << target_index;
            }
        }
        inverted
    }

    fn bisection_swap(bitstring: u16) -> u16 {
        let mut swapped = 0;
        let to_swap = &BitSlice::<_, Lsb0>::from_element(&bitstring);
        for (byte_index, byte_chunk) in to_swap.chunks_exact(8).enumerate() {
            let byte_offset = (byte_index as u16) * 8;
            for (nybble_index, nybble_chunk) in byte_chunk.chunks_exact(4).enumerate() {
                let nybble_offset = byte_offset + (1 - nybble_index as u16) * 4;
                for (bit_pair_index, bit_pair) in nybble_chunk.chunks_exact(2).enumerate() {
                    let bit_pair_offset = nybble_offset + (bit_pair_index as u16) * 2;
                    let b0 = bit_pair[0];
                    let b1 = bit_pair[1];
                    let b0_idx = bit_pair_offset + 1;
                    let b1_idx = bit_pair_offset;
                    if b0 {
                        swapped |= 1 << b0_idx;
                    }
                    if b1 {
                        swapped |= 1 << b1_idx;
                    }
                }
            }
        }
        swapped
    }

    fn generate_symmetries_for_current_bit_pattern(&self) -> [u16; 4] {
        let original = self.current_bit_pattern;
        let inverted_nybbles = Self::inverted_nybbles(original);
        let bisection_swapped = Self::bisection_swap(original);
        let bisection_swapped_inverted_nybbles = Self::inverted_nybbles(bisection_swapped);
        [
            original,
            inverted_nybbles,
            bisection_swapped,
            bisection_swapped_inverted_nybbles,
        ]
    }

    pub fn increment_bit(&mut self, flip_bit: Option<u32>) {
        if let Some(flip_bit) = flip_bit {
            self.current_bit_pattern ^= 1 << flip_bit;
            return;
        };
        #[cfg(debug_assertions)]
        {
            self.current_bit_pattern = self.current_bit_pattern.saturating_add(1);
        }
        #[cfg(not(debug_assertions))]
        {
            self.current_bit_pattern += 1;
        }
    }

    #[inline(never)]
    pub fn check_for_symmetries(&mut self) -> Option<BTreeSet<u16>> {
        // generates the simulated bitstring by flipping the bit at the given bit index.
        // if the current bitstring is not a symmetry, swap symmetries generates all symmetries
        // of the current bitstring that will be simulated, and returns the number of symetries
        // so that the amplitude of the bitstring can be multiplied.
        //
        // note: this can probably be made faster by just generating symmetries on the fly rather than
        // tracking them in a hashset, and understanding if a symmetry has already been encountered based
        // on the smallest bit.
        if self.symmetries.contains(&self.current_bit_pattern) {
            return None;
        }
        let symmetries = self.generate_symmetries_for_current_bit_pattern();
        let deduplicate = BTreeSet::from_iter(symmetries);
        for symmetry in deduplicate.iter() {
            self.symmetries.insert(*symmetry);
        }
        Some(deduplicate)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    fn test_symmetries(symmetry_group: Vec<u16>) {
        let all_symmetries = BTreeSet::from_iter(symmetry_group.iter());
        for symmetry in symmetry_group.iter() {
            let mut ss = SwapSymmetries::new();
            let mut other_symmetries = all_symmetries.clone();
            other_symmetries.remove(&symmetry);
            for _ in 0..*symmetry {
                ss.increment_bit(None);
            }
            ss.check_for_symmetries().expect("to have symmetries");
            let mut ss_symmetries_without_self = ss.symmetries.clone();
            ss_symmetries_without_self.remove(symmetry);
            assert_eq!(ss_symmetries_without_self.len(), other_symmetries.len(), "no extra symmetries - {:?} = {:?}", ss.symmetries, other_symmetries);
            for symmetry in other_symmetries {
                assert!(
                    ss.symmetries.contains(symmetry),
                    "did not contain symmetry 0b{symmetry:b} | {:?}",
                    ss.symmetries.iter().map(|s| format!("0b{s:b}")).collect_vec()
                );
            }
        }
    }

    // p1:=0123456789ABCDEF; // let p1 be the original 16-digit Boolean pattern we are looking at. Let's say Sergey's check was applied (changing the order of Sergey's and this check may be something to play with) and this pattern needs to be used to simulate the corresponding Clifford. The digits 0..F here denote positions of Boolean values 0/1 in the given 16-tuple
    // p2:=32107654BA98FEDC; // a first swapped pattern (each slice of four bits is reversed)
    // if (p1>p2) {p1=NextInGrayCodeOrder(p1); break;} // break out of the entire loop and consider a next p1, this is a best-case scenario since no Clifford computation is required for this p1, and we also found it quickly. Here, I am assuming the canonical representative is the pattern least in the lexicographic order
    // p3:=54761032DCFE98BA; // a new swapped pattern (recursive bisection swap
    // if (p1>p3) {p1=NextInGrayCodeOrder(p1); break;} // again, we do not need to consider this p1 since it is represented by a different pattern---maybe p3, but maybe other
    // p4:=67452301EFCDAB89; // (subdivision swap of the above) last swapping pattern, there could be 4 more if pattern-dependent Z's can be handled---Sergey please clarify
    // if (p1>p4) {p1=NextInGrayCodeOrder(p1); break;} // p4 takes care of p1, so no need for p1
    // Cardinality:=NumberofDifferentPatternsInTheSet(p1,p2,p3,p4); // merge sort and count? As far as I can tell, the Cardinality should always be a power of 2. If we find ourselves on this line, p1 cannot be discarded as it is the smallest in its equivalence class and thus it is the canonical representative, and we must execute a Clifford circuit over it
    // Answer = Answer + Cardinality*SimulateClifford(p1); // call the Clifford simulator and remember to multiply the answer returned by the count of the number of patterns p1's equivalence class

    #[test]
    fn test_one() {
        test_symmetries(vec![0b1, 0b1000, 0b100000, 0b1000000]);
        test_symmetries( vec![0b1000, 0b1, 0b100000, 0b1000000]);
        test_symmetries( vec![0b0111, 0b1110, 0b10110000, 0b11010000]);
        test_symmetries( vec![0b10010110]);
    }
}
