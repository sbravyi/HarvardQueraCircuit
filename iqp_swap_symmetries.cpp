
#include "iqp_swap_symmetries.h"
#include <set>
#include <vector>
using namespace std;

IQPSwapSymmetries::IQPSwapSymmetries()
    : tracked_symmetries{ {} }
{
}

// k = 2 -> invert order of bits in the x register.

unsigned long inverted_nybbles(unsigned long bitstring) {
    unsigned long inverted = 0;
    for (unsigned long i = 0; i < 64; i++) {
        unsigned long bit_index = 1UL << i;
        if (bitstring & bit_index) {
            unsigned long target_index = (3 - ((i) % 4)) + (i / 4) * 4;
            inverted |= 1 << target_index;
        }
    }
    return inverted;
}

unsigned long bisection_swap(unsigned long bitstring) {
    unsigned long swapped = 0;
    const unsigned long bytes_in_long = 64/8;
    const unsigned long nybbles_in_byte = 8/4;
    const unsigned long pairs_in_nybble = 4/2;
    for (unsigned long byte_idx = 0; byte_idx < bytes_in_long; byte_idx++) {
        unsigned long byte_offset = byte_idx * 8;
        for (unsigned long nybble_index = 0; nybble_index < nybbles_in_byte; nybble_index++) {
            unsigned long target_nybble_offset = byte_offset + (1UL - nybble_index) * 4;
            unsigned long original_nybble_offset = byte_offset + nybble_index * 4;
            for (unsigned long bit_pair_index = 0; bit_pair_index < pairs_in_nybble; bit_pair_index++) {
                unsigned long target_bit_pair_offset = target_nybble_offset + (bit_pair_index) * 2;
                unsigned long original_bit_pair_offset = original_nybble_offset + (bit_pair_index) * 2;
                unsigned long b0 = bitstring & (0b1 << original_bit_pair_offset);
                unsigned long b1 = bitstring & (0b10 << original_bit_pair_offset);
                if (b0) {
                    swapped |= 0b10 << target_bit_pair_offset;
                }
                if (b1) {
                    swapped |= 0b1 << target_bit_pair_offset;
                }
            }
        }
    }
    return swapped;
}

vector<unsigned long> generate_symmetries(unsigned long bitstring) {
    std::vector<unsigned long> symmetries;
    symmetries.reserve(4);
    symmetries.push_back(bitstring); // identity
    unsigned long inverted = inverted_nybbles(bitstring);
    if (inverted != bitstring) {
        symmetries.push_back(inverted);
    }
    unsigned long swapped = bisection_swap(bitstring);
    if ((swapped != bitstring) && (swapped != inverted)) {
        symmetries.push_back(swapped);
    }
    if (symmetries.size() == 3) { // (inverted_nybbles * bisection_swap) is only a distinct symmetry if neither permutation is identity
        unsigned long bisection_invert = inverted_nybbles(swapped);
        symmetries.push_back(bisection_invert);
    }
    return symmetries;
}

unsigned IQPSwapSymmetries::is_symmetry_or_should_multiply_amplitude(unsigned long bitstring) {
    auto it = tracked_symmetries.find(bitstring);
    if (it != tracked_symmetries.end()) {
        return 0;
    }
    auto symmetries = generate_symmetries(bitstring);
    for (unsigned long symmetry: symmetries) {
        tracked_symmetries.insert(symmetry);
    }
    return symmetries.size();
}
