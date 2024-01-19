#include <iostream>
#include <bitset>
#include <vector>
#include <set>
using namespace std;



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

std::string print_symms(const vector<unsigned long>& symms) {
    std::string result;
    std::for_each(symms.begin(), symms.end(),
            [&result](unsigned long s) {
                std::bitset<sizeof(short unsigned) * 8> binaryRepresentation(s);
                std::string binaryString = binaryRepresentation.to_string();
                result += (",0b" + binaryString +  "");
            });
    return result;
}


void run_test(unsigned long bitstring, unsigned expected_symmetries) {
    auto symms = generate_symmetries(bitstring);
    unsigned n_symmetries = symms.size();
    std::bitset<sizeof(unsigned long) * 8> binaryRepresentation(bitstring);
    std::string binaryString = binaryRepresentation.to_string();
    if (n_symmetries == expected_symmetries) {
        std::cout << "SUCCESS!: 0b" + binaryString + " got " + std::to_string(expected_symmetries) + " symmetries!"<< std::endl;
        std::cout << "Symmetries: " << print_symms(symms) << std::endl;
    } else {
        std::cout << "FAILURE!: 0b" + binaryString + " expected " + std::to_string(expected_symmetries) + " symmetries but got " + std::to_string(n_symmetries) + " instead" << std::endl;
    }
}

bool areVectorsEqual(const std::vector<unsigned long>& vector1, const std::vector<unsigned long>& vector2) {
    // Check if the sizes of the vectors are equal
    if (vector1.size() != vector2.size()) {
        return false;
    }

    // Sort the vectors to ensure elements are in the same order
    std::vector<unsigned long> sortedVector1 = vector1;
    std::vector<unsigned long> sortedVector2 = vector2;
    std::sort(sortedVector1.begin(), sortedVector1.end());
    std::sort(sortedVector2.begin(), sortedVector2.end());

    // Compare the sorted vectors element by element
    return sortedVector1 == sortedVector2;
}

void test_symmetry_group(const vector<unsigned long>& group) {
    for (unsigned long s: group) {
        auto symms = generate_symmetries(s);
        if (!areVectorsEqual(symms, group)) {
            std::cout << "FAILURE! [group test]: " << print_symms(symms)  << " != " << print_symms(group) << std::endl;
            return;
        }
    }
    std::cout << "SUCCESS! [group test]: " << print_symms(group) << endl;
}

void test_u16_nonoverlaping_cover() {
    uint16_t maxUnsigned16Bit = std::numeric_limits<uint16_t>::max();
    std::set<std::vector<unsigned long>> all_symms;
    for (unsigned x = 0; x < maxUnsigned16Bit; x++) {
        auto symms = generate_symmetries(x);
        all_symms.insert(symms);
        unsigned n_seen = 0;
        for (std::vector<unsigned long> x: all_symms) {
            
        }
    }
}

int main() {
    test_symmetry_group({0b1, 0b1000, 0b100000, 0b1000000});
    test_symmetry_group({0b10010110});
    test_u16_nonoverlaping_cover();
    return 0;
}