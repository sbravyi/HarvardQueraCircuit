#ifndef SWAP_SYMMETRIES_H
#define SWAP_SYMMETRIES_H
#include <set>

class IQPSwapSymmetries {
private:
    std::set<unsigned> tracked_symmetries;
public:
    IQPSwapSymmetries();
    unsigned is_symmetry_or_should_multiply_amplitude(unsigned long bitstring);
};

#endif