// gf2.cpp
#ifndef GF2_HPP_
#define GF2_HPP_

#include <cstddef>
#include <cstdint>
#include <vector>

// NOTE: all matrices use column major order because 
// that's what linear algebra uses.

namespace gf2 {
    typedef uint64_t Word;

    // Number of bits per Word
    const uint32_t BLOCK = 8 * sizeof(Word);

    // GF(2) matrix stored efficiently with each column
    // as a vector of Words. 
    //
    // Assumes that data.size() >= cols
    // (only the first ``cols'' columns are used) and that
    // for each column, column.size() * BLOCK >= rows
    // (only the first ``rows'' bits are used)
    struct Matrix {
        size_t cols;
        size_t rows;
        std::vector<std::vector<Word>> data;
    };

    Matrix transpose(const Matrix &matrix);

    // Returns, in matrix form, a list of vectors 
    // spanning the nullspace of the matrix.
    Matrix nullspace(const Matrix &matrix);
};

#endif // GF2_HPP_
