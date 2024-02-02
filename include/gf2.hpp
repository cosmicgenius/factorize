// gf2.cpp
#ifndef GF2_HPP_
#define GF2_HPP_

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// All matrices are F2 matrices encoded by bits 
// in integer words

namespace gf2 {
    typedef uint64_t Word;
    typedef uint32_t HalfWord;

    // Number of bits per Word
    constexpr uint32_t BLOCK = 8 * sizeof(Word),
                       HALF_BLOCK = 8 * sizeof(HalfWord);

    // GF(2) matrix stored efficiently with each column
    // as a vector of Words. 
    //
    // Assumes that data.size() >= cols
    // (only the first ``cols'' columns are used) and that
    // for each column, column.size() * BLOCK >= rows
    // (only the first ``rows'' bits are used)
    struct MatrixCM {
        size_t cols;
        size_t rows;
        std::vector<std::vector<Word>> col_data;
    };

    // For each row, store the indices corresponding to 
    // the columns with a 1, where the indices are guaranteed to 
    // be sorted in increasing order
    struct SparseMatrix {
        size_t rows;
        size_t cols;
        std::vector<std::vector<uint32_t>> row_data;
    };

    // A polynomial in M_{m, n}(F2)[X], or alternatively a matrix of 
    // polynomials in F2[X]
    // Stored as a vector of matrices, with element i representing the 
    // coefficient of X^i. The matrices are row-major static arrays
    // with each row being some int type (so the matrices can only be so wide).
    template<typename TRow, size_t TNumRows>
    struct BlockPoly {
        std::vector<std::array<TRow, TNumRows>> data;

        // Multiply in place by square matrix/linear poly
        BlockPoly<TRow, TNumRows>& mul_eq(const std::array<TRow, 8 * sizeof(TRow)> &con,
                                          const std::array<TRow, 8 * sizeof(TRow)> &lin);
        // Multiply by non square matrix/linear poly, must create new matrix series
        template<typename TOtherRow> 
        BlockPoly<TOtherRow, TNumRows> mul(const std::array<TOtherRow, 8 * sizeof(TRow)> &con,
                                           const std::array<TOtherRow, 8 * sizeof(TRow)> &lin)
                                       const;
    };

    template<typename TRow> 
    bool is_singular(std::array<TRow, 8 * sizeof(TRow)> matrix);

    MatrixCM transpose(const MatrixCM &matrix);
    SparseMatrix transpose(const SparseMatrix &sparse_matrix);

    // Returns, in matrix form, a list of vectors 
    // spanning the nullspace of the matrix.
    MatrixCM nullspace(const MatrixCM &matrix);

    // Returns a block vector of gf2::HALF_BLOCK vectors that are 
    // in the nullspace of the matrix sparse_matrix. With high probability,
    // the vectors (columns of the returned block vector) are not 0.
    std::vector<gf2::HalfWord> sparse_nullspace(SparseMatrix sparse_matrix);
};

#endif // GF2_HPP_
