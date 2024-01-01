// gf2.cpp

#include "../include/gf2.hpp"
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

gf2::Matrix gf2::transpose(const gf2::Matrix &matrix) {
    const size_t &cols = matrix.cols, &rows = matrix.rows;
    const std::vector<std::vector<gf2::Word>> &data = matrix.data;

    std::vector<std::vector<gf2::Word>> transposed;
    transposed.reserve(rows);

    for (uint32_t r = 0; r < rows; r++) {
        std::vector<gf2::Word> cur_row;
        cur_row.reserve(cols / gf2::BLOCK + (cols % gf2::BLOCK != 0));

        // Write c = BLOCK * cq + cr for cr < BLOCK
        // and similarly for r
        uint32_t cq, cr;
        const uint32_t rq = r / gf2::BLOCK;
        const gf2::Word rr_indicator = (gf2::Word(1) << (r % gf2::BLOCK));

        for (cq = 0; cq < cols / gf2::BLOCK; cq++) {
            gf2::Word t = 0;
            for (cr = 0; cr < gf2::BLOCK; cr++) {
                if (data[BLOCK * cq + cr][rq] & rr_indicator) {
                    t |= (gf2::Word(1) << cr);
                }
            }
            cur_row.push_back(t);
        }

        // Extra section from if cols % BLOCK != 0
        if (cols % gf2::BLOCK) {
            gf2::Word t = 0;
            for (cr = 0; cr < cols % gf2::BLOCK; cr++) {
                if (data[BLOCK * cq + cr][rq] & rr_indicator) {
                    t |= (gf2::Word(1) << cr);
                }
            }
            cur_row.push_back(t);
        }
        transposed.push_back(cur_row);
    }
    return Matrix{rows, cols, std::move(transposed)};
    //return Matrix{rows, cols, transposed};
}

// See https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
// for the algorithm (Gaussian Elimination), except transposed
// It's basically this https://en.wikipedia.org/wiki/Kernel_(linear_algebra)#Computation_by_Gaussian_elimination
//
// The algorithm transposes so that we can quickly xor rows
gf2::Matrix gf2::nullspace(const gf2::Matrix &matrix) {
    gf2::Matrix transposed = gf2::transpose(matrix);
    const size_t &cols = matrix.cols, &rows = matrix.rows;
    std::vector<std::vector<gf2::Word>> &transposed_data = transposed.data;

    // For each linearly indepedent new row, we mark the 
    // unique column where it owns the pivot, or -1 if no such column
    std::vector<int32_t> marked_col(rows, -1);
    // Whether a column is marked
    std::vector<bool> is_marked(cols, false);
    
    uint32_t rv_len = cols / gf2::BLOCK + (cols % gf2::BLOCK!= 0);
    for (uint32_t r = 0; r < rows; r++) {
        int32_t pivot = -1;
        for (uint32_t c = 0; c < cols; c++) {
            if (transposed_data[r][c / gf2::BLOCK] & (gf2::Word(1) << (c % gf2::BLOCK))) {
                pivot = c;
                break;
            }
        }
        if (pivot == -1) continue;

        marked_col[r] = pivot;
        is_marked[pivot] = true;

        uint32_t pivot_q = pivot / gf2::BLOCK;
        gf2::Word pivot_indicator = (gf2::Word(1) << (pivot % gf2::BLOCK));

        // Look over all other rows r2, and if there is a 1 
        // in the pivot column, then xor rom r into row r2
        for (uint32_t r2 = 0; r2 < rows; r2++) {
            if (r2 == r) continue;
            if (transposed_data[r2][pivot_q] & pivot_indicator) {
                for (uint32_t cq = 0; cq < rv_len; cq++) {
                    transposed_data[r2][cq] ^= transposed_data[r][cq];
                }
            }
        }
    }
    std::vector<std::vector<gf2::Word>> kernel;

    // Each unmarked column can be written as a linear combination
    // of the marked columns
    for (uint32_t c = 0; c < cols; c++) {
        if (!is_marked[c]) {
            std::vector<gf2::Word> vec(cols / gf2::BLOCK + (cols % gf2::BLOCK != 0), 0);

            uint32_t cq = c / gf2::BLOCK;
            gf2::Word cr_indicator = (gf2::Word(1) << (c % gf2::BLOCK));
            vec[cq] |= cr_indicator;

            for (uint32_t r = 0; r < rows; r++) {
                if (transposed_data[r][cq] & cr_indicator) {
                    vec[marked_col[r] / gf2::BLOCK] |= (gf2::Word(1) << (marked_col[r] % gf2::BLOCK));
                }
            }
            kernel.push_back(std::move(vec));
        }
    }
    return Matrix{kernel.size(), cols, std::move(kernel)};
}
