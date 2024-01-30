// gf2.cpp

#include "../include/gf2.hpp"
#include <algorithm>
#include <bitset>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

template struct gf2::BlockPoly<gf2::HalfWord, gf2::HALF_BLOCK>;
template struct gf2::BlockPoly<gf2::Word, gf2::HALF_BLOCK>;

gf2::MatrixCM gf2::transpose(const gf2::MatrixCM &matrix) {
    const size_t &cols = matrix.cols, &rows = matrix.rows;
    const std::vector<std::vector<gf2::Word>> &col_data = matrix.col_data;

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
                if (col_data[BLOCK * cq + cr][rq] & rr_indicator) {
                    t |= (gf2::Word(1) << cr);
                }
            }
            cur_row.push_back(t);
        }

        // Extra section from if cols % BLOCK != 0
        if (cols % gf2::BLOCK) {
            gf2::Word t = 0;
            for (cr = 0; cr < cols % gf2::BLOCK; cr++) {
                if (col_data[BLOCK * cq + cr][rq] & rr_indicator) {
                    t |= (gf2::Word(1) << cr);
                }
            }
            cur_row.push_back(t);
        }
        transposed.push_back(std::move(cur_row));
    }
    return gf2::MatrixCM{rows, cols, std::move(transposed)};
    //return Matrix{rows, cols, transposed};
}

gf2::SparseMatrix gf2::transpose(const gf2::SparseMatrix &sparse_matrix) {
    const size_t &rows = sparse_matrix.rows, &cols = sparse_matrix.cols;
    const std::vector<std::vector<uint32_t>> &data = sparse_matrix.row_data;
    
    std::vector<std::vector<uint32_t>> transposed(cols);
    for (uint32_t r = 0; r < rows; r++) {
        for (const uint32_t c : data[r]) {
            transposed[c].push_back(r);
        }
    }
    return gf2::SparseMatrix{rows, cols, std::move(transposed)};
}

// See https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
// for the algorithm (Gaussian Elimination), except transposed
// It's basically this https://en.wikipedia.org/wiki/Kernel_(linear_algebra)#Computation_by_Gaussian_elimination
//
// The algorithm transposes so that we can quickly xor rows
gf2::MatrixCM gf2::nullspace(const gf2::MatrixCM &matrix) {
    gf2::MatrixCM transposed = gf2::transpose(matrix);
    const size_t &cols = matrix.cols, &rows = matrix.rows;
    std::vector<std::vector<gf2::Word>> &transposed_data = transposed.col_data;

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
    return gf2::MatrixCM{kernel.size(), cols, std::move(kernel)};
}

// Multiply in place by square matrix/linear poly
template<typename TRow, size_t TNumRows>
gf2::BlockPoly<TRow, TNumRows>& gf2::BlockPoly<TRow, TNumRows>::mul_eq(
        const std::array<TRow, 8 * sizeof(TRow)> &con,
        const std::array<TRow, 8 * sizeof(TRow)> &lin) {
    uint32_t N = this->data.size();
    this->data.push_back(std::array<TRow, TNumRows>{});

    // This is the bottleneck multiplication step, so we optimize by caching 
    // (ok it's not a cache, it's just precomputation, but don't worry about that)
    // the xors for every possible combination of the con/lin rows for a byte of 
    // rows
    //
    // This is (nominally) an 8x speed improvement, less for small N since constant precomputation,
    // but N > 10,000 generally, so this is fine 
    std::array<std::array<TRow, 256>, sizeof(TRow)> con_byte_cache, lin_byte_cache;
    for (uint32_t r_byte = 0; r_byte < sizeof(TRow); r_byte++) {
        std::fill(con_byte_cache[r_byte].begin(), con_byte_cache[r_byte].end(), 0);
        std::fill(lin_byte_cache[r_byte].begin(), lin_byte_cache[r_byte].end(), 0);
        for (uint32_t byte_val = 0; byte_val < (1 << 8); byte_val++) {
            TRow indicator = 1;
            for (uint32_t byte_idx = 0; byte_idx < 8; byte_idx++, indicator <<= 1) {
                if (byte_val & indicator) {
                    con_byte_cache[r_byte][byte_val] ^= con[8 * r_byte + byte_idx];
                    lin_byte_cache[r_byte][byte_val] ^= lin[8 * r_byte + byte_idx];
                }
            }
        }
    }
    constexpr TRow byte_mask = (1 << 8) - 1;

    for (int deg = N - 1; deg >= 0; deg--) {
        // data[deg + 1] += data[deg] * lin 
        for (uint32_t r = 0; r < TNumRows; r++) {
            /*TRow c_indicator = 1;
            for (uint32_t c = 0; c < 8 * sizeof(TRow); c++, c_indicator <<= 1) {
                if (this->data[deg][r] & c_indicator) this->data[deg + 1][r] ^= lin[c];
            }*/
            TRow val = this->data[deg][r];
            for (uint32_t c_byte = 0; c_byte < sizeof(TRow); c_byte++, val >>= 8) {
                this->data[deg + 1][r] ^= lin_byte_cache[c_byte][val & byte_mask];
            }
        }
        // data[deg] *= con
        for (uint32_t r = 0; r < TNumRows; r++) {
            /*TRow old_row = this->data[deg][r];
            this->data[deg][r] = 0;

            TRow c_indicator = 1;
            for (uint32_t c = 0; c < 8 * sizeof(TRow); c++, c_indicator <<= 1) {
                if (old_row & c_indicator) this->data[deg][r] ^= con[c];
            }*/
            TRow val = this->data[deg][r];
            this->data[deg][r] = 0;
            for (uint32_t c_byte = 0; c_byte < sizeof(TRow); c_byte++, val >>= 8) {
                this->data[deg][r] ^= con_byte_cache[c_byte][val & byte_mask];
            }
        }
    }

    return *this;
}

// Multiply by non square matrix/linear poly, must create new matrix series
template<typename TRow, size_t TNumRows>
template<typename TOtherRow> 
gf2::BlockPoly<TOtherRow, TNumRows> gf2::BlockPoly<TRow, TNumRows>::mul(
        const std::array<TOtherRow, 8 * sizeof(TRow)> &con,
        const std::array<TOtherRow, 8 * sizeof(TRow)> &lin) const {
    uint32_t N = this->data.size();
    gf2::BlockPoly<TOtherRow, TNumRows> res;
    res.data.reserve(N + 1);

    for (uint32_t i = 0; i <= N; i++) {
        res.data.push_back(std::array<TOtherRow, TNumRows>{});
    }

    for (int deg = N - 1; deg >= 0; deg--) {
        // res.data[deg + 1] += data[deg] * lin 
        for (uint32_t r = 0; r < TNumRows; r++) {
            TRow c_indicator = 1;
            for (uint32_t c = 0; c < 8 * sizeof(TRow); c++, c_indicator <<= 1) {
                if (this->data[deg][r] & c_indicator) res.data[deg + 1][r] ^= lin[c];
            }
        }
        // res.data[deg] = data[deg] * con
        for (uint32_t r = 0; r < TNumRows; r++) {
            TRow c_indicator = 1;
            for (uint32_t c = 0; c < 8 * sizeof(TRow); c++, c_indicator <<= 1) {
                if (this->data[deg][r] & c_indicator) res.data[deg][r] ^= con[c];
            }
        }
    }
    return res;
}

// Four different mulitiplication implementations for multiplying the following 
// three data types: 
//  A. sparse_matrix, stored as vector of vector of indices
//  B. block vector, stored as vector of ints
//  C. single block matrix, stored as array of ints
// We support:
//  1. A   * B,
//  2. B^T * B,
//  3. B   * C
//  4. C   * C,

// Does res += sparse_matrix * v, where v and res are stored in row-major order,
// Assumes v and res are N x gf2::HALF_BLOCK, and sparse_matrix is N x N,
// for some fixed N.
void half_word_sparse_mul_vec(
        std::vector<gf2::HalfWord> &res, 
        const std::vector<gf2::HalfWord> &v, 
        const gf2::SparseMatrix &sparse_matrix) {
    const uint32_t N = sparse_matrix.cols;
    const std::vector<std::vector<uint32_t>> &data = sparse_matrix.row_data;

    // Multiplies half words at a time
    for (uint32_t r = 0; r < N; r++) {
        for (const uint32_t i : data[r]) {
            res[r] ^= v[i];
        }
    }
}

// Returns lhs^T rhs for lhs and rhs N x gf2::HALF_BLOCK, 
// everything in row-major order. Assumes that lhs.size() == rhs.size()
std::array<gf2::HalfWord, gf2::HALF_BLOCK> half_word_vec_dot_vec(
        const std::vector<gf2::HalfWord> &lhs, const std::vector<gf2::HalfWord> &rhs) {
    const uint32_t N = lhs.size();

    std::array<gf2::HalfWord, gf2::HALF_BLOCK> res{};
    // Multiplies half words at a time
    for (uint32_t r = 0; r < gf2::HALF_BLOCK; r++) {
        gf2::HalfWord indicator = gf2::HalfWord(1) << r;
        for (uint32_t i = 0; i < N; i++) {
            if (lhs[i] & indicator) {
                res[r] ^= rhs[i];
            }
        }
    }
    return res;
}

// Multiply N x gf2::HALF_BLOCK by gf2::HALF_BLOCK x gf2:;HALF_BLOCK,
// and adds it to res, i.e. does res += lhs * rhs.
void half_word_vec_mul_sq(
        std::vector<gf2::HalfWord> &res,
        const std::vector<gf2::HalfWord> &lhs,
        const std::array<gf2::HalfWord, gf2::HALF_BLOCK> &rhs) {
    
    const size_t N = lhs.size();
    for (uint32_t r = 0; r < N; r++) {
        for (uint32_t i = 0; i < gf2::HALF_BLOCK; i++) {
            if (lhs[r] & (gf2::HalfWord(1) << i)) {
                res[r] ^= rhs[i];
            }
        }
    }
}

// Multiply two square gf2::HALF_BLOCK x gf2::HALF_BLOCK matrices
std::array<gf2::HalfWord, gf2::HALF_BLOCK> half_word_sq_mul_sq(
    const std::array<gf2::HalfWord, gf2::HALF_BLOCK> &lhs, 
    const std::array<gf2::HalfWord, gf2::HALF_BLOCK> &rhs) {
    
    std::array<gf2::HalfWord, gf2::HALF_BLOCK> res{};
    for (uint32_t r = 0; r < gf2::HALF_BLOCK; r++) {
        for (uint32_t i = 0; i < gf2::HALF_BLOCK; i++) {
            if (lhs[r] & (gf2::HalfWord(1) << i)) {
                res[r] ^= rhs[i];
            }
        }
    }
    return res;
}

// Gaussian elimination
template<typename TRow> 
bool gf2::is_singular(std::array<TRow, 8 * sizeof(TRow)> matrix) {
    for (uint32_t r = 0; r < 8 * sizeof(TRow); r++) {
        TRow indicator = 1;
        int32_t pivot = -1;
        for (size_t c = 0; c < 8 * sizeof(TRow); c++, indicator <<= 1) {
            if (matrix[r] & indicator) {
                pivot = c;
                break;
            }
        }
        if (pivot == -1) return true;
        for (uint32_t r2 = 0; r2 < 8 * sizeof(TRow); r2++) {
            if (r2 == r) continue;
            matrix[r2] ^= matrix[r];
        }
    }
    return false;
}

template bool gf2::is_singular<gf2::HalfWord>(std::array<gf2::HalfWord, gf2::HALF_BLOCK> matrix);

// Returns P given E[t], and modifies delta
std::pair<std::array<gf2::Word, gf2::BLOCK>, std::array<gf2::Word, gf2::BLOCK>>
generate_P(const std::array<gf2::Word, gf2::HALF_BLOCK> &E_con_RM, 
           std::array<size_t, gf2::BLOCK> &delta) {
    // We do Gaussian elimination on the columns to find P, 
    // so everything is column major
    std::array<gf2::Word, gf2::BLOCK> P_init_CM{};
    std::array<gf2::HalfWord, gf2::BLOCK> E_con_CM{};

    // Convert E_con_RM, the "constant term" that we wish 
    // to annihilate, to column major
    for (uint32_t r = 0; r < gf2::HALF_BLOCK; r++) {
        gf2::HalfWord r_indicator = gf2::HalfWord(1) << r;
        gf2::Word c_indicator = 1;
        for (uint32_t c = 0; c < gf2::BLOCK; c++, c_indicator <<= 1) {
            if (E_con_RM[r] & c_indicator) {
                E_con_CM[c] |= r_indicator;
            }
        }
    }

    // P starts at identity 
    for (uint32_t c = 0; c < gf2::BLOCK; c++) {
        P_init_CM[c] |= gf2::Word(1) << c;
    }

    // Bubble sort the deltas, this is slow but delta.size() == gf2::BLOCK == 64, so 
    // it doesn't matter. Also, we expect the deltas to be close to uniform.
    bool sorted = false;
    while (!sorted) {
        sorted = true;
        for (uint32_t c = 0; c < gf2::BLOCK - 1; c++) {
            if (delta[c] > delta[c + 1]) {
                std::swap(delta[c], delta[c + 1]);
                std::swap(P_init_CM[c], P_init_CM[c + 1]);
                std::swap(E_con_CM[c], E_con_CM[c + 1]);
                sorted = false;
            }
        }
    }
    
    // Modified Gaussian elimination with only subtracting to later columns
    std::array<bool, gf2::BLOCK> has_pivot{};
    for (uint32_t r = 0; r < gf2::HALF_BLOCK; r++) {
        gf2::HalfWord r_indicator = gf2::HalfWord(1) << r;
        int32_t pivot = -1;
        for (uint32_t c = 0; c < gf2::BLOCK; c++) {
            if ((E_con_CM[c] & r_indicator) && !has_pivot[c]) {
                pivot = c;
                has_pivot[c] = true;
                break;
            }
        }
        if (pivot == -1) continue;
        for (uint32_t c2 = pivot + 1; c2 < gf2::BLOCK; c2++) {
            if (E_con_CM[c2] & r_indicator) {
                E_con_CM[c2] ^= E_con_CM[pivot];
                P_init_CM[c2] ^= P_init_CM[pivot];
            }
        }
    }
    
    // Columns with pivot need to be multiplied by X
    // We also transpose here to get the final P
    std::array<gf2::Word, gf2::BLOCK> P_final_con_RM{}, P_final_lin_RM{};
    for (uint32_t r = 0; r < gf2::BLOCK; r++) {
        gf2::Word r_indicator = gf2::Word(1) << r;
        gf2::Word c_indicator = 1;
        for (uint32_t c = 0; c < gf2::BLOCK; c++, c_indicator <<= 1) {
            if (P_init_CM[c] & r_indicator) {
                if (has_pivot[c]) P_final_lin_RM[r] |= c_indicator;
                else              P_final_con_RM[r] |= c_indicator;
            }
        }
    }
    // Multiplying by X incurs a delta increment
    for (uint32_t c = 0; c < gf2::BLOCK; c++) {
        if (has_pivot[c]) delta[c]++;
    }

    return {P_final_con_RM, P_final_lin_RM};
}

// Block Wiedemann, see the excellent explication here:
// https://members.loria.fr/EThome/files/fastbw.pdf
//
// Assumes cols >= rows, and if sparse_matrix is square, that it is singular
std::vector<gf2::HalfWord> gf2::sparse_nullspace(gf2::SparseMatrix sparse_matrix) {
    if (sparse_matrix.cols < sparse_matrix.rows) {
        throw std::invalid_argument("Matrix must have at least as many columns as rows");
    }

    // Make it square by padding
    sparse_matrix.rows = sparse_matrix.cols;
    sparse_matrix.row_data.resize(sparse_matrix.rows);

    /*for (uint32_t r = 0; r < sparse_matrix.rows; r++) {
        std::cout << sparse_matrix.row_data[r].size() << " ";
    }
    std::cout << std::endl;*/
    const uint32_t N = sparse_matrix.cols;
    constexpr uint32_t n = gf2::HALF_BLOCK;

    uint32_t L = 2 * N / n + 10;
    std::vector<gf2::HalfWord> X(N), Z(N), V(N, 0), nV(N, 0);

    std::random_device dev;
    static_assert(n == 32, "Randomization requires 32 bit");
    std::mt19937 rng(dev()); 
    std::uniform_int_distribution<std::mt19937::result_type> dist;

    for (uint32_t i = 0; i < N; i++) {
        X[i] = dist(rng);
        Z[i] = dist(rng);
    }

    half_word_sparse_mul_vec(V, Z, sparse_matrix); // v = B z
    gf2::BlockPoly<gf2::HalfWord, gf2::HALF_BLOCK> A;
    A.data.reserve(L);
    for (uint32_t k = 0; k < L; k++) {
        A.data.push_back(half_word_vec_dot_vec(X, V)); // A[k] = x^T v

        // v = B v
        half_word_sparse_mul_vec(nV, V, sparse_matrix);
        std::swap(nV, V);
        std::fill(nV.begin(), nV.end(), 0);
    }

    std::array<gf2::HalfWord, gf2::HALF_BLOCK> F_first_block, AF_first_block_lin;
    // Generate until 
    do {
        for (uint32_t i = 0; i < gf2::HALF_BLOCK; i++) {
            F_first_block[i] = dist(rng);
        }
       
        AF_first_block_lin = half_word_sq_mul_sq(A.data[1], F_first_block);
    } while (gf2::is_singular(AF_first_block_lin));

    // Initialize the matrices 
    std::array<gf2::Word, gf2::HALF_BLOCK> F_con, F_lin;
    for (uint32_t i = 0; i < gf2::HALF_BLOCK; i++) {
        F_con[i] = gf2::Word(F_first_block[i]);
        // Second block should be the identity matrix times X
        F_lin[i] = gf2::Word(1) << (32 + i);
    }

    // E = AF - G, where G is some constant s.t. X | E
    gf2::BlockPoly<gf2::Word, gf2::HALF_BLOCK> E = A.mul(F_con, F_lin); 
    std::fill(E.data[0].begin(), E.data[0].end(), 0);
    //E.data.erase(E.data.begin());

    gf2::BlockPoly<gf2::Word, gf2::HALF_BLOCK> F;
    F.data.push_back(F_con);
    F.data.push_back(F_lin);

    // Degree bound for block Wiedemann
    std::array<size_t, gf2::BLOCK> delta;
    std::fill(delta.begin(), delta.end(), 1);
    size_t t = 1;

    // Main iteration
    double gen_p_time = 0.0, mul_p_time = 0.0;
    for (;t < L ; t++) {
        std::chrono::system_clock::time_point tStart = std::chrono::system_clock::now();

        std::pair<std::array<gf2::Word, gf2::BLOCK>, std::array<gf2::Word, gf2::BLOCK>> 
            P = generate_P(E.data[t], delta); 

        gen_p_time += std::chrono::duration_cast<std::chrono::microseconds>
                     (std::chrono::system_clock::now() - tStart).count() / 1'000'000.0;
        tStart = std::chrono::system_clock::now();

        // Now that P has been found, 
        // multiply AF = G + E (i.e. F and E) by P on the right,  
        // ignore G (not stored)
        F.mul_eq(P.first, P.second);
        E.mul_eq(P.first, P.second);
        mul_p_time += std::chrono::duration_cast<std::chrono::microseconds>
                     (std::chrono::system_clock::now() - tStart).count() / 1'000'000.0;
    }
    std::cout << std::fixed << std::setprecision(3)
              << "gen P time: " << gen_p_time << std::endl
              << "mul P time: " << mul_p_time << std::endl;


    // If t - delta[c] > N / n, then column c of F is a
    // linear generator, which may produce a kernel vector
    std::vector<uint32_t> lin_gen_c;
    size_t max_delta = 0;
    for (uint32_t c = 0; c < gf2::BLOCK; c++) {
        if (n * t > n * delta[c] + N) {
            max_delta = std::max(max_delta, delta[c]);
            lin_gen_c.push_back(c);
            if (lin_gen_c.size() == gf2::HALF_BLOCK) break;
        }
    }

    gf2::BlockPoly<gf2::HalfWord, gf2::HALF_BLOCK> F_to_check;
    std::array<size_t, gf2::HALF_BLOCK> deg_to_check{};
    F_to_check.data.reserve(F.data.size());
    // Compress the linear generators together into a half word
    for (size_t deg = 0; deg < F.data.size(); deg++) {
        const std::array<gf2::Word, gf2::HALF_BLOCK> &F_coeff = F.data[deg];
        std::array<gf2::HalfWord, gf2::HALF_BLOCK> F_to_check_coeff{};
        gf2::HalfWord c_indicator = 1;
        for (uint32_t idx = 0; idx < lin_gen_c.size(); idx++, c_indicator <<= 1) {
            const uint32_t c = lin_gen_c[idx];
            for (uint32_t r = 0; r < gf2::HALF_BLOCK; r++) {
                if (F_coeff[r] & (gf2::Word(1) << c)) {
                    F_to_check_coeff[r] |= c_indicator;
                    deg_to_check[idx] = std::max(deg_to_check[idx], deg);
                }
            }
        }
        F_to_check.data.push_back(F_to_check_coeff);
    }

    size_t max_deg = *std::max_element(deg_to_check.begin(), deg_to_check.end());

    /*for (uint32_t idx = 0; idx < lin_gen_c.size(); idx++) {
        std::cout << deg_to_check[idx] << " ";
    }
    std::cout << std::endl;
    std::cout << max_deg << " " << max_delta << std::endl;*/

    // For each linear generator column of F_to_check, 
    // we generate a vector (which corresponds to a column of W).
    // Then, we repeatedly left multiply W by sparse_matrix; 
    // if a column turns from nonzero to zero, then we have a kernel vector
    std::vector<gf2::HalfWord> W(N, 0), U(N, 0), res(N, 0); // w, u = 0, 0
    std::vector<bool> fail(lin_gen_c.size(), false); // whether the generated vector is zero (which is a failure)

    for (size_t deg = 0; deg <= max_deg + max_delta + 10; deg++) {
        half_word_sparse_mul_vec(U, W, sparse_matrix); // u = B w
        if (deg <= max_deg) {
            half_word_vec_mul_sq(U, Z, F_to_check.data[deg]); // u = B w + z F_to_check[deg]
        }

        gf2::HalfWord not_identically_zero = 0; // Whether a column of U is not identically zero
        for (uint32_t r = 0; r < N; r++) {
            not_identically_zero |= U[r];
        }
        //std::cout << std::bitset<gf2::HALF_BLOCK>(not_identically_zero) << std::endl;

        gf2::HalfWord c_indicator = 1;
        
        for (uint32_t idx = 0; idx < lin_gen_c.size(); idx++, c_indicator <<= 1) {
            // This column has reached the deg to stop generating the vector
            if (deg == deg_to_check[idx]) {
                fail[idx] = !(not_identically_zero & c_indicator);
            }
            // this means the column of F[deg] is 
            if (deg > deg_to_check[idx]) {
                // If this column of U is 0, then the corresponding column of W
                // is the desired kernel vector, so copy it to res
                if (!fail[idx] && !(not_identically_zero & c_indicator)) {
                    for (uint32_t r = 0; r < N; r++) {
                        if (W[r] & c_indicator) res[r] |= c_indicator;
                    }
                }
            }
        }
        std::swap(U, W);
        std::fill(U.begin(), U.end(), 0);
    }
    return res;
}

