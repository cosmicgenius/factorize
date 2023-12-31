// test.cpp

#include "../include/util.hpp"
#include "../include/gf2.hpp"
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    clock_t tStart = clock();

    assert(util::primes_less_than(100).size() == 25);
    assert(util::primes_less_than(10'000).size() == 1229);
    assert(util::primes_less_than(10'000'000).size() == 664'579);
    assert(util::primes_less_than(100'000'000).size() == 5'761'455);

    std::cout << "util::primes_less_than: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;
    tStart = clock();

    uint32_t primes[] = {11, 101, 1009, 10007, 100003, 1000003};
    
    for (const uint32_t& p : primes) {
        for (uint32_t x = 0; x < p; x++) {
            uint64_t x64 = x;
            uint32_t y = (x64 * x64) % p;
            uint32_t y_rt = util::square_root_modulo_prime(y, p);
            assert(y_rt == x || y_rt + x == p);
        }
    }
    std::cout << "util::square_root_modulo_prime: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;
    tStart = clock();

    uint32_t N = 1 << 12;
    std::vector<std::vector<gf2::Word>> data(N);
    srand (1);
    for (uint32_t i = 0; i < N; i++) {
        std::vector<gf2::Word> &row = data[i];
        row = std::vector<gf2::Word>(N / gf2::BLOCK);

        for (gf2::Word &word : row) {
            word = gf2::Word(rand()) * gf2::Word(0xAC5AAB21);
        }
    }
    gf2::Matrix matrix1{N, N, data};
    gf2::Matrix matrix2 = gf2::transpose(matrix1);
    gf2::Matrix matrix3 = gf2::transpose(matrix2);

    /*std::vector<gf2::Matrix*> matrices = {&matrix1, &matrix2, &matrix3};
    int idx = 0;
    for (gf2::Matrix* matrix : matrices) {
        std::cout << "matrix" << ++idx << std::endl;
        for (uint32_t i = 0; i < N; i++) {
            for (uint32_t j = 0; j < N; j++) {
                std::cout << bool(matrix->data[j][i / gf2::BLOCK] & (gf2::Word(1) << (i % gf2::BLOCK)));
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }*/

    for (uint32_t i = 0; i < N; i++) {
        for (uint32_t j = 0; j < N / gf2::BLOCK; j++) {
            assert(matrix1.data[i][j] == matrix3.data[i][j]);
        }
    }
    std::cout << "gf2::transpose: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;

    const gf2::Matrix kernel = gf2::nullspace(matrix1);
    for (const std::vector<gf2::Word> &vec : kernel.data) {
        std::vector<gf2::Word> sum(N / gf2::BLOCK, 0);
        for (uint32_t c = 0; c < kernel.rows; c++) {
            if (vec[c / gf2::BLOCK] & (gf2::Word(1) << (c % gf2::BLOCK))) {
                for (uint32_t rq = 0; rq < N / gf2::BLOCK; rq++) {
                    sum[rq] ^= matrix1.data[c][rq];
                }   
            }
        }
        
        for (uint32_t rq = 0; rq < N / gf2::BLOCK; rq++) {
            assert(sum[rq] == 0);
        }   
    }

    std::cout << "gf2::nullspace: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;
}
