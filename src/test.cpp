// test.cpp

#include "../include/util.hpp"
#include "../include/gf2.hpp"
#include <cassert>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>

struct Four {
    uint32_t a;
    uint32_t b;
    uint32_t c;
    uint32_t d;
};

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

    // See src/sieve_handler.cpp for explanation of how the magic nums are chosen

    tStart = clock();
    /*std::vector<uint32_t> primes2{3, 5, 7, 11, 13, 17, 19, 23, 29, 
        101, 211, 1009, 2003, 10007, 100003, 1000003};*/
    std::vector<uint32_t> primes2 = util::primes_less_than(1'300'000);
    primes2.erase(primes2.begin()); // don't want the prime 2
    std::vector<std::pair<uint32_t, uint32_t>> magic_num;
    
    for (const uint32_t p : primes2) {
        uint32_t k = util::modular_inv_mod_prime<int32_t>((static_cast<uint64_t>(1) << 32) % p, p);
        uint32_t m = ((static_cast<uint64_t>(k) << 32) - 1) / p;
        uint32_t b = ((static_cast<uint64_t>(1) << 32) - 1) / p;
        magic_num.emplace_back(-m, b);
    }
    std::cout << "fdiv preproc: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;

    auto take = [] (uint32_t n) { return n % 74 == 38; };
    
    constexpr int MAXN = 1 << 18;
    constexpr int CACHE_NUM = 64;

    tStart = clock();
    int num_div_1 = 0;
    time_t now = time(NULL);
    srand(now);
    for (int i = 0; i < MAXN; i++) {
        uint32_t n = rand();
        if (!take(n)) continue;
        for (const uint32_t p : primes2) {
            if (n % p == 0) num_div_1++;
        }
    }
    std::cout << "normal div: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;

    tStart = clock();
    int num_div_2 = 0;
    srand(now);
    uint32_t n[CACHE_NUM];
    int num_n = 0;
    while (num_n < MAXN) {
        int num_valid = 0;
        while (num_valid < CACHE_NUM && num_n < MAXN) {
            int nv = rand();
            if (take(nv)) n[num_valid++] = nv;
            num_n++;
        }
        for (const uint32_t p : primes2) {
            for (int j = 0; j < num_valid; j++) {
                if (n[j] % p == 0) num_div_2++;
            }
        }
    }
    std::cout << "cache normal div: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;

    tStart = clock();
    int num_div_3 = 0;
    srand(now);
    for (int i = 0; i < MAXN; i++) {
        uint32_t n = rand();
        if (!take(n)) continue;
        for (int j = 0, sz = primes2.size(); j < sz; j++) {
            if (n * magic_num[j].first <= magic_num[j].second) {
                num_div_3++;
            }
        }
    }

    std::cout << "fdiv: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;

    tStart = clock();
    int num_div_4 = 0;
    srand(now);
    num_n = 0;
    while (num_n < MAXN) {
        int num_valid = 0;
        while (num_valid < CACHE_NUM && num_n < MAXN) {
            int nv = rand();
            if (take(nv)) n[num_valid++] = nv;
            num_n++;
        }
        for (const std::pair<uint32_t, uint32_t> &mb : magic_num) {
            for (int j = 0; j < num_valid; j++) {
                if (n[j] * mb.first <= mb.second) num_div_4++;
            }
        }
    }
    assert(num_div_1 == num_div_2);
    assert(num_div_2 == num_div_3);
    assert(num_div_3 == num_div_4);
    std::cout << "cache fdiv: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;

    now = time(NULL);
    tStart = clock();

    srand(now);
    constexpr int M = 100'000, rep = 1'000;
    uint32_t tot1 = 0;
    std::vector<uint32_t> a(M), b(M), c(M), d(M);
    //uint32_t a[M], b[M], c[M], d[M];
    for (int i = 0; i < M; i++) {
        a[i] = rand();
        b[i] = rand();
        c[i] = rand();
        d[i] = rand();
    }

    for (int r = 0; r < rep; r++) {
        for (int i = 0; i < M; i++) {
            tot1 ^= (a[i] % b[i]) * (c[i] % d[i]);
            a[i]++;
            c[i]++;
        }
    }
    std::cout << "separate: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;
    tStart = clock();

    srand(now);
    uint32_t tot2 = 0;
    std::vector<Four> abcd(M);
    //Four abcd[M];
    for (int i = 0; i < M; i++) {
        abcd[i] = {static_cast<uint32_t>(rand()),
            static_cast<uint32_t>(rand()),
            static_cast<uint32_t>(rand()),
            static_cast<uint32_t>(rand())};
    }

    for (int r = 0; r < rep; r++) {
        for (Four &f : abcd) {
            tot2 ^= (f.a % f.b) * (f.c % f.d);
            f.a++;
            f.c++;
        }
    }
    std::cout << "pair: " << std::fixed << std::setprecision(3)
              << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
              << std::endl;
    assert(tot1 == tot2);
}
