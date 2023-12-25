// main.cpp

#include "../include/util.hpp"
#include <chrono>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <time.h>

const int PRIMALITY_REPS = 40;

mpz_class find_nontrivial_factor(const mpz_class &n) {
    if (mpz_probab_prime_p(n.get_mpz_t(), PRIMALITY_REPS)) {
        return mpz_class(1);
    }

    std::pair<mpz_class, uint32_t> p = util::factor_perfect_power(n);
    if (p.second != 1) {
        return p.first;
    }

    gmp_randclass rand_state (gmp_randinit_mt);
    rand_state.seed(std::chrono::system_clock::now().time_since_epoch().count());

    mpz_class g;
    do {
        mpz_class turtle = 1 + rand_state.get_z_range(n - 2);
        mpz_class hare(turtle);
        g = 1;

        while (g == 1) {
            turtle *= turtle;
            turtle += 1;
            turtle %= n;

            hare *= hare;
            hare += 1;
            hare %= n;
            hare *= hare;
            hare += 1;
            hare %= n;

            g = gcd(abs(turtle - hare), n);
        }
    } while (g == n);

    // std::cout << " -> " << mpz_t::gcd(g, n) << "]" << std::flush;

    return g;
}

int main() {
    while (true) {
        std::cout << "n: ";

        std::string n_str;
        std::cin >> n_str;
        mpz_class n(n_str);

        clock_t tStart = clock();
        util::print_prime_fact(n, [](const mpz_class &b) { return find_nontrivial_factor(b); });

        std::cout << "Time taken: " << std::fixed << std::setprecision(3)
                  << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s"
                  << std::endl;
    }
}
