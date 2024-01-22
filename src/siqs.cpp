// siqs.cpp

#include "../include/util.hpp"
#include "../include/sieve_handler.hpp"
#include <algorithm>
#include <bits/chrono.h>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <map>
#include <vector>

const uint32_t PRIMALITY_REPS = 40;
const mpz_class POLLARD_CUTOFF = mpz_class(1) << 64;
const uint32_t LOOK_BEHIND = 16;

mpz_class mini_pollard(mpz_class N) {
    gmp_randclass rand_state (gmp_randinit_mt);
    rand_state.seed(std::chrono::system_clock::now().time_since_epoch().count());

    mpz_class g;
    do {
        mpz_class turtle = 1 + rand_state.get_z_range(N - 2);
        mpz_class hare(turtle);
        g = 1;

        while (g == 1) {
            turtle *= turtle;
            turtle += 1;
            turtle %= N;

            hare *= hare;
            hare += 1;
            hare %= N;
            hare *= hare;
            hare += 1;
            hare %= N;

            g = gcd(abs(turtle - hare), N);
        }
    } while (g == N);
    return g;
}

mpz_class find_nontrivial_factor(const mpz_class &N) {
    if (mpz_probab_prime_p(N.get_mpz_t(), PRIMALITY_REPS)) {
        return mpz_class(1);
    }

    std::pair<mpz_class, uint32_t> pr = util::factor_perfect_power(N);
    if (pr.second != 1) {
        return pr.first;
    }

    std::cout << "Factoring " << N << std::endl; 

    if (N < POLLARD_CUTOFF) {
        mpz_class d = mini_pollard(N);
        std::cout << "Found small factor " << d << ". Returning early" << std::endl;
        return d;
    }

    SieveHandler sieve_handler(N);

    std::cout << "Initializing basic bounds" << std::endl;
    
    sieve_handler.InitHeuristics();

    std::cout << "Knuth-Shroeppel multiplier = " << sieve_handler.get_KS_multiplier_() << std::endl;
    std::cout << "Base size = " << sieve_handler.get_base_size_() << std::endl;
    std::cout << "Sieve radius = " << sieve_handler.get_sieve_radius_() << std::endl;
    std::cout << "Partial prime bound = " << sieve_handler.get_partial_prime_bound_() << std::endl;

    uint32_t small_d = sieve_handler.CheckSmallPrimes();
    if (small_d != 1) {
        std::cout << "Found small factor " << small_d << ". Returning early" << std::endl;
        return small_d;
    }

    std::cout << "Precomputing prime functions" << std::endl;
    sieve_handler.PrecomputePrimeFunctions();

    std::cout << "Initializing sievers" << std::endl;
    sieve_handler.InitSievers();

    size_t prev_num_res = 0, cur_num_res = 0;
    
    // Stores number of results from the past LOOK_BEHIND
    // poly groups to give an average number of results being
    // produced per poly group 
    std::vector<uint32_t> look_behind_num_res;
    uint32_t look_behind_diff = 0;

    mpz_class d = 1; 
    
    while (d == 1) {
        std::cout << "Sieving..." << std::endl;

        // A function that is run sequentially after 
        // any polygroup finishes 
        auto pretty_print = [&] (uint32_t polygroup) {
            cur_num_res = sieve_handler.get_sieve_results_size_();
            look_behind_diff = cur_num_res;

            look_behind_num_res.push_back(cur_num_res);
            if (look_behind_num_res.size() > LOOK_BEHIND) {
                look_behind_diff -= look_behind_num_res[0];
                look_behind_num_res.erase(look_behind_num_res.begin());
            }

            std::cout << std::fixed << std::setprecision(2);

            std::cout << "Polygrp " << (polygroup + 1) << " with "
                      << "num_critical=" << sieve_handler.get_num_critical_() << ": "
                      << cur_num_res - prev_num_res
                      << " new res, "
                      << cur_num_res
                      << " / "
                      << sieve_handler.get_results_target_()
                      << " total "
                      << "(" << 1000.0 * cur_num_res / sieve_handler.get_results_target_() / 10.0 << "%)," 
                      << " averaging " 
                      << round(1000.0 * look_behind_diff / look_behind_num_res.size()) / 1000.0
                      << " in the last few, "
                      << sieve_handler.get_partial_sieve_results_size_()
                      << " partials. " << '\r' << std::flush;
            prev_num_res = cur_num_res;
        };

        sieve_handler.Sieve(pretty_print);

        std::cout << std::endl;
        std::cout << "Stopping sieving after finding " << cur_num_res << " > " 
                  << sieve_handler.get_results_target_() << " results" << std::endl;
        std::cout << "Generating matrix" << std::endl;

        sieve_handler.GenerateMatrix();
        std::pair<size_t, size_t> matrix_dim = sieve_handler.get_matrix_dim();

        std::cout << "Generated matrix with " << matrix_dim.first << " results and " 
                  << matrix_dim.second << " primes. Trying to extract a divisor" << std::endl;

        d = sieve_handler.TryExtractDivisor();

        if (d == 1) std::cout << "Failed to find a nontrivial divisor. Restarting" << std::endl;
    }

    std::cout << "Found nontrivial divisor d=" << d << std::endl; 
    Timer timer = sieve_handler.get_timer_();
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Polygrp init time: " << timer.init_grp_time << "s" << std::endl;
    std::cout << "Single poly init time: " << timer.init_poly_time << "s" << std::endl;
    std::cout << "Set heights time: " << timer.set_height_time << "s" << std::endl;
    std::cout << "Check heights time: " << timer.check_time << "s" << std::endl;
    std::cout << "Finding kernel time: " << timer.kernel_time << "s" << std::endl;
    return d;
}

int main() {
    std::cout << "n: ";

    std::string n_str;
    std::cin >> n_str;
    mpz_class n(n_str);

    std::chrono::system_clock::time_point tStart = std::chrono::system_clock::now();
    util::print_prime_fact(n, [](const mpz_class &b) { return find_nontrivial_factor(b); });

    std::cout << "Total time taken: " 
              << std::chrono::duration_cast<std::chrono::microseconds>
                 (std::chrono::system_clock::now() - tStart).count() / 1'000'000.0
              << "s"
              << std::endl;
}

