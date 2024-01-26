// sieve_handler.hpp
#ifndef SIEVE_HANDLER_HPP_
#define SIEVE_HANDLER_HPP_

#include "siever.hpp"
#include "gf2.hpp"
#include <cstdint>
#include <functional>
#include <gmpxx.h>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

class SieveHandler {
private:
    // true N is the N to be factored, N is the value 
    // after multiplier
    mpz_class N_, true_N_ = 0;
    uint32_t base_size_; 
    uint32_t sieve_radius_;
    uint32_t partial_prime_bound_; 
    uint32_t result_target_; 

    std::vector<PrimeSize> all_small_primes_;
    std::vector<PrimeSize> factor_base_;
    // Store square root of N mod p and log(p) for each p based on fb index
    std::vector<PrimeSize> fb_nsqrt_;
    std::vector<LogType> fb_logp_;
    std::vector<std::pair<uint32_t, uint32_t>> fb_magic_num_p_;

    uint32_t num_critical_;
    //uint32_t num_noncritical_;
    
    uint32_t critical_fb_lower_;
    uint32_t critical_fb_upper_;

    uint32_t total_sieved_;
    
    std::vector<SieveResult> sieve_results_;
    std::unordered_map<uint32_t, SieveResult> partial_sieve_results_;

    std::mutex res_mutex_;

    // Columns correspond to results,
    // rows correspond to primes, with the entry
    // corresponding to whether the prime appears an odd 
    // number of times in the prime factorization of the 
    // smooth result
    gf2::Matrix results_matrix_;

    // Matrix indices after pruning
    //std::vector<uint32_t> relevant_fb_idxs_; // does not actually need to be stored 
    //                                            because the SieveResult struct stores 
    //                                            actual fb indices
    std::vector<uint32_t> relevant_res_idxs_;

    // Lists of sievers and threads where 
    // thread i runs siever i
    std::vector<Siever> sievers_;
    std::vector<std::thread> threads_;

    uint32_t polygroup_ = 0;
    Timer timer_;

public:
    explicit SieveHandler(mpz_class N);

    void InitHeuristics();
    
    // Guard against primes in the factor base 
    // dividing N by checking and returning early if so
    PrimeSize CheckSmallPrimes();

    void PrecomputePrimeFunctions();
    void InitSievers();

    // Sieves until there are enough results
    // Takes a void(uint32_t) function that runs with argument polygroup_ 
    // every time a polygroup finishes sieving 
    // (relevant for multithreaded sieving)
    void Sieve(const std::function<void(uint32_t)> &on_finish_polygrp);

    void GenerateMatrix();
    mpz_class TryExtractDivisor();

    uint32_t get_KS_multiplier_() const;

    uint32_t get_base_size_() const;
    uint32_t get_sieve_radius_() const;
    uint32_t get_partial_prime_bound_() const;

    uint32_t get_results_target_() const;

    uint32_t get_num_critical_() const;
    //uint32_t get_num_noncritical_() const;

    size_t get_sieve_results_size_() const;
    size_t get_partial_sieve_results_size_() const;

    std::pair<size_t, size_t> get_matrix_dim() const;

    Timer get_timer_() const;
};

#endif // SIEVE_HANDLER_HPP_
