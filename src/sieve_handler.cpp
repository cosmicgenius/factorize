// sieve_handler.cpp

#include "../include/util.hpp"
#include "../include/siever.hpp"
#include "../include/gf2.hpp"
#include "../include/sieve_handler.hpp"
#include <cmath>
#include <cstdint>
#include <gmpxx.h>
#include <iostream>
#include <vector>

const double LOG2 = 0.69314718056;

// heuristically obtained 
const double BASE_SIZE_MULTIPLIER = 0.105;
const double SIEVE_RADIUS_POWER = 1.5;
const double SIEVE_RADIUS_MULTIPLIER = 0.7;

const uint32_t A_FACTOR_TARGET = 3000;

//const uint32_t THREADS = 1;

SieveHandler::SieveHandler(mpz_class N) : N_(N) {}

void SieveHandler::InitHeuristics() {
    // heuristically obtained approximate base size
    double ln_n = std::log(this->N_.get_d());
    double approx_base_size = exp(sqrt(ln_n * std::log(ln_n) * BASE_SIZE_MULTIPLIER));
    
    // we want to (Eratosthenes) sieve 
    // until we have around 2 * approx_base_size primes 
    // so that we get around approx_base_size primes in the factor base
    // that have n as a qr (asymptotically half by Chebatorev)
    //
    // this bound comes from a prime number theorem bound
    uint32_t prime_bound = (uint32_t)(2 * approx_base_size *
                 (std::log(2 * approx_base_size) + std::log(std::log(2 * approx_base_size)) - 1));

    this->all_small_primes_ = util::primes_less_than(prime_bound);
    
    this->factor_base_.clear();
    for (const PrimeSize &p : this->all_small_primes_) {
        if (p == 2 || mpz_kronecker_ui(this->N_.get_mpz_t(), p) == 1) {
            factor_base_.push_back(p);
            //std::cout << p << " ";
        }
    }
    //std::cout << std::endl;

    this->base_size_ = (uint32_t)this->factor_base_.size();
    this->sieve_radius_ = (uint32_t)(pow(this->base_size_, SIEVE_RADIUS_POWER) * SIEVE_RADIUS_MULTIPLIER);
    this->partial_prime_bound_ = (uint32_t)(100 * prime_bound);

    this->result_target_ = this->base_size_ + uint32_t(ln_n);

    this->total_sieved_ = 0;
}

PrimeSize SieveHandler::CheckSmallPrimes() {
    for (const PrimeSize &p : this->all_small_primes_) {
        if (this->N_ % p == 0) return p;
    }
    return 1;
}

void SieveHandler::PrecomputePrimeFunctions() {
    this->fb_nsqrt_ = std::vector<PrimeSize> (this->base_size_);
    this->fb_logp_ = std::vector<LogType> (this->base_size_);

    for (uint32_t i = 0; i < this->base_size_; i++) {
        const PrimeSize &p = this->factor_base_[i];
        this->fb_logp_[i] = LogType(std::log(double(p)) * LOG_PRECISION);

        try {
          this->fb_nsqrt_[i] = util::square_root_modulo_prime(this->N_, p);
        } catch (const std::domain_error &err) {
          continue;
        }
    }
}

void SieveHandler::InitSievers() {
    // For siqs, we want to sieve with polynomials of the form
    // (ax + b)^2 - N where a is the product of primes in the factor base
    // (and we can make b in [0, a))
    // To make this number small for x in [-sieve_radius, sieve_radius] (sieving range)
    // we should pick a ~ sqrt(2N) / sieve_radius is optimal
    mpz_class a_target = sqrt(2 * this->N_) / this->sieve_radius_;

    // We want each critical prime to be around A_FACTOR_TARGET, so we need 
    // around log(a_target) / log(A_FACTOR_TARGET) critical primes
    // We use ceiling for safety, since if we go smaller there are more primes for us to choose
    this->num_critical_ = ceil(std::log(a_target.get_d()) / std::log(A_FACTOR_TARGET));
    this->num_noncritical_ = this->base_size_ - this->num_critical_;

    double true_a_target = std::exp(std::log(a_target.get_d()) / this->num_critical_);
    
    this->critical_fb_lower_ = std::lower_bound(this->factor_base_.begin(), this->factor_base_.end(), 
            uint32_t(true_a_target / 1.5)) - this->factor_base_.begin(),
    this->critical_fb_upper_ = std::lower_bound(this->factor_base_.begin(), this->factor_base_.end(), 
            uint32_t(true_a_target * 1.5)) - this->factor_base_.begin();

    //std::cout << true_a_target << " " << this->critical_fb_lower_ << " " << this->critical_fb_upper_ << std::endl;

    this->sievers_.reserve(1);
    this->sievers_.emplace_back(this->N_, a_target, this->base_size_,
            this->sieve_radius_, this->partial_prime_bound_,
            this->num_critical_, this->num_noncritical_, 
            this->critical_fb_lower_, this->critical_fb_upper_, 
            this->factor_base_, this->fb_nsqrt_, this->fb_logp_,
            this->total_sieved_, this->sieve_results_, this->partial_sieve_results_);
}

bool SieveHandler::Sieve() {
    for (Siever &siever : this->sievers_) {
        siever.SievePolynomialGroup();
    }

    return this->sieve_results_.size() > this->result_target_;
}

mpz_class SieveHandler::TryExtractDivisor() {
    const size_t cols = this->sieve_results_.size(), rows = this->base_size_ + 1;
    std::vector<std::vector<gf2::Word>> matrix_data;
    matrix_data.reserve(cols);

    // TODO: (repeatedly) prune this array for results that are the only 
    // ones to use a certain prime?
    for (const SieveResult &result : this->sieve_results_) {
        std::vector<gf2::Word> res_vec(rows / gf2::BLOCK + (rows % gf2::BLOCK!= 0), 0);
        if (result.sgn) res_vec[0] |= 1; // Set first bit to 1 for the prime ``-1''

        for (const uint32_t &fb_idx : result.prime_fb_idxs) {
            const uint32_t q = (fb_idx + 1) / gf2::BLOCK, r = (fb_idx + 1) % gf2::BLOCK;
            res_vec[q] ^= (gf2::Word(1) << r);
        }
        matrix_data.push_back(res_vec);
    }
    
    gf2::Matrix kernel = gf2::nullspace(gf2::Matrix{cols, rows, matrix_data});
    for (const std::vector<gf2::Word> &vec : kernel.data) {
        // lhs is the value that comes as a product of squares
        // rhs is the square produced from multipyling smooth values
        // TODO: Is mod N faster?
        mpz_class lhs_numer_rt = 1, lhs_denom_rt = 1, rhs_rt = 1;
        bool rhs_sgn = false;
        std::vector<uint32_t> prime_exp(this->base_size_, 0);

        for (uint32_t res_idx = 0; res_idx < cols; res_idx++) {
            if (vec[res_idx / gf2::BLOCK] & (gf2::Word(1) << (res_idx % gf2::BLOCK))) {
                // Result should be a part of the product
                const SieveResult &result = this->sieve_results_[res_idx];

                lhs_numer_rt *= result.numer;
                lhs_denom_rt *= result.denom;
                rhs_sgn ^= result.sgn;
                for (const uint32_t &fb_idx : result.prime_fb_idxs) {
                    prime_exp[fb_idx]++;
                }
            }
        }
        //if (rhs_sgn) std::cout << "BAD sign" << std::endl;
        for (uint32_t fb_idx = 0; fb_idx < this->base_size_; fb_idx++) {
            //if (prime_exp[fb_idx] % 2 != 0) std::cout << "BAD " << fb_idx << std::endl;
            rhs_rt *= util::pow(this->factor_base_[fb_idx], prime_exp[fb_idx] / 2);
        }

        // It seems like this is a worthwhile mod N, but maybe not
        lhs_denom_rt %= this->N_;
        rhs_rt %= this->N_;

        mpz_class g = gcd(abs(lhs_numer_rt - lhs_denom_rt * rhs_rt), this->N_);

        if (g > 1 && g < this->N_) {
            return g;
        }
    }
    return 1;
}

uint32_t SieveHandler::get_base_size_() const { return this->base_size_; }
uint32_t SieveHandler::get_sieve_radius_() const { return this->sieve_radius_; }
uint32_t SieveHandler::get_partial_prime_bound_() const { return this->partial_prime_bound_; }

uint32_t SieveHandler::get_results_target_() const { return this->result_target_; }

uint32_t SieveHandler::get_num_critical_() const { return this->num_critical_; }
uint32_t SieveHandler::get_num_noncritical_() const { return this->num_noncritical_; }

size_t SieveHandler::get_sieve_results_size_() const { return this->sieve_results_.size(); }
size_t SieveHandler::get_partial_sieve_results_size_() const { return this->partial_sieve_results_.size(); }

