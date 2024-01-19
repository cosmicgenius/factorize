// siever.cpp

#include "../include/util.hpp"
#include "../include/siever.hpp"
#include <algorithm>
#include <bits/chrono.h>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <ostream>
#include <unordered_map>
#include <vector>


// See https://www.mersenneforum.org/showthread.php?t=7212#post99513
// the idea is that tiny primes are quite wasteful, and contribute mostly uniform heights
// (~ TINY_PRIME_PADDING) at least to smooth values
// and thus can be ignored
//
// Ok turns out this is mostly useless because what it really does 
// is a tiny SetHeights speed up and combined with more large prime forgiveness
const PrimeSize TINY_PRIME_BOUND = 2;
const double TINY_PRIME_PADDING = 2;

// TODO: does this denom optimization even help?
SieveResult::SieveResult(const mpz_class &&numer, const mpz_class &&denom, 
        const bool sgn, const std::vector<uint32_t> &&prime_fb_idxs) 
    : numer(std::move(numer)), denom(std::move(denom)), sgn(sgn), prime_fb_idxs(std::move(prime_fb_idxs)) {}

// TODO: does this even help?
SieveResult& SieveResult::operator=(SieveResult &&rhs) {
    if (this != &rhs) {
        this->numer = std::move(rhs.numer);
        this->denom = std::move(rhs.denom);
        this->sgn = rhs.sgn;
        this->prime_fb_idxs = std::move(rhs.prime_fb_idxs);
    }
    return *this;
}
SieveResult::SieveResult(SieveResult &&rhs) noexcept {
    *this = std::move(rhs);
}

Siever::Siever(const mpz_class &N, const mpz_class &a_target, 
        const uint32_t &base_size, const uint32_t &sieve_radius, const uint32_t &partial_prime_bound,
        const uint32_t &num_critical, const uint32_t &num_noncritical, 
        const uint32_t &critical_fb_lower, const uint32_t &critical_fb_upper, 
        const std::vector<PrimeSize> &factor_base,
        const std::vector<PrimeSize> &fb_nsqrt, const std::vector<LogType> &fb_logp,
        uint32_t &total_sieved,
        std::vector<SieveResult> &sieve_results, std::unordered_map<uint32_t, SieveResult> &partial_sieve_results,
        Timer &timer) :
        N_(N), a_target_(a_target), 
        base_size_(base_size), sieve_radius_(sieve_radius), partial_prime_bound_(partial_prime_bound),
        num_critical_(num_critical), num_noncritical_(num_noncritical),
        critical_fb_lower_(critical_fb_lower), critical_fb_upper_(critical_fb_upper),
        factor_base_(factor_base), fb_nsqrt_(fb_nsqrt), fb_logp_(fb_logp),
        total_sieved_(total_sieved), 
        sieve_results_(sieve_results), partial_sieve_results_(partial_sieve_results),
        timer_(timer) {
    double a_d = a_target.get_d(),
           N_over_a_d = mpz_class(N / a_target).get_d();

    this->sieve_height_ = std::vector<LogType>(2 * sieve_radius + 1);

    // Populate sieve_height_target_ with log|a (x - sieve_radius) ** 2 - N / a| - log(partial_prime_bound),
    // which is the approximate log of the value of the polynomial.
    this->sieve_height_target_ = std::vector<LogType>(2 * sieve_radius + 1, 0);

    double log_partial_prime_bound = std::log(double(partial_prime_bound));
    for (uint32_t x = 0; x <= 2 * sieve_radius; x++) {
        double logpoly = std::log(std::abs(
                    a_d * (double(x) - sieve_radius) * (double(x) - sieve_radius) - N_over_a_d
                ));
        if (logpoly > log_partial_prime_bound + TINY_PRIME_PADDING) {
            this->sieve_height_target_[x] = LogType(
                    (logpoly - log_partial_prime_bound - TINY_PRIME_PADDING)
                * LOG_PRECISION);
        }
    }
}

void Siever::InitCriticalPrimes() {
    // We find critical primes through a method described here:
    // https://www.mersenneforum.org/showthread.php?p=535652
    //
    // Essentially we randomly pick num_critical - 1 primes in the prime base in the correct range,
    // then hand pick the last prime to be as close as possible
    //
    // Technically, this method ever so slightly skews the last prime to be
    // larger due to 
    // 1. the distribution of primes, and 
    // 2. the distribution of rand, 
    // but it probably doesn't matter.

    this->a_ = 1;
    this->critical_fb_idxs_.clear();

    for (uint32_t _i = 0; _i < this->num_critical_ - 1; _i++) {
        uint32_t picked_fb_idx = 0;

        // if there are repeats, throw out and keep selecting
        do {
            picked_fb_idx = rand() % (this->critical_fb_upper_ - this->critical_fb_lower_) 
                         + this->critical_fb_lower_;
        } while (this->critical_fb_idxs_.find(picked_fb_idx) !=
                 this->critical_fb_idxs_.end());

        this->critical_fb_idxs_.insert(picked_fb_idx);
        this->a_ *= this->factor_base_[picked_fb_idx];
    }

    // Optimal size of the last prime based on the target size of a,
    // and the primes picked so far
    PrimeSize last_prime_predict = mpz_class(this->a_target_ / this->a_).get_ui();
    uint32_t last_prime_fb_idxs = std::lower_bound(this->factor_base_.begin(), 
            this->factor_base_.end(), last_prime_predict) - this->factor_base_.begin();

    for (uint32_t fb_idx = last_prime_fb_idxs; fb_idx < this->factor_base_.size(); fb_idx++) {
        if (this->critical_fb_idxs_.find(fb_idx) == this->critical_fb_idxs_.end()) {
            this->critical_fb_idxs_.insert(fb_idx);
            this->a_ *= this->factor_base_[fb_idx];
            goto done;
        }
    }
    for (int32_t fb_idx = last_prime_fb_idxs - 1; fb_idx >= 0; fb_idx--) {
        if (this->critical_fb_idxs_.find(fb_idx) == this->critical_fb_idxs_.end()) {
            this->critical_fb_idxs_.insert(fb_idx);
            this->a_ *= this->factor_base_[fb_idx];
            goto done;
        }
    }
done:
    this->noncritical_fb_idxs_.clear();
    this->noncritical_fb_idxs_.reserve(this->num_noncritical_);
    // get noncritical indexes by testing each for membership in critical_fb_idxs
    // this is slow but num_critical ~ 10, so it doesn't really matter
    for (uint32_t fb_idx = 0; fb_idx < this->base_size_; fb_idx++) {
        if (this->factor_base_[fb_idx] > TINY_PRIME_BOUND &&
                this->critical_fb_idxs_.find(fb_idx) == this->critical_fb_idxs_.end()) {
            this->noncritical_fb_idxs_.push_back(fb_idx);
        }
    }
}

void Siever::InitCriticalDeltas() {
    uint32_t crt_idx = 0;
    this->crt_indicator_ = std::vector<mpz_class>(this->num_critical_, 0);
    for (const uint32_t &critical_fb_idx : this->critical_fb_idxs_) {
        const PrimeSize p = this->factor_base_[critical_fb_idx];
        mpz_class product_others = this->a_ / p;

        // The idea is that product_others * gamma will be a square root 
        // of N mod p, but 0 mod all other critical primes.
        uint32_t gamma = this->fb_nsqrt_[critical_fb_idx] 
                * util::modular_inv_mod_prime(mpz_class(product_others % p).get_ui(), p) % p;
        if (gamma > p / 2) gamma = p - gamma;

        mpz_class indicator = product_others * gamma;
        this->crt_indicator_[crt_idx++] = indicator;
    }

    this->a_inv_ = std::vector<PrimeSize>(this->num_noncritical_);
    this->soln_delta_ = std::vector<std::vector<PrimeSize>>(this->num_noncritical_);

    // idx indexes over noncritical primes
    for (uint32_t idx = 0; idx < this->num_noncritical_; idx++) {
        const uint32_t noncritical_fb_idx = this->noncritical_fb_idxs_[idx];
        const PrimeSize q = this->factor_base_[noncritical_fb_idx];

        this->a_inv_[idx] = 
            util::modular_inv_mod_prime(mpz_class(this->a_ % q).get_ui(), q);
        
        this->soln_delta_[idx].resize(this->num_critical_);
        for (crt_idx = 0; crt_idx < this->num_critical_; crt_idx++) {
            uint64_t indicator_64 = mpz_class(this->crt_indicator_[crt_idx] % q).get_ui();
            this->soln_delta_[idx][crt_idx] = (indicator_64 << 1) * this->a_inv_[idx] % q;
        }   
    }
}

// All idx's index over the noncritical primes
void Siever::InitNextPoly() {
    if (this->poly_ == 0) {
        this->b_ = 0;
        for (uint32_t crt_idx = 0; crt_idx < this->num_critical_; crt_idx++) {
            this->b_ += this->crt_indicator_[crt_idx];
        }
        uint32_t idx = 0;
        this->soln_ = std::vector<std::pair<PrimeSize, PrimeSize>>(this->num_noncritical_);
        for (const uint32_t noncritical_fb_idx : this->noncritical_fb_idxs_) {
            const PrimeSize q = this->factor_base_[noncritical_fb_idx];
            const mpz_class &rt = this->fb_nsqrt_[noncritical_fb_idx];
            //const LogType &logq = this->fb_logp_[noncritical_fb_idx];
            this->soln_[idx].first = mpz_class(
                    (((-rt - this->b_ % q) * this->a_inv_[idx] + this->sieve_radius_) % q + q) % q
                ).get_ui();
            this->soln_[idx].second = mpz_class(
                    (((rt - this->b_ % q) * this->a_inv_[idx] + this->sieve_radius_) % q + q) % q
                ).get_ui();
            
            idx++;
        }

        this->poly_++;
        return;
    } 

    // Position of Gray code difference bit
    uint32_t nu = 0;
#ifdef __GNUC__
    nu = __builtin_ctz(this->poly_);
#else
    uint32_t poly2 = poly;
    poly2 = poly2 | (poly2 << 1);
    poly2 = poly2 | (poly2 << 2);
    poly2 = poly2 | (poly2 << 4);
    poly2 = poly2 | (poly2 << 8);
    poly2 = poly2 | (poly2 << 16);
    nu = std::bitset<32>(~poly2).count();
#endif
    uint32_t gray_code_indicator = (this->poly_ >> (nu + 1));
    if (gray_code_indicator % 2 == 0) {
        this->b_ -= (this->crt_indicator_[nu] << 1);
    } else {
        this->b_ += (this->crt_indicator_[nu] << 1);
    }

    uint32_t idx = 0;
    for (const uint32_t noncritical_fb_idx : this->noncritical_fb_idxs_) {
        const PrimeSize q = this->factor_base_[noncritical_fb_idx];
        if (gray_code_indicator % 2 == 0) {
            this->soln_[idx].first  = 
                (this->soln_[idx].first  + this->soln_delta_[idx][nu]) % q;
            this->soln_[idx].second = 
                (this->soln_[idx].second + this->soln_delta_[idx][nu]) % q;
        } else {
            this->soln_[idx].first  = 
                (this->soln_[idx].first  + q - this->soln_delta_[idx][nu]) % q;
            this->soln_[idx].second = 
                (this->soln_[idx].second + q - this->soln_delta_[idx][nu]) % q;
        }
        idx++;
    }
    this->poly_++;
}

void Siever::SetHeights() {
    std::fill(this->sieve_height_.begin(), this->sieve_height_.end(), 0);
    uint32_t sieve_diameter = 2 * this->sieve_radius_;

    // Sieve
    for (uint32_t idx = 0; idx < this->num_noncritical_; idx++) {
        const uint32_t noncritical_fb_idx = this->noncritical_fb_idxs_[idx];
        const uint32_t q = this->factor_base_[noncritical_fb_idx];
        const LogType &logq = this->fb_logp_[noncritical_fb_idx];

        // We assume that q = 2 is taken care of by the tiny prime bound
        const std::pair<PrimeSize, PrimeSize> &solns = this->soln_[idx];

        for (uint32_t x = solns.first; x <= sieve_diameter; x += q) {
            this->sieve_height_[x] += logq;
        }
        for (uint32_t x = solns.second; x <= sieve_diameter; x += q) {
            this->sieve_height_[x] += logq;
        }
    }
}

void Siever::CheckHeights() {
    uint32_t sieve_diameter = 2 * this->sieve_radius_;

    // Check results
    // We have (ax+b) ** 2 - N = a(a x ** 2 + 2b x + c) for the integer c = (b ** 2 - N) / a
    // We also calculate rt = ax + b for extracting the final square equality mod N
    mpz_class c = (this->b_ * this->b_ - this->N_) / this->a_;
    for (uint32_t x = 0; x <= sieve_diameter; x++) {
        if (this->sieve_height_[x] > this->sieve_height_target_[x]) {
            // list of fb indices for primes dividing (ax+b) ** 2 - N, so by default, 
            // it should contain all of the critical primes 
            std::vector<uint32_t> prime_fb_idxs(critical_fb_idxs_.begin(), critical_fb_idxs_.end());

            int32_t x_int = int32_t(x) - this->sieve_radius_;
            mpz_class rt = this->a_ * x_int + this->b_;
            mpz_class polyval = this->a_ * x_int * x_int + this->b_ * x_int * 2 + c;
            bool sgn = polyval < 0;

            CheckSmoothness(x, polyval, prime_fb_idxs);

            if (polyval == 1) {
                this->sieve_results_.emplace_back(std::move(rt), 1, sgn, std::move(prime_fb_idxs));
            } else if (polyval.fits_uint_p()
                    && polyval.get_ui() <= this->partial_prime_bound_) {

                // We don't check this, but all partials are going to be primes (or their negations),
                // simply because our partial_prime_bound_ is not that much (only like 100x) 
                // bigger than the factor base primes, so anything that small will either be smooth 
                // or prime
                uint32_t partial = polyval.get_si();
                
                InsertPartial(partial, sgn, rt, prime_fb_idxs);
            } 
            /*else {
                std::cout << "rejected abs(" << polyval << ") > " << this->partial_prime_bound_ << std::endl;
            }*/
        }
    }
}

void Siever::CheckSmoothness(const int32_t x, mpz_class &polyval, std::vector<uint32_t> &prime_fb_idxs) {
    // Remove the ``-1'' prime
    polyval = abs(polyval);

    // TODO: the original implementation in cosmicgenius/algorithm
    // trial divides twice, first to check for smoothness, then to 
    // actually get primes, presumably because there are too many 
    // false positives/unpaired partials and constructing the 
    // vector for these is wasteful. Is this justified?
    
    // this is slightly inefficient 
    // (should be possible to simultaneously check for divisibility
    // and get the quotient if divisible), but we are dividing 
    // like ~ 10 times, so does it really matter?
    
    // This isn't faster
    // Use raw c impl to be faster
    // mpz_t n, q, _r, one;
    // mpz_init(n);
    // mpz_init(q);
    // mpz_init(_r);
    // mpz_init(one);
    // uint32_t true_res;

    // mpz_set(n, polyval.get_mpz_t());
    // mpz_set_ui(one, 1);

    // Three separate checks for tiny primes, critical primes, and noncritical primes
    for (uint32_t fb_idx = 0; fb_idx < this->base_size_; fb_idx++) {
        const PrimeSize p = this->factor_base_[fb_idx];
        if (p > TINY_PRIME_BOUND) break;
        
        bool divisible = (polyval % p == 0);
        while (divisible) {
            polyval /= p;
            prime_fb_idxs.push_back(fb_idx);
            divisible = (polyval % p == 0);
        }
        if (polyval == 1) return;
    }

    for (const uint32_t &fb_idx : this->critical_fb_idxs_) {
        const PrimeSize p = this->factor_base_[fb_idx];

        bool divisible = (polyval % p == 0);
        while (divisible) {
            polyval /= p;
            prime_fb_idxs.push_back(fb_idx);
            divisible = (polyval % p == 0);
        }
        if (polyval == 1) return;
    }

    uint32_t idx = 0;
    for (const uint32_t &fb_idx : this->noncritical_fb_idxs_) {
        const PrimeSize p = this->factor_base_[fb_idx];
        bool divisible = (std::abs(x - int32_t(this->soln_[idx].first )) % p == 0) ||
                         (std::abs(x - int32_t(this->soln_[idx].second)) % p == 0);
        
        // Outer check for divisibility
        while (divisible) {
        // while (true_res == 0) {
            // mpz_swap(n, q); 
            polyval /= p;
            prime_fb_idxs.push_back(fb_idx);
            divisible = (polyval % p == 0);
            // true_res = mpz_fdiv_qr_ui(q, _r, n, p);
        }
        if (polyval == 1) return;
        // if (mpz_cmp(n, one) == 0) goto done_checking;
        idx++;
    }
}

void Siever::InsertPartial(const uint32_t partial, const bool sgn,
        const mpz_class &rt, const std::vector<uint32_t> &prime_fb_idxs) {
    std::unordered_map<uint32_t, SieveResult>::const_iterator partial_res_pr = 
        this->partial_sieve_results_.find(partial);
    if (partial_res_pr != this->partial_sieve_results_.end()) {
        const SieveResult& res = partial_res_pr->second;

        std::vector<PrimeSize> combined_fb_idx(res.prime_fb_idxs.begin(), res.prime_fb_idxs.end());
        combined_fb_idx.reserve(prime_fb_idxs.size() + res.prime_fb_idxs.size());

        // No point moving primitives
        combined_fb_idx.insert(combined_fb_idx.end(), prime_fb_idxs.begin(), prime_fb_idxs.end());

        // Merge the two prime lists

        // Combine the two partial results into one of the form 
        // (rt1 * rt2) ** 2 = partial ** 2 * (something smooth) mod N
        //
        // TODO: check whether it is good to mod this product by N
        this->sieve_results_.emplace_back(res.numer * rt, mpz_class(partial), 
                sgn ^ res.sgn, std::move(combined_fb_idx));

        // The other partial result (the one stored in the map) can be kept
        // because e1 + e2, e1 + e3, ..., e1 + en forms a basis for 
        // (e1 + e2 + ... + en) perp over GF(2)n 
    } else {
        this->partial_sieve_results_.emplace(partial, 
                SieveResult(std::move(rt), 1, sgn, std::move(prime_fb_idxs)));
    }
}


void Siever::SievePoly() {
    SetHeights();
    this->timer_.set_height_time += UpdateTime();
    CheckHeights();
    this->timer_.check_time += UpdateTime();
}

std::chrono::system_clock::time_point time() {
    return std::chrono::system_clock::now();
}

double Siever::UpdateTime() {
    std::chrono::system_clock::time_point time_now = time();
    double res = 
        std::chrono::duration_cast<std::chrono::microseconds>(time_now - this->time_prev_).count() 
        / 1'000'000.0;
    this->time_prev_ = time_now;

    return res;
}

void Siever::SievePolynomialGroup() {
    this->time_prev_ = time();
    InitCriticalPrimes();
    InitCriticalDeltas();
    this->timer_.init_grp_time += UpdateTime();

    // We only need half of these since we do 
    // not need both -b and b
    uint32_t num_polys = (1 << (this->num_critical_ - 1));

    for (this->poly_ = 0; this->poly_ < num_polys;) {
        //std::cout << "poly " << this->poly_ << std::endl;
        InitNextPoly();
        this->timer_.init_poly_time += UpdateTime();
        SievePoly();

        this->total_sieved_++;
    }
}
