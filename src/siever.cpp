// siever.cpp

#include "../include/util.hpp"
#include "../include/siever.hpp"
#include <algorithm>
#include <bits/chrono.h>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <iterator>
#include <mutex>
#include <ostream>
#include <unordered_map>
#include <vector>


// See https://www.mersenneforum.org/showthread.php?t=7212#post99513
// the idea is that tiny primes are quite wasteful, and contribute mostly uniform heights
// (~ TINY_PRIME_PADDING) at least to smooth values
// and thus can be ignored
const PrimeSize TINY_PRIME_BOUND = 256;
const double TINY_PRIME_PADDING = 2.0;

Timer& Timer::operator+=(const Timer &rhs) {
    this->init_grp_time += rhs.init_grp_time;
    this->init_first_poly_time += rhs.init_first_poly_time;
    this->init_next_poly_time += rhs.init_next_poly_time;
    this->set_height_time += rhs.set_height_time;
    this->wait_res_time += rhs.wait_res_time;
    this->check_time += rhs.check_time;
    this->trial_divide_time += rhs.trial_divide_time;
    this->insert_time += rhs.insert_time;
    this->flush_time += rhs.flush_time;
    this->kernel_time += rhs.kernel_time;
    return *this;
}
Timer Timer::operator+(const Timer &rhs) const {
    Timer res(*this);
    return res += rhs;
}

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

PossibleResult::PossibleResult(const int32_t x_shift, const mpz_class &&rt, const mpz_class &&polyval, 
        const bool sgn, const std::vector<uint32_t> &&prime_fb_idxs)
    : x_shift(x_shift), rt(std::move(rt)), polyval(std::move(polyval)), sgn(sgn), prime_fb_idxs(std::move(prime_fb_idxs)) {}

PossibleResult& PossibleResult::operator=(PossibleResult &&rhs) {
    if (this != &rhs) {
        this->x_shift = rhs.x_shift;
        this->rt = std::move(rhs.rt);
        this->polyval = std::move(rhs.polyval);
        this->sgn = rhs.sgn;
        this->prime_fb_idxs = std::move(rhs.prime_fb_idxs);
    }
    return *this;
}
PossibleResult::PossibleResult(PossibleResult &&rhs) noexcept {
    *this = std::move(rhs);
}

Siever::Siever(const mpz_class N, const mpz_class a_target, 
        const uint32_t base_size, const uint32_t sieve_radius, const uint32_t partial_prime_bound,
        const uint32_t num_critical, //const uint32_t &num_noncritical, 
        const uint32_t critical_fb_lower, const uint32_t critical_fb_upper, 
        const std::vector<PrimeSize> factor_base,
        const std::vector<PrimeSize> fb_nsqrt, const std::vector<LogType> fb_logp,
        const std::vector<std::pair<uint32_t, uint32_t>> fb_magic_num_p,
        uint32_t &total_sieved,
        std::vector<SieveResult> &sieve_results, std::unordered_map<uint32_t, SieveResult> &partial_sieve_results,
        std::mutex &res_mutex) :
        N_(N), a_target_(a_target), 
        base_size_(base_size), sieve_radius_(sieve_radius), partial_prime_bound_(partial_prime_bound),
        num_critical_(num_critical), //num_noncritical_(num_noncritical),
        critical_fb_lower_(critical_fb_lower), critical_fb_upper_(critical_fb_upper),
        factor_base_(factor_base), fb_nsqrt_(fb_nsqrt), fb_logp_(fb_logp), fb_magic_num_p_(fb_magic_num_p),
        total_sieved_(total_sieved), 
        sieve_results_(sieve_results), partial_sieve_results_(partial_sieve_results),
        res_mutex_(res_mutex) {
    double a_d = a_target.get_d(),
           N_d = N.get_d();
    double N_over_a_d = N_d / a_d;

    this->sieve_height_ = std::vector<LogType>(2 * sieve_radius + 1);
    this->big_div_ = std::vector<bool>(2 * sieve_radius + 1);

    // OLD:
    // Populate sieve_height_target_ with log|a (x - sieve_radius) ^ 2 - N / a| - log(partial_prime_bound),
    // which is the approximate log of the value of the polynomial.

    // NEW:
    // sieve_height_target_ is no longer an array, but now just a single value. For the vast majority of x, 
    // it is around the value we get from x = 0, i.e. log(N/a) - log(partial_prime_bound)
    double log_partial_prime_bound = std::log(double(partial_prime_bound));
    /*for (uint32_t x = 0; x <= 2 * sieve_radius; x++) {
        double logpoly = std::log(std::abs(
                    a_d * (double(x) - sieve_radius) * (double(x) - sieve_radius) - N_over_a_d
                ));
        if (logpoly > log_partial_prime_bound + TINY_PRIME_PADDING) {
            this->sieve_height_target_[x] = LogType(
                    (logpoly - log_partial_prime_bound - TINY_PRIME_PADDING)
                    * LOG_PRECISION);
        }
    }*/
    this->sieve_height_target_neg_ = LogType(
            -(std::log(N_over_a_d) - log_partial_prime_bound - (TINY_PRIME_PADDING + std::log2(this->base_size_)))
            * LOG_PRECISION);
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
    this->noncritical_fb_idxs_.reserve(this->base_size_ - this->num_critical_);
    // get noncritical indexes by testing each for membership in critical_fb_idxs
    // this is slow but num_critical ~ 10, so it doesn't really matter
    for (uint32_t fb_idx = 0; fb_idx < this->base_size_; fb_idx++) {
        if (this->factor_base_[fb_idx] > TINY_PRIME_BOUND &&
                this->critical_fb_idxs_.find(fb_idx) == this->critical_fb_idxs_.end()) {
            this->noncritical_fb_idxs_.push_back(fb_idx);
        }
    }
    this->num_noncritical_ = this->noncritical_fb_idxs_.size();
}

void Siever::InitCriticalDeltas() {
    uint32_t crt_idx = 0;
    this->crt_indicator_ = std::vector<mpz_class>(this->num_critical_, 0);
    for (const uint32_t &critical_fb_idx : this->critical_fb_idxs_) {
        const PrimeSize p = this->factor_base_[critical_fb_idx];
        mpz_class product_others = this->a_ / p;

        // The idea is that product_others * gamma will be a square root 
        // of N mod p, but 0 mod all other critical primes.
        int32_t gamma = this->fb_nsqrt_[critical_fb_idx] 
                * util::modular_inv_mod_prime<int32_t>(mpz_class(product_others % p).get_si(), p) % p;
        if (gamma > p / 2) gamma = p - gamma;

        mpz_class indicator = product_others * gamma;
        this->crt_indicator_[crt_idx++] = indicator;
    }

    this->a_inv_ = std::vector<PrimeSize>(this->num_noncritical_);
    this->soln_delta_ = std::vector<std::vector<PrimeSize>>(this->num_noncritical_);
    this->soln_delta_neg_ = std::vector<std::vector<PrimeSize>>(this->num_noncritical_);

    // idx indexes over noncritical primes
    for (uint32_t idx = 0; idx < this->num_noncritical_; idx++) {
        const uint32_t noncritical_fb_idx = this->noncritical_fb_idxs_[idx];
        const PrimeSize q = this->factor_base_[noncritical_fb_idx];

        this->a_inv_[idx] = 
            util::modular_inv_mod_prime<int32_t>(mpz_class(this->a_ % q).get_si(), q);
        
        this->soln_delta_[idx].resize(this->num_critical_);
        this->soln_delta_neg_[idx].resize(this->num_critical_);
        for (crt_idx = 0; crt_idx < this->num_critical_; crt_idx++) {
            uint64_t indicator_64 = mpz_class(this->crt_indicator_[crt_idx] % q).get_ui();
            this->soln_delta_[idx][crt_idx] = (indicator_64 << 1) * this->a_inv_[idx] % q;
            this->soln_delta_neg_[idx][crt_idx] = q - this->soln_delta_[idx][crt_idx];
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
        this->soln_ = std::vector<std::pair<PrimeSize, PrimeSize>>(this->num_noncritical_);
        for (uint32_t idx = 0; idx < this->num_noncritical_; idx++) {
            const uint32_t noncritical_fb_idx = this->noncritical_fb_idxs_[idx];
            const PrimeSize q = this->factor_base_[noncritical_fb_idx];
            const int64_t rt = this->fb_nsqrt_[noncritical_fb_idx];
            //const LogType &logq = this->fb_logp_[noncritical_fb_idx];
            
            int64_t b_mod_q = static_cast<int64_t>(mpz_class(this->b_ % q).get_si());
            this->soln_[idx].first  = (((-rt - b_mod_q) * this->a_inv_[idx] + this->sieve_radius_) % q + q) % q;
                /*mpz_class(
                    (((-rt - this->b_ % q) * this->a_inv_[idx] + this->sieve_radius_) % q + q) % q
                ).get_ui();*/
            this->soln_[idx].second = ((( rt - b_mod_q) * this->a_inv_[idx] + this->sieve_radius_) % q + q) % q;
        }

        this->timer_.init_first_poly_time += UpdateTime();
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

    if (gray_code_indicator % 2 == 0) {
        for (uint32_t idx = 0; idx < this->num_noncritical_; idx++) {
            const PrimeSize q = this->factor_base_[this->noncritical_fb_idxs_[idx]];
            this->soln_[idx].first  = 
                (this->soln_[idx].first  + this->soln_delta_[idx][nu]) % q;
            this->soln_[idx].second = 
                (this->soln_[idx].second + this->soln_delta_[idx][nu]) % q;
        }
    } else {
        for (uint32_t idx = 0; idx < this->num_noncritical_; idx++) {
            const PrimeSize q = this->factor_base_[this->noncritical_fb_idxs_[idx]];
            this->soln_[idx].first  = 
                (this->soln_[idx].first  + this->soln_delta_neg_[idx][nu]) % q;
            this->soln_[idx].second = 
                (this->soln_[idx].second + this->soln_delta_neg_[idx][nu]) % q;
        }
    }
    this->timer_.init_next_poly_time += UpdateTime();
    this->poly_++;
}

void Siever::SetHeights() {
    //std::fill(this->sieve_height_.begin(), this->sieve_height_.end(), 0);
    std::fill(this->sieve_height_.begin(), this->sieve_height_.end(), 
            this->sieve_height_target_neg_);
    std::fill(this->big_div_.begin(), this->big_div_.end(), false);
    int32_t sieve_diameter = 2 * this->sieve_radius_;

    // Sieve small primes 
    uint32_t idx = 0;
    for (; idx < this->num_noncritical_; idx++) {
        const uint32_t noncritical_fb_idx = this->noncritical_fb_idxs_[idx];
        const PrimeSize q = this->factor_base_[noncritical_fb_idx];
        const LogType logq = this->fb_logp_[noncritical_fb_idx];

        if (q > sieve_diameter) break;

        // We assume that q = 2 is taken care of by the tiny prime bound
        //const std::pair<PrimeSize, PrimeSize> &solns = this->soln_[idx];
        int32_t x = this->soln_[idx].first, y = this->soln_[idx].second;
        /*const int len = sieve_diameter / q;
        LogType *s1 = this->sieve_height_.data() + this->soln_[idx].first,
                *s2 = this->sieve_height_.data() + this->soln_[idx].second;

        for (int i = 0; i < len; i++, s1 += q, s2 += q) {
            *s1 += logq;
            *s2 += logq;
        }

        if (this->soln_[idx].first + len * q <= sieve_diameter) *s1 += logq;
        if (this->soln_[idx].second + len * q <= sieve_diameter) *s2 += logq;*/

        if (x > y) std::swap(x, y);
        for (; y <= sieve_diameter; x += q, y += q) {
            this->sieve_height_[x] += logq;
            this->sieve_height_[y] += logq;
        }

        if (x <= sieve_diameter) this->sieve_height_[x] += logq;
        /*for (uint32_t x = solns.first; x <= sieve_diameter; x += q) {
            this->sieve_height_[x] += logq;
        }
        for (uint32_t x = solns.second; x <= sieve_diameter; x += q) {
            this->sieve_height_[x] += logq;
        }*/
    }

    // Sieve big primes
    for (; idx < this->num_noncritical_; idx++) {
        const LogType logq = this->fb_logp_[this->noncritical_fb_idxs_[idx]];

        int32_t x = this->soln_[idx].first, y = this->soln_[idx].second;
        if (x <= sieve_diameter) {
            this->sieve_height_[x] += logq;
            this->big_div_[x] = true;
        }
        if (y <= sieve_diameter) {
            this->sieve_height_[y] += logq;
            this->big_div_[y] = true;
        }
    }
}

void InsertPartial(std::vector<SieveResult> &sieve_results, 
        std::unordered_map<uint32_t, SieveResult> &partial_sieve_results,
        const uint32_t partial, const bool sgn,
        const mpz_class &rt, const std::vector<uint32_t> &prime_fb_idxs) {

    std::unordered_map<uint32_t, SieveResult>::const_iterator partial_res_pr = 
        partial_sieve_results.find(partial);

    // If there is already a partial, then merge the partials to create a 
    // full result and add it 
    if (partial_res_pr != partial_sieve_results.end()) {
        const SieveResult& res = partial_res_pr->second;

        // Merge the two prime lists
        std::vector<uint32_t> combined_fb_idx(res.prime_fb_idxs.begin(), res.prime_fb_idxs.end());
        combined_fb_idx.reserve(prime_fb_idxs.size() + res.prime_fb_idxs.size());

        // No point moving primitives
        combined_fb_idx.insert(combined_fb_idx.end(), prime_fb_idxs.begin(), prime_fb_idxs.end());

        // Combine the two partial results into one of the form 
        // (rt1 * rt2) ** 2 = partial ** 2 * (something smooth) mod N
        //
        // TODO: check whether it is good to mod this product by N
        sieve_results.emplace_back(res.numer * rt, mpz_class(partial), 
                sgn ^ res.sgn, std::move(combined_fb_idx));

        // The other partial result (the one stored in the map) can be kept
        // because e1 + e2, e1 + e3, ..., e1 + en forms a basis for 
        // (e1 + e2 + ... + en) perp over GF(2)n 
    } else {
        partial_sieve_results.emplace(partial, 
                SieveResult(std::move(rt), 1, sgn, std::move(prime_fb_idxs)));
    }
}

void Siever::CheckPossibleCache() {
    std::vector<bool> stop(this->cache_size_, false);
    std::vector<uint32_t> big_div_to_check;
    big_div_to_check.reserve(this->cache_size_);

    for (int i = 0; i < this->cache_size_; i++) {
        PossibleResult &pos_res = this->possible_res_cache_[i];

        if (this->big_div_[pos_res.x_shift])
        big_div_to_check.push_back(i);
    }

    //std::vector<int> big_div(this->cache_size_, 0);
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
    /*// Use raw c impl to be faster
    mpz_t n, q, _r, one;
    mpz_init(n);
    mpz_init(q);
    mpz_init(_r);
    mpz_init(one);
    uint32_t true_res; // res as a uint32_t, which is what we care about

    mpz_set(n, polyval.get_mpz_t());
    mpz_set_ui(one, 1);*/

    // Four separate checks for tiny primes, small noncritical primes, big noncritical primes, 
    // and critical primes
    for (uint32_t fb_idx = 0; fb_idx < this->base_size_; fb_idx++) {
        const PrimeSize p = this->factor_base_[fb_idx];
        if (p > TINY_PRIME_BOUND) break;
        
        for (int i = 0; i < this->cache_size_; i++) {
            if (stop[i]) continue;

            PossibleResult &pos_res = this->possible_res_cache_[i];
            bool divisible = (pos_res.polyval % p == 0);
            
            while (divisible) {
                pos_res.polyval /= p;
                pos_res.prime_fb_idxs.push_back(fb_idx);
                if (pos_res.polyval == 1) {
                    stop[i] = true;
                    break;
                }
                divisible = (pos_res.polyval % p == 0);
            }
        }
    }

    uint32_t idx = 0;
    for (; idx < this->num_noncritical_; idx++) {
        const uint32_t noncritical_fb_idx = this->noncritical_fb_idxs_[idx];
        const PrimeSize p = this->factor_base_[noncritical_fb_idx];
        const std::pair<uint32_t, uint32_t> magic_num = this->fb_magic_num_p_[noncritical_fb_idx];

        if (p > 2 * this->sieve_radius_) break;

        for (int i = 0; i < this->cache_size_; i++) {
            if (stop[i]) continue;

            PossibleResult &pos_res = this->possible_res_cache_[i];
            /*bool divisible = ((x - this->soln_[idx].first ) % p == 0) ||
                             ((x - this->soln_[idx].second) % p == 0);*/
            bool divisible = (static_cast<uint32_t>(std::abs(pos_res.x_shift - this->soln_[idx].first))
                                * magic_num.first <= magic_num.second) ||
                             (static_cast<uint32_t>(std::abs(pos_res.x_shift - this->soln_[idx].second))
                                * magic_num.first <= magic_num.second);
            
            while (divisible) {
                pos_res.polyval /= p;
                pos_res.prime_fb_idxs.push_back(noncritical_fb_idx);
                if (pos_res.polyval == 1) {
                    stop[i] = true;
                    break;
                }
                divisible = (pos_res.polyval % p == 0);
            }
            /*// Outer check for divisibility
            if (divisible) {
                true_res = 0;
                mpz_fdiv_q_ui(q, n, p);
                while (true_res == 0) {
                    mpz_swap(n, q); 
                    prime_fb_idxs.push_back(fb_idx);
                    true_res = mpz_fdiv_qr_ui(q, _r, n, p);
                }
            }
            if (mpz_cmp(n, one) == 0) goto done_checking;*/
        }
    }

    // big primes
    for (; idx < this->num_noncritical_; idx++) {
        const uint32_t noncritical_fb_idx = this->noncritical_fb_idxs_[idx];
        const PrimeSize p = this->factor_base_[noncritical_fb_idx];
        const std::pair<uint32_t, uint32_t> magic_num = this->fb_magic_num_p_[noncritical_fb_idx];

        for (const uint32_t i : big_div_to_check) {
            if (stop[i]) continue;
            PossibleResult &pos_res = this->possible_res_cache_[i];

            bool divisible = (static_cast<uint32_t>(std::abs(pos_res.x_shift - this->soln_[idx].first))
                                * magic_num.first <= magic_num.second) ||
                             (static_cast<uint32_t>(std::abs(pos_res.x_shift - this->soln_[idx].second))
                                * magic_num.first <= magic_num.second);
            
            while (divisible) {
                pos_res.polyval /= p;
                pos_res.prime_fb_idxs.push_back(noncritical_fb_idx);
                if (pos_res.polyval == 1) {
                    stop[i] = true;
                    break;
                }
                divisible = (pos_res.polyval % p == 0);
            }
        }
    }

    /*for (int i = 0; i < this->cache_size_; i++) {
        std::cout << big_div[i] << " ";
    }
    std::cout << std::endl;*/

    for (const uint32_t &fb_idx : this->critical_fb_idxs_) {
        const PrimeSize p = this->factor_base_[fb_idx];

        for (int i = 0; i < this->cache_size_; i++) {
            if (stop[i]) continue;

            PossibleResult &pos_res = this->possible_res_cache_[i];
            bool divisible = (pos_res.polyval % p == 0);
            
            while (divisible) {
                pos_res.polyval /= p;
                pos_res.prime_fb_idxs.push_back(fb_idx);
                if (pos_res.polyval == 1) {
                    stop[i] = true;
                    break;
                }
                divisible = (pos_res.polyval % p == 0);
            }
        }
        /*true_res = mpz_fdiv_qr_ui(q, _r, n, p);
        while (true_res == 0) {
            mpz_swap(n, q);
            prime_fb_idxs.push_back(fb_idx);
            true_res = mpz_fdiv_qr_ui(q, _r, n, p);
        }
        if (mpz_cmp(n, one) == 0) goto done_checking;*/
    }
/*done_checking:
    polyval = mpz_get_ui(n);
    mpz_clear(n);
    mpz_clear(q);
    mpz_clear(_r);
    mpz_clear(one);*/
    this->timer_.trial_divide_time += UpdateTime();

    for (int i = 0; i < this->cache_size_; i++) {
        PossibleResult &pos_res = this->possible_res_cache_[i];
        if (pos_res.polyval == 1) {
            pos_res.prime_fb_idxs.insert(pos_res.prime_fb_idxs.end(), 
                    this->critical_fb_idxs_.begin(), this->critical_fb_idxs_.end());
            this->temp_sieve_results_.emplace_back(std::move(pos_res.rt), 1, 
                    pos_res.sgn, std::move(pos_res.prime_fb_idxs));
        } else if (pos_res.polyval.fits_uint_p()
                && pos_res.polyval.get_ui() <= this->partial_prime_bound_) {
            // We don't check this, but all partials are going to be primes (or their negations),
            // simply because our partial_prime_bound_ is not that much (only like 100x) 
            // bigger than the factor base primes, so anything that small will either be smooth 
            // or prime
            uint32_t partial = pos_res.polyval.get_ui();
            pos_res.prime_fb_idxs.insert(pos_res.prime_fb_idxs.end(), 
                    this->critical_fb_idxs_.begin(), this->critical_fb_idxs_.end());
            
            InsertPartial(this->temp_sieve_results_, this->temp_partial_sieve_results_, 
                    partial, pos_res.sgn, pos_res.rt, pos_res.prime_fb_idxs);
        } 
        /*else {
            std::cout << "rejected abs(" << polyval << ") > " << this->partial_prime_bound_ << std::endl;
        }*/
    }
    this->timer_.insert_time += UpdateTime();
}


void Siever::CheckHeights() {
    int32_t sieve_diameter = 2 * this->sieve_radius_;

    // Check results
    // We have (ax+b) ** 2 - N = a(a x ** 2 + 2b x + c) for the integer c = (b ** 2 - N) / a
    // We also calculate rt = ax + b for extracting the final square equality mod N
    mpz_class c = (this->b_ * this->b_ - this->N_) / this->a_;

    for (int32_t x = 0; x <= sieve_diameter;) {
        this->cache_size_ = 0;
        while (x <= sieve_diameter) {
            //if (this->sieve_height_[x] > this->sieve_height_target_[x]) {
            if (this->sieve_height_[x] > 0) {
                this->timer_.wait_res_time += UpdateTime();
                int32_t x_int = x - this->sieve_radius_;

                // We need to set this inline because mpz_class has 
                // no move constructor >:(
                PossibleResult &pos_res = this->possible_res_cache_[this->cache_size_++];
                pos_res.x_shift = x;
                pos_res.rt = this->a_ * x_int + this->b_;
                pos_res.polyval = pos_res.rt + this->b_;
                pos_res.polyval *= x_int;
                pos_res.polyval += c;
                /*pos_res.polyval = this->a_;
                pos_res.polyval *= x_int;
                pos_res.polyval += this->b_ << 1;
                pos_res.polyval *= x_int;
                pos_res.polyval += c; */
                //pos_res.sgn = (sgn(pos_res.polyval) == -1);
                pos_res.sgn = pos_res.polyval < 0;

                pos_res.polyval = abs(pos_res.polyval);

                //mpz_class rt = this->a_ * x_int + this->b_;
                //mpz_class polyval = this->a_ * (x_int * x_int) + this->b_ * (x_int << 1) + c;
                //bool sgn = polyval < 0;

                // list of fb indices for primes dividing (ax+b)^2 - N, so by default, 
                // it should contain all of the critical primes 
                //std::vector<uint32_t> prime_fb_idxs;
                //prime_fb_idxs.reserve(20);

                pos_res.prime_fb_idxs = std::vector<uint32_t>();
                //this->possible_res_cache_[cache_size].prime_fb_idxs.reserve(20);

                //this->possible_res_cache_.emplace_back(x, std::move(rt), std::move(abs(polyval)), 
                        //sgn, std::move(prime_fb_idxs));

                this->timer_.check_time += UpdateTime();
                if (this->cache_size_ == POSSIBLE_RES_CACHE_SIZE) {
                    x++;
                    break;
                }
            }
            x++;
        }

        CheckPossibleCache();
    }
}

void Siever::SievePoly() {
    SetHeights();
    this->timer_.set_height_time += UpdateTime();
    CheckHeights();
}

std::chrono::system_clock::time_point time() {
    return std::chrono::system_clock::now();
}

double Siever::UpdateTime() {
    std::chrono::system_clock::time_point time_now = time();
    double res = 
        std::chrono::duration_cast<std::chrono::nanoseconds>(time_now - this->time_prev_).count() 
        / 1'000'000'000.0;
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
        SievePoly();

        this->total_sieved_++;
    }
}

void Siever::FlushResults() {
    UpdateTime();
    std::move(this->temp_sieve_results_.begin(), this->temp_sieve_results_.end(),
            std::back_inserter(this->sieve_results_));

    for (const std::pair<const uint32_t, SieveResult> &partial_p : 
            this->temp_partial_sieve_results_) {
        InsertPartial(this->sieve_results_, this->partial_sieve_results_, 
                partial_p.first, 
                partial_p.second.sgn, partial_p.second.numer, partial_p.second.prime_fb_idxs);
    }

    this->temp_sieve_results_.clear();
    this->temp_partial_sieve_results_.clear();

    this->timer_.flush_time += UpdateTime();
}

Timer Siever::get_timer_() const { return this->timer_; }
void Siever::reset_timer_() { this->timer_ = Timer(); }
