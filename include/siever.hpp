// siever.hpp
#ifndef SIEVER_HPP_
#define SIEVER_HPP_

#include <array>
#include <chrono>
#include <cstdint>
#include <gmpxx.h>
#include <mutex>
#include <set>
#include <unordered_map>
#include <vector>

typedef int32_t PrimeSize;

// Fixed point type for storing logs that stores 6
// bits behind the point in a 16 bit integer.
typedef int32_t LogType;
const LogType LOG_PRECISION = 1 << 12;

struct Timer {
    // We keep track of these times for prof
    double init_grp_time = 0;
    double init_first_poly_time = 0;
    double init_next_poly_time = 0;
    double set_height_time = 0;
    double wait_res_time = 0;
    double check_time = 0;
    double trial_divide_time = 0;
    double insert_time = 0;
    double flush_time = 0;
    double kernel_time = 0;

    Timer& operator+=(const Timer &rhs);
    Timer operator+(const Timer &rhs) const;
};

// Result from sieving, where numer and denom are values such that 
// such that (numer / denom) ** 2 = polyval mod N, and sgn and prime_fb_idxs
// denotes the fb indices of either a full factorization of polyval 
// or a partial factorization thereof
//
// It is not guaranteed that prime_fb_idxs is sorted
//
// We will often see denom = 1 and numer = rt, for polyval = rt ** 2 - N,
// but this is not guaranteed.
struct SieveResult {
    mpz_class numer;
    mpz_class denom;
    bool sgn; // true is negative (i.e. presence of the ``-1'' prime)
    std::vector<uint32_t> prime_fb_idxs;
    
    SieveResult(const mpz_class &&numer, const mpz_class &&denom, 
            const bool sgn, const std::vector<uint32_t> &&prime_fb_idxs);

    // Moves
    SieveResult& operator=(SieveResult &&rhs);
    SieveResult(SieveResult &&rhs) noexcept;

    // Notably, a copy constructor is NOT provided 
    // (there should never be a reason to copy a result instead 
    // of simply moving it)
    SieveResult& operator=(const SieveResult &rhs) = delete;
    SieveResult(const SieveResult &rhs) = delete;
};

constexpr int POSSIBLE_RES_CACHE_SIZE = 8;

// Possible result tied to a single poly (a, b)
struct PossibleResult {
    int32_t x_shift; // The value of x + sieve_radius
    mpz_class rt; // The value of ax+b
    mpz_class polyval; // The value of ((ax+b)^2 - N) / a that is being factored, 
                       // polyval can be destroyed
    bool sgn; // sgn of (ax+b)^2 - N
    std::vector<uint32_t> prime_fb_idxs; // populated by prime factors of polyval
    
    PossibleResult(const int32_t x_shift, const mpz_class &&rt, const mpz_class &&polyval, 
            const bool sgn, const std::vector<uint32_t> &&prime_fb_idxs);

    PossibleResult() = default;

    // Moves 
    PossibleResult& operator=(PossibleResult &&rhs);
    PossibleResult(PossibleResult &&rhs) noexcept;

    // See above
    PossibleResult& operator=(const PossibleResult &rhs) = delete;
    PossibleResult(const PossibleResult &rhs) = delete;
};

class Siever {
private:
    const mpz_class N_;
    const mpz_class a_target_;

    const uint32_t base_size_;
    const uint32_t sieve_radius_;
    const uint32_t partial_prime_bound_;
    
    const uint32_t num_critical_;
    uint32_t num_noncritical_;
    const uint32_t critical_fb_lower_;
    const uint32_t critical_fb_upper_;

    const std::vector<PrimeSize> factor_base_;
    const std::vector<PrimeSize> fb_nsqrt_;
    const std::vector<LogType> fb_logp_;
    const std::vector<std::pair<uint32_t, uint32_t>> fb_magic_num_p_;

    uint32_t &total_sieved_;

    std::vector<SieveResult> &sieve_results_;
    std::unordered_map<uint32_t, SieveResult> &partial_sieve_results_;

    // Mutex for safe accessing of the above result objects
    std::mutex &res_mutex_;

    std::vector<SieveResult> temp_sieve_results_;
    std::unordered_map<uint32_t, SieveResult> temp_partial_sieve_results_;

    uint32_t poly_ = 0;
    mpz_class a_, b_;

    // fb_idxs are ONLY for indexing into the factor base 
    // and containers marked ``fb''
    std::set<uint32_t> critical_fb_idxs_;

    // We want b^2 = N mod a, which we created by summing 
    // \pm crt_indicators that are 0 modulo all critical primes except for one prime p,
    // for which it is a square root of N
    //
    // These correspond to B_l in https://math.dartmouth.edu/~carlp/implementing.pdf
    std::vector<mpz_class> crt_indicator_;

    std::vector<uint32_t> noncritical_fb_idxs_;
    // For each noncritical prime q, store a ** -1, every possible 
    // value of 2 * crt_indicator * a ** -1 (to change solutions), and the two solutions 
    // to (ax + b)^2 = N mod q
    std::vector<PrimeSize> a_inv_;
    std::vector<std::vector<PrimeSize>> soln_delta_;
    std::vector<std::vector<PrimeSize>> soln_delta_neg_;
    std::vector<std::pair<PrimeSize, PrimeSize>> soln_;

    // The value sieve_height[x] will be populated with 
    // the log of the (squarefree) fb-smooth part of 1/a * ((a(x - sieve_radius) + b)^2 - N)
    std::vector<LogType> sieve_height_;
    std::vector<bool> big_div_; // whether a certain polyval is divisible by a "big" (i.e. > sieve_diameter)
                                // noncritical prime

    // Heights to reach to trigger trial division check for smoothness
    // Value is approximately log|1/a * ((ax)^2 - N)| - log(partial_prime_bound) 
    // ~ log|a * x^2 - N / a| - log(partial_prime_bound) 
    //std::vector<LogType> sieve_height_target_;
    LogType sieve_height_target_neg_;

    // (Re)chooses critical primes that multiply to a, shared by poly group
    void InitCriticalPrimes();

    // Precomputes deltas for faster transitions between polynomials
    // in this same poly group
    void InitCriticalDeltas();

    // Setups up data (e.g. b, soln) for the next polynomial (current polynomial tracked by 
    // this->poly_). Undefined behaviour if this->poly_ >= 2 ** this->num_critical_
    void InitNextPoly();

    // Sieves polynomial for smooth values 
    void SievePoly();

    void SetHeights();

    // Checks heights and throws results/partials in temp
    void CheckHeights();
    
    int cache_size_;
    std::array<PossibleResult, POSSIBLE_RES_CACHE_SIZE> possible_res_cache_;
    void CheckPossibleCache();

    Timer timer_;
    std::chrono::system_clock::time_point time_prev_;

    // Return time from now to time_prev_ and updates time_prev_
    double UpdateTime();

public:
    Siever(const mpz_class N, const mpz_class a_target, 
            const uint32_t base_size, const uint32_t sieve_radius, const uint32_t large_prime_bound, 
            const uint32_t num_critical, //const uint32_t num_noncritical, 
            const uint32_t critical_fb_lower, const uint32_t critical_fb_upper, 
            const std::vector<PrimeSize> factor_base, 
            const std::vector<PrimeSize> fb_nsqrt, const std::vector<LogType> fb_logp,
            const std::vector<std::pair<uint32_t, uint32_t>> fb_magic_num_p,
            uint32_t &total_sieved,
            std::vector<SieveResult> &sieve_results, std::unordered_map<uint32_t, SieveResult> &partial_sieve_results_,
            std::mutex &res_mutex);

    Siever() = delete;

    // No copying, default moving
    Siever(const Siever &rhs) = delete;
    Siever(Siever &&rhs) = default;

    // Sieves an entire poly group (chosen randomly) for smooth values
    void SievePolynomialGroup();

    // Flush temp and put the results into the main results
    // Requires lock
    void FlushResults();

    Timer get_timer_() const;
    void reset_timer_();
};

#endif // SIEVER_HPP_
