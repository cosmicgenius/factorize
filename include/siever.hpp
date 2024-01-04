// siever.hpp
#ifndef SIEVER_HPP_
#define SIEVER_HPP_

#include <cstdint>
#include <gmpxx.h>
#include <unordered_map>
#include <set>
#include <vector>

typedef uint32_t PrimeSize;

// Fixed point type for storing logs that stores 6
// bits behind the point in a 16 bit integer.
const uint16_t LOG_PRECISION = 1 << 6;
typedef uint16_t LogType;

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
};

class Siever {
private:
    const mpz_class N_;
    const mpz_class a_target_;

    const uint32_t base_size_;
    const uint32_t sieve_radius_;
    const uint32_t partial_prime_bound_;
    
    const uint32_t num_critical_;
    const uint32_t num_noncritical_;
    const uint32_t critical_fb_lower_;
    const uint32_t critical_fb_upper_;

    const std::vector<PrimeSize> &factor_base_;
    const std::vector<PrimeSize> &fb_nsqrt_;
    const std::vector<LogType> &fb_logp_;

    uint32_t &total_sieved_;

    std::vector<SieveResult> &sieve_results_;
    std::unordered_map<uint32_t, SieveResult> &partial_sieve_results_;

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
    std::vector<std::pair<PrimeSize, PrimeSize>> soln_;

    // The value sieve_height[x] will be populated with 
    // the log of the (squarefree) fb-smooth part of 1/a * ((a(x - sieve_radius) + b)^2 - N)
    std::vector<LogType> sieve_height_;

    // Heights to reach to trigger trial division check for smoothness
    // Value is approximately log|1/a * ((ax)^2 - N)| - log(partial_prime_bound) 
    // ~ log|a * x^2 - N / a| - log(partial_prime_bound) 
    std::vector<LogType> sieve_height_target_;

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
    void CheckHeights();
    void CheckSmoothness(const int32_t x, mpz_class &polyval, std::vector<uint32_t> &prime_fb_idxs);
    void InsertPartial(const uint32_t partial, const bool sgn,
            const mpz_class &rt, const std::vector<uint32_t> &prime_fb_idxs);

public:
    explicit Siever(const mpz_class &N, const mpz_class &a_target, 
            const uint32_t &base_size, const uint32_t &sieve_radius, const uint32_t &large_prime_bound, 
            const uint32_t &num_critical, const uint32_t &num_noncritical, 
            const uint32_t &critical_fb_lower, const uint32_t &critical_fb_upper, 
            const std::vector<PrimeSize> &factor_base, 
            const std::vector<PrimeSize> &fb_nsqrt, const std::vector<LogType> &fb_logp,
            uint32_t &total_sieved,
            std::vector<SieveResult> &sieve_results, std::unordered_map<uint32_t, SieveResult> &partial_sieve_results_);

    // Sieves an entire poly group (chosen randomly) for smooth values
    void SievePolynomialGroup();
};

#endif // SIEVER_HPP_
