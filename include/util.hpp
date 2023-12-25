// util.hpp
#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <gmpxx.h>
#include <functional>
#include <map>

namespace util {

// 2-adic valuation of n.
// Returns a pair containing \\nu_2(n) and n / 2 ** \\nu_2(n)
std::pair<uint32_t, mpz_class> two_adic_val(mpz_class n);

// If n = a ** b where a is minimal, returns {a, b}.
// If n is not a perfect power, then b = 1, so this returns {n, 1}.
std::pair<mpz_class, uint32_t> factor_perfect_power(const mpz_class &n);

// Power of two order, i.e. the least n such that t ** (2 ** n) is 1 mod p
uint32_t pow2_ord(const mpz_class &t, const mpz_class &p);

// Prints the prime factorization of n (in p1 ^ e1 * p2 * ..., where e.g. e2 = 1, form)
// given a function to find nontrivial factors of n.
void print_prime_fact(const mpz_class &n,
              const std::function<mpz_class(const mpz_class &)> &&find_nontrivial_factor);

} // namespace util

#endif // UTIL_HPP_
