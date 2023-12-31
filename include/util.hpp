// util.hpp
#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <cstdint>
#include <gmpxx.h>
#include <functional>
#include <map>
#include <vector>

namespace util {

// 2-adic valuation of n.
// Returns a pair containing \\nu_2(n) and n / 2 ** \\nu_2(n)
std::pair<uint32_t, mpz_class> two_adic_val(mpz_class n);

// If n = a ** b where a is minimal, returns {a, b}.
// If n is not a perfect power, then b = 1, so this returns {n, 1}.
std::pair<mpz_class, uint32_t> factor_perfect_power(const mpz_class &n);

// Returns a vector of all primes less than n.
std::vector<uint32_t> primes_less_than(const uint32_t &n);

// Returns a ** exp mod m
uint32_t pow(const uint32_t &a, uint32_t exp, const uint32_t &m);

// Returns a ** exp
mpz_class pow(const mpz_class &a, uint32_t exp);

// Finds some x such that x^2 = n (mod p)
//
// Throws std::domain_error if n is not a quadratic residue mod p.
uint32_t square_root_modulo_prime(const mpz_class &n, const uint32_t &p);

// Prints the prime factorization of n (in p1 ^ e1 * p2 * ..., where e.g. e2 = 1, form)
// given a function to find nontrivial factors of n.
void print_prime_fact(const mpz_class &n,
              const std::function<mpz_class(const mpz_class &)> &&find_nontrivial_factor);

// Finds a ** -1 mod p.
// Assumes that p is prime. The use of int32_t instead of uint32_t is
// to allow use of Extended Euclidean algorithm.
int32_t modular_inv_mod_prime(const int32_t &a, const int32_t &p);


} // namespace util

#endif // UTIL_HPP_
