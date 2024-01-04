// util.cpp
#include "../include/util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <map>
#include <math.h>
#include <ostream>
#include <stack>
#include <string>
#include <vector>
#include <utility>

std::pair<uint32_t, mpz_class> util::two_adic_val(mpz_class n) {
    for (uint32_t r = 0;; r++) {
        // if last bit is 1, it is odd, and we return
        if (mpz_tstbit(n.get_mpz_t(), 0)) return {r, n};
        n >>= 1;
    }
}

std::pair<mpz_class, uint32_t> util::factor_perfect_power(const mpz_class &n) {
    uint32_t rmax = mpz_sizeinbase(n.get_mpz_t(), 2) - 1;
    mpz_t root;
    mpz_init(root);
    for (uint32_t r = rmax; r >= 2; r--) {
        if (mpz_root(root, n.get_mpz_t(), r)) {
            // root is assgined the correct value by the function mpz_root
            mpz_class ans(root);
            mpz_clear(root);
            return {ans, r};
        }
    }

    mpz_clear(root);
    return {n, 1};
}

// Wheel augmented sieve of Eratosthenes with prime basis 2, 3, 5, 7.
std::vector<uint32_t> util::primes_less_than(const uint32_t &n) {
    if (n <= 2) return {};
    if (n <= 3) return {2};
    if (n <= 5) return {2, 3};
    if (n <= 7) return {2, 3, 5};

    uint32_t wheel_residues[48]; // totient(210) = 48
    std::map<uint32_t, uint32_t> wr_back; // opposite direction of wheel_residues
    int min_res = -1; // smallest residue bigger than n % 210

    int cur_residue = 0;
    for (uint32_t i = 1; i < 210; i++) {
        if (i % 2 != 0 && i % 3 != 0 && i % 5 != 0 && i % 7 != 0) {
            wr_back[i] = cur_residue;

            if (min_res == -1 && i > n % 210) min_res = i;
            wheel_residues[cur_residue++] = i;
        }
    }

    if (min_res == -1) min_res = 1;

    uint32_t rootn = (uint32_t)sqrt(n);
    uint32_t num_residues = 48 * (n / 210) + wr_back[min_res];

    // is_prime[n] represents whether the nth wheel residue
    // is prime, with is_prime[0] corresponding to 1, is_prime[1]
    // corresponding to 11, then 13, etc.
    std::vector<bool> is_prime(num_residues, true);

    is_prime[0] = false;

    auto get_residue = [&wheel_residues] (uint32_t idx) 
        -> uint32_t { return 210 * (idx / 48) + wheel_residues[idx % 48]; };

    auto get_index = [&wr_back] (uint32_t res) 
        -> uint32_t { return 48 * (res / 210) + wr_back[res % 210]; };

    for (uint32_t idx_p = 1; idx_p < num_residues; idx_p++) {
        uint32_t res_p = get_residue(idx_p);
        if (res_p > rootn) break;

        if (is_prime[idx_p]) {
            for (uint32_t idx_q = idx_p; idx_q < num_residues; idx_q++) {
                uint32_t prod_pq = res_p * get_residue(idx_q);
                if (prod_pq > n) break;
                is_prime[get_index(prod_pq)] = false;
            }
        }
    }

    std::vector<uint32_t> primes = {2, 3, 5, 7};
    for (uint32_t idx = 0; idx < num_residues; idx++) {
        if (is_prime[idx]) {
          uint32_t res = get_residue(idx);
          if (res >= n) break;
          primes.push_back(res);
        }
    }
    return primes;
}

// Exponentiation by squaring
uint32_t util::pow(const uint32_t &a, uint32_t exp, const uint32_t &m) {
    // Require 64 bits since we are squaring before doing mod m (< 2 ** 32)
    uint64_t ans = 1, pow = a;

    while (exp) {
        if (exp & 1) ans = (ans * pow) % m;
        exp >>= 1;
        pow = (pow * pow) % m;
    }
    return ans;
}

mpz_class util::pow(const mpz_class &a, uint32_t exp) {
    mpz_class ans = 1, pow = a;

    while (exp) {
        if (exp & 1) ans *= pow;
        exp >>= 1;
        pow *= pow;
    }
    return ans;
}

uint32_t pow2_ord(const uint32_t &t, const uint32_t &p) {
    // Require 64 bits since we are squaring before doing mod p (< 2 ** 32)
    uint64_t t_mod_p = t % p;
    if (t_mod_p == 0) return 0;

    for (uint32_t o = 0;; o++) {
        if (t_mod_p == 1) return o;
        t_mod_p = t_mod_p * t_mod_p % p;
    }
}

// Tonelli Shanks algorithm.
uint32_t util::square_root_modulo_prime(const mpz_class &n, const uint32_t &p) {
    uint64_t n_mod_p = mpz_class(n % p).get_ui();

    // Trivial cases
    if (n_mod_p == 0) return 0;
    if (p == 2) return 1;

    if (mpz_kronecker_ui(n.get_mpz_t(), p) == -1) {
        throw std::domain_error(
        n.get_str() + " is not a quadratic residue modulo " + std::to_string(p));
    }

    // If p = 3 mod 4, then (n^((p+1)/4))^2 = n * n^((p-1)/2) = n mod p
    if (p % 4 == 3) return util::pow(n_mod_p, ((p + 1) / 4), p);

    // Tonelli Shanks implementation 
    std::pair<uint32_t, mpz_class> res = util::two_adic_val(p - 1);
    uint32_t S = res.first;
    uint32_t Q = res.second.get_ui();

    uint32_t z = 2;
    while (z == 0 || mpz_kronecker_ui(mpz_class(z).get_mpz_t(), p) == 1) {
        z = rand() % p; // we don't really care about uniformity here
                        // (we just need a single nqr), so this is fine
    }

    uint32_t M;

    // Require 64 bits since we are squaring before doing mod p (< 2 ** 32)
    uint64_t c, R, t;
    uint64_t L = util::pow(n_mod_p, (Q - 1) >> 1, p);
    M = S;
    c = util::pow(z, Q, p);
    R = L * n_mod_p % p;
    t = L * R % p;

    while (t != 1) {
        //std::cout << M << " " << c << " " << R << " " << t << std::endl;

        uint32_t o = pow2_ord(t, p);
        uint64_t b = util::pow(c, util::pow(2, M - o - 1, p - 1), p);

        M = o;
        c = b * b % p;
        t = t * c % p;
        R = R * b % p;
    }

    return R;
}

void get_prime_factors(std::map<mpz_class, uint32_t> &rop,
                       const mpz_class &n,
                       const std::function<mpz_class(const mpz_class &)> 
                           &&find_nontrivial_factor) {
    rop.clear();

    // stack holding (possible) composites to be factored, with multiplicity
    std::stack<std::pair<mpz_class, uint32_t>> st;

    mpz_class one = 1;

    // checks last bit for parity
    if (mpz_tstbit(n.get_mpz_t(), 0) == 0) {
        std::pair<uint32_t, mpz_class> p = util::two_adic_val(n);
        rop[2] = p.first;

        if (p.second == one) return;
        st.push({p.second, 1});
    } else {
        st.push({n, 1});
    }

    while (!st.empty()) {
        // we are examining top.first ** top.second
        std::pair<mpz_class, uint32_t> top = st.top();
        st.pop();

        // write top.first = p.first ** p.second for p.first minimal
        std::pair<mpz_class, uint32_t> p = util::factor_perfect_power(top.first);

        // get a nontrivial factor of p.first
        mpz_class f = find_nontrivial_factor(p.first);

        // if p.first is prime, we add a contribution of p.first ** (top.second * p.second)
        if (f == one) {
            rop[p.first] += top.second * p.second;
        } else {
            st.push({f, top.second * p.second});
            st.push({top.first / f, top.second * p.second});
        }
    }
}

void util::print_prime_fact(const mpz_class &n,
              const std::function<mpz_class(const mpz_class &)> &&find_nontrivial_factor) {
  if (n == 1) {
    std::cout << "n = 1" << std::endl;
    return;
  }

  std::map<mpz_class, uint32_t> prime_fact;
  //get_prime_factors(prime_fact, n, [](const mpz_class &b) { return find_nontrivial_factor(b); });
  get_prime_factors(prime_fact, n, std::move(find_nontrivial_factor));

  auto power_text = [](std::pair<mpz_class, uint32_t> p) {
    if (p.second == 1) {
      return p.first.get_str();
    }
    return p.first.get_str() + " ^ " + std::to_string(p.second);
  };

  std::cout << "the prime factorization of n is n = ";

  for (auto it = prime_fact.begin(); it != prime_fact.end(); it++) {
    if (it == prime_fact.begin()) {
      std::cout << power_text(*it);
    } else {
      std::cout << " * " << power_text(*it);
    }
  }
  std::cout << std::endl;
}

// Extended Euclidean Algorithm
int32_t util::modular_inv(const int32_t &a, const int32_t &m) {
    int u = a % m, v = m;
    int x1 = 1, x2 = 0;

    while (u != 1) {
      int q = v / u, r = v % u;
      int x = x2 - q * x1;

      v = u;
      u = r;
      x2 = x1;
      x1 = x;
    }
    return (x1 + m) % m;
}

