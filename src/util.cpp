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
#include <ostream>
#include <stack>
#include <vector>
#include <utility>

std::pair<uint32_t, mpz_class> util::two_adic_val(mpz_class n) {
    for (uint32_t r = 0;; r++) {
        // if last bit is 1, it is odd, and we return
        if (mpz_tstbit(n.get_mpz_t(), 0)) return {r, n};
        n /= 2;
    }
}

/*
 * Assuming n > 0 and r \ne 0, sets rop to the 
 * value of floor(n ^ (1 / r))
 *
 * Finds the exact value via binary search
 * after estimating via bounding by two powers of two
 * (there may be a faster way to estimate)
 */
void get_integral_rth_root(mpz_t rop, const mpz_t &n, const uint32_t &r) {
    mpz_t lo, hi, mid, one, tot, pow_val;
    mpz_init(lo);
    mpz_init(hi);
    mpz_init(mid);
    mpz_init_set_ui(one, 1);
    mpz_init(tot);
    mpz_init(pow_val);

    mpz_mul_2exp(lo, one, (mpz_sizeinbase(n, 2) - 1) / r);
    mpz_mul_2exp(hi, lo, 1);
        
    while (mpz_cmp(hi, lo) > 0) {
        mpz_add(tot, lo, hi);
        mpz_tdiv_q_2exp(mid, tot, 1);
        mpz_pow_ui(pow_val, mid, r);

        if (mpz_cmp(pow_val, n) < 0) {
            mpz_add(lo, mid, one);
        } else {
            mpz_set(hi, mid);
        }
    }
    mpz_set(rop, lo);

    mpz_clear(lo);
    mpz_clear(hi);
    mpz_clear(mid);
    mpz_clear(one);
    mpz_clear(tot);
    mpz_clear(pow_val);
}

std::pair<mpz_class, uint32_t> util::factor_perfect_power(const mpz_class &n) {
    uint32_t rmax = mpz_sizeinbase(n.get_mpz_t(), 2) - 1;
    mpz_t n_t, root, pow;
    mpz_init(n_t);
    mpz_init(root);
    mpz_init(pow);

    mpz_set(n_t, n.get_mpz_t());
    for (uint32_t r = rmax; r >= 2; r--) {
        get_integral_rth_root(root, n_t, r);
        mpz_pow_ui(pow, root, r);
        if (mpz_cmp(pow, n_t) == 0) {
            mpz_class ans(root);

            mpz_clear(n_t);
            mpz_clear(root);
            mpz_clear(pow);

            return {ans, r};
        }
    }

    mpz_clear(n_t);
    mpz_clear(root);
    mpz_clear(pow);

    return {n, 1};
}

/*uint32_t pow2_ord(const mpz_class &t, const mpz_class &p) {
  mpz_class t_mod_p = t % p;
  if (t_mod_p == 0) {
    return 0;
  }

  for (uint32_t o = 0;; o++) {
    if (t_mod_p == 1) {
      return o;
    }
    t_mod_p = t_mod_p * t_mod_p % p;
  }
}*/

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
