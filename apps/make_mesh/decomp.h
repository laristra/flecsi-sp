#pragma once

#include <array>
#include <vector>
#include <algorithm>

namespace cartmesh {

using size_t = std::size_t;

inline std::vector<size_t> sieve(size_t n)
{
  std::vector<bool> prime(n+1, true);
  std::vector<size_t> ret;

  for (size_t p = 2; p*p <= n; p++) {
    if (prime[p]) {
      for (size_t i = 2*p; i <= n; i += p) {
        prime[i] = false;
      }
    }
  }

  for (size_t p = 2; p <= n; p++) {
    if (prime[p])
      ret.push_back(p);
  }

  return ret;
}


inline std::vector<size_t> factor(size_t np)
{
  std::vector<size_t> facs;
  auto primes = sieve(np);

  size_t p = 0;
  for (size_t i = 0; i < primes.size(); i++) {
    if (np % primes[p] == 0) {
      facs.push_back(primes[p]);
      np = np / primes[p];
    } else {
      p++;
    }
  }

  std::sort(facs.begin(), facs.end());
  std::reverse(facs.begin(), facs.end());

  return facs;
}


template<short D>
std::array<size_t, D> decomp(std::array<size_t, D> n, size_t np)
{
  std::array<size_t,D> parts;
  parts.fill(1);

  auto facs = factor(np);

  // greedy decomposition
  for (auto fac : facs) {
    auto maxind = std::distance(n.begin(),
                                std::max_element(n.begin(), n.end()));
    parts[maxind] *= fac;
    n[maxind] /= fac;
  }

  return parts;
}

}
