#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <gmpxx.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

#include <gmpxx.h>
extern "C" {
#include "dgs/dgs/dgs.h"
}
#include "encryption/labeled.h"
#include "hash/chameleon.h"
#include "params.h"
#include "sphf/sphf.h"
#include "trapdoor/mp12.h"

NTL_CLIENT
using namespace std;

void init(ZZ modulus) {
  ZZ_p::init(modulus);
  SetNumThreads(4);
}

int main(void) {
  // const uint32_t n = DIMENSION;
  // const uint32_t k = Q_BITS;
  // const uint32_t q = MODULUS;
  // const double_t rounding = RANDOMIZED_ROUNDING_PARAMETER;
  // const double_t sigma = GAUSS_WIDTH;
  const uint32_t n = 128;
  const uint32_t k = 12;
  const uint32_t q = 2053;
  const double_t rounding = RANDOMIZED_ROUNDING_PARAMETER;
  const double_t sigma = 2 * sqrt(n);  // set sigma = 2*sqrt(n) as in BBDQ17
                                       // to ensure that the scheme is
                                       // correct

  init(ZZ(q));

  Trapdoor::MP12* trapdoor = new Trapdoor::MP12(n, k, q, sigma, rounding);

  // Encryption
  vec_ZZ_p f_coeffs;
  f_coeffs.SetLength(n);
  for (uint32_t i = 0; i < n; i++) {
    f_coeffs[i] = F_COEFFS_128_2053[i];
  }
  Encryption::Labeled enc{n, k, q, sigma, rounding, f_coeffs, trapdoor};
  mat_ZZ_p A0;
  mat_ZZ_p A1;
  mat_ZZ R;
  enc.keygen(A0, A1, R);
  mat_ZZ_p Au;
  vec_ZZ_p u;
  vec_ZZ_p c;

  u = enc.taggen(Au, A0, A1);

  double_t t = GetWallTime();
  enc.tag_encrypt(c, Au, A0, A1, 1);
  cout << GetWallTime() - t << endl;

  uint8_t mu;
  enc.decrypt(mu, u, c, R, A0, A1);

  cout << (int)mu << endl;

  delete trapdoor;

  return 0;
}