#include "random.h"

NTL_CLIENT
using namespace std;

ZZ_p random_Zp() {
  const ZZ bnd = ZZ_p::modulus();
  return conv<ZZ_p>(randint(bnd));
}

ZZ randint(const ZZ bnd) {
  if (bnd <= 1) {
    return ZZ(0);
  }
  ZZ z;
  size_t nbchar = NumBits(bnd) / (sizeof(unsigned char) * CHAR_BIT) + 1;
  unsigned char* rnd = new unsigned char[nbchar];
  ZZ mask(2);
  mask <<= NumBits(bnd);
  mask -= 1;

  do {
    fastrandombytes(rnd, nbchar);
    ZZFromBytes(z, rnd, nbchar);
    z &= mask;
  } while (z >= bnd);
  delete[] rnd;
  return z;
}

// generate a uniformly random matrix in Z_p
void generate_random_matrix(Mat<ZZ_p>& A) {
  long n = A.NumRows();
  long m = A.NumCols();
  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      A[i][j] = random_Zp();
}

vec_ZZ rand_vec_bin(uint32_t len) {
  vec_ZZ m;
  m.SetLength(len);
  unsigned char r;
  for (uint32_t i = 0; i < len; i++) {
    fastrandombytes(&r, 1);
    m[i] = (r % 2);
  }
  return m;
}

vec_ZZ_p rand_vec_Zp(uint32_t len) {
  vec_ZZ_p m;
  m.SetLength(len);
  for (uint32_t i = 0; i < len; i++) {
    m[i] = random_Zp();
  }
  return m;
}

uint8_t real_random(double_t prob) {
  std::default_random_engine gen;
  std::uniform_real_distribution<double_t> dis(0.0, 1.0);

  if (dis(gen) < prob) {
    return 1;
  }
  return 0;
}