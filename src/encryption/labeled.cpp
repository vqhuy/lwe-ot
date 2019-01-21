#include "encryption/labeled.h"
#include <NTL/BasicThreadPool.h>
#include <algorithm>
#include <cmath>
#include "prng/random.h"

Encryption::Labeled::Labeled(const uint32_t dimension,
                             const uint32_t modulus_len,
                             const uint32_t modulus,
                             const double_t sigma,
                             const double_t rounding,
                             const vec_ZZ_p f_coeffs,
                             Trapdoor::MP12* mp12) {
  _n = dimension;
  _q = modulus;
  _k = modulus_len;
  _sigma = sigma;
  _r = rounding;
  _m = dimension * (modulus_len + 2);
  _f_coeffs = f_coeffs;
  _trapdoor = mp12;

  // compute other things
  _t = _sigma * sqrt(_m) * _r;
  _B = 2 * _t * sqrt(_m);
  _B_prime = _q / sqrt(_m);

  // compute encoding vector [0, ..., 0, q/2]
  _encoding.SetLength(_m);
  _encoding[_m - 1] = ceil(modulus / 2.0);
}

vec_ZZ Encryption::Labeled::encode(uint8_t mu) {
  return mu * _encoding;
}

vec_ZZ_p Encryption::Labeled::random_element_in_U() {
  vec_ZZ_p u;
  u.SetLength(_n);
  uint32_t deg = conv<uint32_t>(randint(ZZ(_n)));
  ZZ_p coeff;
  while ((coeff = random_Zp()) == 0)
    ;
  u[deg] = coeff;
  return u;
}

mat_ZZ_p Encryption::Labeled::convolution(vec_ZZ_p u) {
  mat_ZZ_p hu;
  hu.SetDims(_n, _n);

  vec_ZZ_p hu_prev;
  hu_prev.SetLength(_n);

  hu[0] = u;  // u * 1
  for (uint32_t i = 1; i < _n; i++) {
    hu_prev = hu[i - 1];
    std::rotate(hu_prev.begin(), hu_prev.end() - 1, hu_prev.end());
    hu_prev[0] = 0;
    ZZ_p last_coeff = hu[i - 1][_n - 1];
    hu[i] = hu_prev - last_coeff * _f_coeffs;
  }

  return hu;
}

void Encryption::Labeled::keygen(mat_ZZ_p& A0, mat_ZZ_p& A1, mat_ZZ& R) {
  // note: A1 = -(A_bar * R)
  _trapdoor->generate_0(R, A0, A1);
}

vec_ZZ_p Encryption::Labeled::taggen(mat_ZZ_p& Au,
                                     const mat_ZZ_p A0,
                                     const mat_ZZ_p A1) {
  vec_ZZ_p u;
  u = this->random_element_in_U();

  // Au = [I | A0 | A1 + h(u)*G]
  // cheating: just compute Au = A1 + h(u)*G
  mul(Au, this->convolution(u), _trapdoor->gadget());
  add(Au, A1, Au);
  return u;
}

void Encryption::Labeled::tag_encrypt(vec_ZZ_p& c,
                                      const mat_ZZ_p Au,
                                      const mat_ZZ_p A0,
                                      const mat_ZZ_p A1,
                                      const uint8_t mu) {
  vec_ZZ_p s;
  vec_ZZ e;
  ZZ e_norm;

  // sample s uniformly
  s = rand_vec_Zp(_n);

  // sample e from discrete Gauss
  e.SetLength(_m);
  do {
    NTL_EXEC_RANGE(_m, first, last)
    for (uint32_t i = first; i < last; i++) {
      e[i] = _trapdoor->sample(_sigma, 0);
      // e[i] = _trapdoor->sample(_t, 0);
    }
    NTL_EXEC_RANGE_END
    InnerProduct(e_norm, e, e);
  } while (SqrRoot(conv<RR>(e_norm)) > _B);

  // encrypt
  vec_ZZ cipher;
  vec_ZZ_p aux;
  cipher = conv<vec_ZZ>(s);
  mul(aux, s, A0);
  cipher.append(conv<vec_ZZ>(aux));
  mul(aux, s, Au);
  cipher.append(conv<vec_ZZ>(aux));
  add(cipher, cipher, e);
  add(cipher, cipher, encode(mu));

  c = conv<vec_ZZ_p>(cipher);
}

void Encryption::Labeled::encrypt(vec_ZZ_p& u,
                                  vec_ZZ_p& c,
                                  const mat_ZZ_p A0,
                                  const mat_ZZ_p A1,
                                  const uint8_t mu) {
  mat_ZZ_p Au;
  u = taggen(Au, A0, A1);
  tag_encrypt(c, Au, A0, A1, mu);
}

void Encryption::Labeled::decrypt(uint8_t& mu,
                                  const vec_ZZ_p u,
                                  const vec_ZZ_p c,
                                  const mat_ZZ R,
                                  const mat_ZZ_p A0,
                                  const mat_ZZ_p A1) {
  vec_ZZ e2;
  vec_ZZ e;
  vec_ZZ_p s;
  _trapdoor->invert(s, e2, R, A0, A1, this->convolution(u), 2 * c);

  e.SetLength(e2.length());
  for (uint32_t i = 0; i < e.length() - 1; i++) {
    if (e2[i] % 2 != 0) {
      panic("decrypt() - cannot get error");
    }
    e[i] = e2[i] / 2;
  }
  e[e.length() - 1] = e2[e.length() - 1] / 2;

  ZZ e_norm;
  InnerProduct(e_norm, e, e);
  if (SqrRoot(conv<RR>(e_norm)) > _B_prime) {
    panic("decrypt() - error too big");
  }
  mu = e2[e.length() - 1] % 2;
}