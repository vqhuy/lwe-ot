#include "sphf/sphf.h"
#include <NTL/BasicThreadPool.h>
#include <algorithm>
#include <cmath>
#include "prng/random.h"

SPHF::NaturalBitPHF::NaturalBitPHF(const uint32_t dimension,
                                   const uint32_t modulus_len,
                                   const uint32_t modulus,
                                   const double_t rounding,
                                   Trapdoor::MP12* mp12) {
  this->_n = dimension;
  this->_k = modulus_len;
  this->_q = modulus;
  this->_r = rounding;
  this->_trapdoor = mp12;

  this->_m = dimension * (modulus_len + 2);
  // calculate sigma according to Corollary 3.2 - BBDQ17
  this->_sigma = _m;  // _sigma = \Theta(m)
}

vec_ZZ SPHF::NaturalBitPHF::hashKG() {
  vec_ZZ hk;
  hk.SetLength(_m);
  NTL_EXEC_RANGE(_m, first, last)
  for (uint32_t i = first; i < last; i++) {
    hk[i] = _trapdoor->sample(_sigma, 0);
  }
  NTL_EXEC_RANGE_END
  return hk;
}

vec_ZZ_p SPHF::NaturalBitPHF::projKG(const vec_ZZ h,
                                     const mat_ZZ_p A0,
                                     const mat_ZZ_p A1,
                                     const mat_ZZ_p H) {
  vec_ZZ_p hp;

  mat_ZZ_p auxA;
  mul(auxA, H, _trapdoor->gadget());
  add(auxA, auxA, A1);

  vec_ZZ_p aux;
  hp = conv<vec_ZZ_p>(h);
  mul(aux, conv<vec_ZZ_p>(h), A0);
  hp.append(aux);
  mul(aux, conv<vec_ZZ_p>(h), auxA);
  hp.append(aux);

  return hp;
}

uint8_t SPHF::NaturalBitPHF::hash(const vec_ZZ h, const vec_ZZ_p c) {
  ZZ inner;
  InnerProduct(inner, h, conv<vec_ZZ>(c));
  return this->rounding(inner);
}

uint8_t SPHF::NaturalBitPHF::projHash(const vec_ZZ_p p,
                                      const vec_ZZ_p c,
                                      const vec_ZZ_p s,
                                      const vec_ZZ e) {
  ZZ inner;
  InnerProduct(inner, conv<vec_ZZ>(p), conv<vec_ZZ>(s));
  return this->rounding(inner);
}

uint8_t SPHF::NaturalBitPHF::rounding(ZZ x) {
  double_t dx = conv<double_t>(x);
  double_t prob = 1.0 / 2.0 + cos((2.0 * M_PI * dx) / _q) / 2.0;
  return real_random(prob);
}