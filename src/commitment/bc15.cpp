#include "commitment/bc15.h"

Commitment::BC15::BC15(const uint32_t dimension,
                       const uint32_t modulus_len,
                       const uint32_t modulus,
                       const double_t rounding,
                       Trapdoor::MP12* mp12) {
  this->_n = dimension;
  this->_k = modulus_len;
  this->_q = modulus;
  this->_r = rounding;
  this->_trapdoor = mp12;
  this->_sigma = 2 * sqrt(_n);  // set sigma = 2*sqrt(n) as in BBDQ17
                                // to ensure that the scheme is correct

  // Encryption
  vec_ZZ_p f_coeffs;
  f_coeffs.SetLength(_n);
  NTL_EXEC_RANGE(_n, first, last)
  for (uint32_t i = first; i < last; i++) {
    f_coeffs[i] = F_COEFFS_512_37889[i];
  }
  NTL_EXEC_RANGE_END
  this->_enc =
      new Encryption::Labeled(_n, _k, _q, _sigma, _r, f_coeffs, _trapdoor);

  // Chameleon Hash
  this->_ch = new Hash::Chameleon(_n, _k, _sigma, _r, _trapdoor);
}

void Commitment::BC15::setupCom() {}

void Commitment::BC15::commit() {}