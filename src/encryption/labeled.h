#ifndef LABELED_H
#define LABELED_H 1

#include <cstdint>
#include "prng/random.h"
#include "trapdoor/mp12.h"

NTL_CLIENT

namespace Encryption {
class Labeled {
 public:
  Labeled(const uint32_t dimension,
          const uint32_t modulus_len,
          const uint32_t modulus,
          const double_t sigma,
          const double_t rounding,
          const vec_ZZ_p f_coeffs,
          Trapdoor::MP12* mp12);

  void keygen(mat_ZZ_p& A0,  // the public key is A = [I | A0 | A1]
              mat_ZZ_p& A1,  // with A1 = A0 * R2 + R1
              mat_ZZ& R      // the secret key is the trapdoor
  );

  vec_ZZ_p taggen(mat_ZZ_p& Au,
                  const mat_ZZ_p A0,
                  const mat_ZZ_p A1);  // generate a fixed tag for tag_encrypt

  void tag_encrypt(vec_ZZ_p& c,        // ciphertext
                   const mat_ZZ_p Au,  // public matrix of the chosen tag h(u)
                   const mat_ZZ_p A0,  // public key
                   const mat_ZZ_p A1,  // (redundant)
                   const uint8_t mu    // m \in {0,1}
  );

  void encrypt(vec_ZZ_p& u,        // to compute the invertible tag H
               vec_ZZ_p& c,        // ciphertext
               const mat_ZZ_p A0,  // public key
               const mat_ZZ_p A1,  //
               const uint8_t mu    // m \in {0,1}
  );

  void decrypt(uint8_t& mu,        // plaintext
               const vec_ZZ_p u,   // tag
               const vec_ZZ_p c,   // ciphertex
               const mat_ZZ R,     // trapdoor - secret key
               const mat_ZZ_p A0,  // public key
               const mat_ZZ_p A1   //
  );

 private:
  uint32_t _n;
  uint32_t _q;
  uint32_t _k;
  uint32_t _m;
  double_t _sigma;
  double_t _r;

  vec_ZZ_p _f_coeffs;

  Trapdoor::MP12* _trapdoor;

  double_t _t;
  double_t _B;
  double_t _B_prime;
  vec_ZZ _encoding;

  vec_ZZ encode(uint8_t mu);

  /**
   * random_element_in_U - return a non-zero random element in U represented as
   * a vector in Z_q^n relative to the standard basis of monomials
   * {1, x, ..., x^{n-1}}
   *
   * a = [a0, a1, ..., a_{n-1}] = a0 + a1*x + ... + a_{n-1}*x^{n-1}
   *
   * take U to be all linear combinations of the monomials 1, x, . . . , x^{n−1}
   * with coefficients in {0, . . . , p − 1}.
   */
  vec_ZZ_p random_element_in_U();

  /**
   * convolution - map an element u in R* to matrix H in Z_q^{n*n} representing
   * the right multiplication by this element on the power basis
   */
  mat_ZZ_p convolution(vec_ZZ_p u);
};
}  // namespace Encryption

#endif  // LABELED_H