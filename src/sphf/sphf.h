#ifndef SPHF_H
#define SPHF_H 1

#include "trapdoor/mp12.h"

namespace SPHF {

class NaturalBitPHF {
 public:
  NaturalBitPHF(const uint32_t dimension,
                const uint32_t modulus_len,
                const uint32_t modulus,
                const double_t rounding,
                Trapdoor::MP12* mp12);

  vec_ZZ hashKG();
  vec_ZZ_p projKG(const vec_ZZ h,
                  const mat_ZZ_p A0,
                  const mat_ZZ_p A1,
                  const mat_ZZ_p H);
  uint8_t hash(const vec_ZZ h, const vec_ZZ_p c);
  uint8_t projHash(const vec_ZZ_p p,
                   const vec_ZZ_p c,
                   const vec_ZZ_p s,
                   const vec_ZZ e);

 private:
  uint32_t _n;
  uint32_t _k;
  uint32_t _q;
  double_t _r;

  uint32_t _m;
  double_t _sigma;

  Trapdoor::MP12* _trapdoor;

  uint8_t rounding(ZZ x);
};
}  // namespace SPHF

#endif  // SPHF_H