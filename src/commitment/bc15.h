#ifndef BC15_H
#define BC15_H 1

#include "encryption/labeled.h"
#include "hash/chameleon.h"

namespace Commitment {
/**
 * Implementation of the SPHF-friendly commitment scheme described in
 * [ACNS:BlaChe15]
 */
class BC15 {
 public:
  BC15(const uint32_t dimension,
       const uint32_t modulus_len,
       const uint32_t modulus,
       const double_t rounding,
       Trapdoor::MP12* mp12);

  void setupCom();
  void commit();

 private:
  uint32_t _n;
  uint32_t _k;
  uint32_t _q;
  double_t _r;

  double_t _sigma;
  Trapdoor::MP12* _trapdoor;

  Encryption::Labeled* _enc;
  Hash::Chameleon* _ch;
};
}  // namespace Commitment

#endif  // BC15_H