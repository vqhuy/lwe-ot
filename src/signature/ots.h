#ifndef OTS_H
#define OTS_H 1

#include "hash/chameleon.h"

namespace Signature {

struct VerifKey {
  Hash::HashKey ck0;
  Hash::HashKey ck1;
  vec_ZZ_p z1;
};

struct PrivKey {
  mat_ZZ tk0;
  mat_ZZ tk1;
  vec_ZZ_p z0;
  vec_ZZ r0;
  vec_ZZ r1;
};

struct Sig {
  vec_ZZ r0_prime;
  vec_ZZ r1;
};

class OTS {
 public:
  OTS(const uint32_t dimension,
      const uint32_t modulus_len,
      Hash::Chameleon* ch);

  void key_gen(VerifKey& vk, PrivKey& sk);
  void sign(const PrivKey sk, const VerifKey vk, const vec_ZZ msg, Sig& sig);
  bool verify(const VerifKey vk, const vec_ZZ msg, const Sig sig);

 private:
  Hash::Chameleon* _ch;
  uint32_t m;
  uint32_t n;
  uint32_t k;
};
}  // namespace Signature

#endif  // OTS_H