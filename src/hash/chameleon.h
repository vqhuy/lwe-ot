#ifndef CHAMELEON_H
#define CHAMELEON_H 1

#include "trapdoor/mp12.h"

namespace Hash {
/**
 * For efficiency, the second component of the hash key A1 is splitted into
 * two matrices A1_0, A1_1 such that A1 = [I | A1_0 | A1_1].
 */
struct HashKey {
  mat_ZZ_p A0;    // A0 \in Z_q^{n*k}
  mat_ZZ_p A1_0;  // A1_0 \in Z_q^{n*n}
  mat_ZZ_p A1_1;  // A1_1 \in Z_q^{n*nk}
};

class Chameleon {
 public:
  Chameleon(const uint32_t dimension,
            const uint32_t modulus_len,
            const double_t sigma,
            const double_t rounding,
            Trapdoor::MP12* mp12);

  void key_gen(HashKey& ck, mat_ZZ& R1);  // use the MP12's LWE-based trapdoor
  void verif_key_gen(void);               // not implemented
  void hash(const Hash::HashKey ck,       // ignore verification key `vk`
            const vec_ZZ msg,             //
            vec_ZZ_p& commit,             //
            vec_ZZ& open);                //
  void coll(const Hash::HashKey ck,       //
            const mat_ZZ tk,              //
            const vec_ZZ_p commit,        //
            const vec_ZZ msg1,            //
            vec_ZZ& open1);               //
  bool verify(const Hash::HashKey ck,     // ignore verification key `vtk`
              const vec_ZZ msg,           //
              const vec_ZZ_p commit,      //
              const vec_ZZ open);         //

  // compute commitment value, used in both `hash` and `verify` function
  void commit(const Hash::HashKey ck,
              const vec_ZZ msg,
              vec_ZZ_p& commit,
              const vec_ZZ open);

 private:
  Trapdoor::MP12* trapdoor;
  uint32_t m;
  uint32_t n;
  uint32_t k;
  double_t s;
  double_t r;
};
}  // namespace Hash

#endif  // CHAMELEON_H
