#ifndef CCA2_H
#define CCA2_H 1

#include "encryption/labeled.h"
#include "signature/ots.h"

namespace Encryption {
class CCA2 {
 public:
  CCA2(Labeled* enc);
  void keygen(mat_ZZ_p& A0,  // the public key is A = [I | A0 | A1]
              mat_ZZ_p& A1,  // with A1 = A0 * R2 + R1
              mat_ZZ& R      // the secret key is the trapdoor
  );

  void encrypt(vec_ZZ_p& u,        // to compute the invertible tag H
               vec_ZZ_p& c,        // ciphertext
               const mat_ZZ_p A0,  // public key
               const mat_ZZ_p A1,  //
               const uint8_t mu    // m \in {0,1}
  );

 private:
  Labeled* _enc;
};
}  // namespace Encryption

#endif  // CCA2_H