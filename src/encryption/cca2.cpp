#include "encryption/cca2.h"

Encryption::CCA2::CCA2(Labeled* enc) {
  this->_enc = enc;
}

void Encryption::CCA2::keygen(mat_ZZ_p& A0, mat_ZZ_p& A1, mat_ZZ& R) {
  _enc->keygen(A0, A1, R);
}

void Encryption::CCA2::encrypt(vec_ZZ_p& u,
                               vec_ZZ_p& c,
                               const mat_ZZ_p A0,
                               const mat_ZZ_p A1,
                               const uint8_t mu) {
  // cheating: since we just use the OTS once so there's no need to construct a
  // map from the verification keys to the tag set. Mapping the verif key vk to
  // a random tag (aka label) is okay.
}