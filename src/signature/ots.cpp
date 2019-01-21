#include "signature/ots.h"

Signature::OTS::OTS(const uint32_t dimension,
                    const uint32_t modulus_len,
                    Hash::Chameleon* ch) {
  this->n = dimension;
  this->k = modulus_len;
  this->m = (modulus_len + 2) * dimension;

  this->_ch = ch;
}

void Signature::OTS::key_gen(VerifKey& vk, PrivKey& sk) {
  vec_ZZ zero;
  zero.SetLength(k);

  _ch->key_gen(vk.ck0, sk.tk0);
  _ch->key_gen(vk.ck1, sk.tk1);

  _ch->hash(vk.ck0, zero, sk.z0, sk.r0);
  _ch->hash(vk.ck1, conv<vec_ZZ>(sk.z0), vk.z1, sk.r1);
}

void Signature::OTS::sign(const PrivKey sk,
                          const VerifKey vk,
                          const vec_ZZ msg,
                          Sig& sig) {
  sig.r1 = sk.r1;
  _ch->coll(vk.ck0, sk.tk0, sk.z0, msg, sig.r0_prime);
}

bool Signature::OTS::verify(const VerifKey vk,
                            const vec_ZZ msg,
                            const Sig sig) {
  vec_ZZ_p z0_prime;
  _ch->commit(vk.ck0, msg, z0_prime, sig.r0_prime);
  return _ch->verify(vk.ck1, conv<vec_ZZ>(z0_prime), vk.z1, sig.r1);
}
