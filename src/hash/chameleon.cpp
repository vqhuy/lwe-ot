#include "hash/chameleon.h"
#include <NTL/BasicThreadPool.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vector.h>
#include <stdexcept>

Hash::Chameleon::Chameleon(const uint32_t dimension,
                           const uint32_t modulus_len,
                           const double_t sigma,
                           const double_t rounding,
                           Trapdoor::MP12* mp12) {
  this->n = dimension;
  this->k = modulus_len;
  this->m = (modulus_len + 2) * dimension;
  this->s = sigma * rounding;  // Gaussian param of hash() and coll()
  this->r = rounding;
  this->trapdoor = mp12;
}

void Hash::Chameleon::key_gen(Hash::HashKey& ck, mat_ZZ& R1) {
  // pick A0 at random
  // cheating: set message length to be k
  ck.A0.SetDims(n, k);
  generate_random_matrix(ck.A0);

  // A1, R1 <- GenTrap()
  trapdoor->generate_1(R1, ck.A1_0, ck.A1_1);
}

void Hash::Chameleon::verif_key_gen(void) {
  throw std::logic_error("Function not yet implemented");
}

void Hash::Chameleon::commit(const Hash::HashKey ck,
                             const vec_ZZ msg,
                             vec_ZZ_p& commit,
                             const vec_ZZ open) {
  vec_ZZ_p r1, r2, r3;
  r1.SetLength(n);
  NTL_EXEC_RANGE(n, first, last)
  for (uint32_t i = first; i < last; i++) {
    r1[i] = conv<ZZ_p>(open[i]);
  }
  NTL_EXEC_RANGE_END
  r2.SetLength(n);
  NTL_EXEC_RANGE(n, first, last)
  for (uint32_t i = first; i < last; i++) {
    r2[i] = conv<ZZ_p>(open[n + i]);
  }
  NTL_EXEC_RANGE_END
  r3.SetLength(n * k);
  NTL_EXEC_RANGE(n * k, first, last)
  for (uint32_t i = first; i < last; i++) {
    r3[i] = conv<ZZ_p>(open[2 * n + i]);
  }
  NTL_EXEC_RANGE_END

  vec_ZZ_p aux;

  mul(aux, ck.A1_0, r2);
  add(commit, r1, aux);
  mul(aux, ck.A1_1, r3);
  add(commit, commit, aux);

  mul(aux, ck.A0, conv<vec_ZZ_p>(msg));
  add(commit, commit, aux);
}

void Hash::Chameleon::hash(const Hash::HashKey ck,
                           const vec_ZZ msg,
                           vec_ZZ_p& commit,
                           vec_ZZ& open) {
  open.SetLength(m);
  NTL_EXEC_RANGE(m, first, last)
  for (uint32_t i = first; i < last; i++) {
    open[i] = trapdoor->sample(s, 0);
  }
  NTL_EXEC_RANGE_END

  this->commit(ck, msg, commit, open);
}

void Hash::Chameleon::coll(const Hash::HashKey ck,
                           const mat_ZZ tk,
                           const vec_ZZ_p commit,
                           const vec_ZZ msg1,
                           vec_ZZ& open1) {
  vec_ZZ_p aux;
  mul(aux, ck.A0, conv<vec_ZZ_p>(msg1));
  vec_ZZ_p u;
  sub(u, commit, aux);
  open1 = trapdoor->sampleD(tk, u);
}

bool Hash::Chameleon::verify(const Hash::HashKey ck,
                             const vec_ZZ msg,
                             const vec_ZZ_p commit,
                             const vec_ZZ open) {
  ZZ norm2;
  InnerProduct(norm2, open, open);
  if (SqrRoot(conv<RR>(norm2)) > s * r * sqrt(m)) {
    return false;
  }

  vec_ZZ_p verif;
  this->commit(ck, msg, verif, open);
  if (!IsZero(verif - commit)) {
    return false;
  }

  return true;
}