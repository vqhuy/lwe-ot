#ifndef MP12_H
#define MP12_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <deque>
#include <memory>
#include <random>
#include <stack>

#include <NTL/BasicThreadPool.h>
#include <NTL/LLL.h>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/matrix.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vector.h>

#include <gmpxx.h>
extern "C" {
#include "dgs/dgs/dgs.h"
}
#include "params.h"
#include "prng/random.h"
#include "utils/utils.h"

NTL_CLIENT

namespace Trapdoor {
/**
 * MP12's trapdoor generator using computational instantiation based on LWE.
 * sampleG: Gaussian preimage sampling based on eprint:2017/308, its
 * implementation is based on the implementation of Bert et al.
 * https://github.com/lbibe/code
 */
class MP12 {
 public:
  MP12(const uint32_t n,
       const uint32_t k,
       const uint32_t q,
       const double_t sigma,
       const double_t r);
  ~MP12(void);

  void preCompute(const mat_ZZ R,
                  const mat_ZZ_p A0,
                  const mat_ZZ_p A1,
                  const double_t s,
                  const uint32_t iter);

  int64_t sample(const double_t s, const double_t c);
  vec_ZZ sampleD(const mat_ZZ R, const vec_ZZ_p u);

  /**
   * TrapGen with tag H equals the identity matrix I.
   * This is the default TrapGen function, and is used in place where the tag H
   * is ignored.
   */
  void generate_1(mat_ZZ& R,     // trapdoor
                  mat_ZZ_p& A0,  // parity-check matrix component,
                  mat_ZZ_p& A1   // with A = [A_bar | A1], A_bar = [I | A0]
  );

  /**
   * TrapGen with tag H equals the zero matrix 0.
   * This is used in the multi-labeled encryption.
   */
  void generate_0(mat_ZZ& R,     // trapdoor
                  mat_ZZ_p& A0,  // parity-check matrix component,
                  mat_ZZ_p& A1   // with A = [A_bar | A1], A_bar = [I | A0]
  );

  void invert(vec_ZZ_p& s,
              vec_ZZ& e,
              const mat_ZZ R,     // trapdoor
              const mat_ZZ_p A0,  // parity-check matrix component
              const mat_ZZ_p A1,  //
              const mat_ZZ_p H,   // invertible tag
              const vec_ZZ_p b    // LWE instance
  );

  mat_ZZ_p gadget(void);

 private:
  // things associated to sampling functions.
  struct encapsulation;
  std::unique_ptr<encapsulation> impl;

  // variables and functions associated to invert()
  mat_ZZ Vk;           // basis of lattice L(g)
  mat_RR gso_Vk_star;  // Gramm-Schmidt data for Vk
  vec_RR gso_c;        // Gramm-Schmidt data for Vk
  void init_invert(void);
  void babai_nearest_plane(ZZ_p& s, const vec_ZZ_p c);
  void invert_G(vec_ZZ_p& s, const vec_ZZ_p gA);
};
}  // namespace Trapdoor

#endif