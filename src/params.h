#ifndef PARAMS_H
#define PARAMS_H

#include <cstdint>

#define MODULUS 2147494753
#define DIMENSION 976
#define GAUSS_WIDTH 41
#define Q_BITS 32

#define MAX_Q_BITS 64                      // Q_BITS = k
#define RANDOMIZED_ROUNDING_PARAMETER 4.5  // r
#define SIGMA_G 5                          // Used in Gauss preimage sampling

extern const uint32_t F_COEFFS_128_2053[];  // toy example

#endif