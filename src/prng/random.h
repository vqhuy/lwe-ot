#ifndef RANDOM_H
#define RANDOM_H 1

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <climits>
#include <cstdint>
#include <random>
#include "fastrandombytes.h"

NTL_CLIENT

ZZ_p random_Zp();
ZZ randint(const ZZ bnd);

void generate_random_matrix(
    Mat<ZZ_p>& A);  // generate a uniformly random matrix in Z_p

vec_ZZ rand_vec_bin(uint32_t len);

vec_ZZ_p rand_vec_Zp(uint32_t len);

uint8_t real_random(double_t prob);

#endif  // RANDOM_H