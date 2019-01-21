#ifndef UTILS_H
#define UTILS_H

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <cstdint>

NTL_CLIENT
using namespace std;

extern "C" {
extern void panic(const char* err);
}

extern "C" {
extern void getBits(uint8_t* bits, ZZ m, uint32_t length);
}

template <class T>
mat_RR cholesky(const T mat);

template <class T>
T mat_concat_vertical(const T mat1, const T mat2);

template <class T>
T mat_concat_horizontal(const T mat1, const T mat2);

template <class T>
RR frobenious_norm(const T mat);

#endif