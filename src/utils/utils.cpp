#include "utils.h"

void panic(const char* err) {
  printf("%s\n", err);
  exit(-1);
}

void getBits(uint8_t* bits, ZZ m, uint32_t length) {
  for (uint32_t i = 0; i < length; ++i) {
    bits[i] = bit(m, i);
  }
}

template <class T>
mat_RR cholesky(const T mat) {
  if (mat.NumCols() != mat.NumRows()) {
    panic("cholesky()");
  }
  mat_RR l = conv<mat_RR>(mat);
  const uint32_t m = mat.NumRows();
  for (uint32_t k = 1; k <= m; k++) {
    for (uint32_t i = 1; i < k; i++) {
      l(i, k) = 0;
    }
    l(k, k) = SqrRoot(l(k, k));
    for (uint32_t i = k + 1; i <= m; i++) {
      l(i, k) = l(i, k) / l(k, k);
    }
    for (uint32_t j = k + 1; j <= m; j++) {
      for (uint32_t i = j; i <= m; i++) {
        l(i, j) = l(i, j) - l(i, k) * l(j, k);
      }
    }
  }
  return l;
}
template mat_RR cholesky<mat_ZZ>(mat_ZZ);
template mat_RR cholesky<mat_RR>(mat_RR);

template <class T>
T mat_concat_vertical(const T mat1, const T mat2) {
  if (mat1.NumCols() != mat2.NumCols()) {
    panic("mat_concat_vertical()");
  }

  uint32_t r1, r2;
  uint32_t c;
  r1 = mat1.NumRows();
  r2 = mat2.NumRows();
  c = mat1.NumCols();

  T res;
  res.SetDims(r1 + r2, c);

  for (uint32_t i = 0; i < r1; i++)
    res[i] = mat1[i];

  for (uint32_t i = 0; i < r2; i++)
    res[r1 + i] = mat2[i];

  return res;
}

template mat_ZZ_p mat_concat_vertical<mat_ZZ_p>(mat_ZZ_p, mat_ZZ_p);
template mat_ZZ mat_concat_vertical<mat_ZZ>(mat_ZZ, mat_ZZ);
template mat_RR mat_concat_vertical<mat_RR>(mat_RR, mat_RR);

template <class T>
T mat_concat_horizontal(const T mat1, const T mat2) {
  if (mat1.NumRows() != mat2.NumRows()) {
    panic("mat_concat_horizontal()");
  }

  uint32_t r;
  uint32_t c1, c2;
  r = mat1.NumRows();
  c1 = mat1.NumCols();
  c2 = mat2.NumCols();

  T res;
  res._mat__rep.SetLength(r);
  res._mat__numcols = c1 + c2;

  for (uint32_t i = 0; i < r; i++) {
    res[i].append(mat1[i]);
    res[i].append(mat2[i]);
  }

  return res;
}

template mat_ZZ_p mat_concat_horizontal<mat_ZZ_p>(mat_ZZ_p, mat_ZZ_p);
template mat_ZZ mat_concat_horizontal<mat_ZZ>(mat_ZZ, mat_ZZ);
template mat_RR mat_concat_horizontal<mat_RR>(mat_RR, mat_RR);

template <class T>
RR frobenious_norm(const T mat) {
  mat_RR B;
  conv(B, mat);

  RR innerp, normB;
  normB = to_RR(0);

  for (int i = 0; i < B.NumRows(); i++) {
    InnerProduct(innerp, B[i], B[i]);
    normB += innerp;
  }

  return SqrRoot(normB);
}

template RR frobenious_norm<mat_RR>(mat_RR);
template RR frobenious_norm<mat_ZZ>(mat_ZZ);