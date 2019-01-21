#include "mp12.h"

using namespace std;
NTL_CLIENT

struct Trapdoor::MP12::encapsulation {
  const uint32_t dimension;    // n
  const uint32_t modulus;      // q
  const uint32_t modulus_len;  // k
  const uint32_t m_bar;        // m_bar = 2 * n
  const uint32_t omega;        // omega = n * k
  const uint32_t m;            // m = m_bar + omega
  const double_t r;            // r
  const double_t sigma;        // Gaussian parameter

  const uint8_t tau = 12;  // tailcut parameter, default is set to 12.
  dgs_rround_dp_t* sampler;
  int64_t sampleZ(const double_t s, const double_t c);

  std::mt19937_64* randomGenerator;
  std::unique_ptr<std::mt19937_64> randomGenerator_t;

  // Variables associated to sampleG
  // const int32_t G_b = 2;
  uint8_t modulus_bits[MAX_Q_BITS];
  double_t G_l[MAX_Q_BITS];
  double_t G_h[MAX_Q_BITS];
  double_t G_d[MAX_Q_BITS];

  // Functions associated to sampleG
  void sampleG_offline(const double_t s);
  int32_t* G_perturb(void);
  vec_ZZ G_sampleD(const double_t s, const double_t* c);
  vec_ZZ G_sampleG(const double_t s, const ZZ_p coset_i);
  vec_ZZ sampleG(const double_t s, const vec_ZZ_p coset);

  // Variables associated to samplePre
  mat_RR param_P;    // sqrt(sigma_P - r^2)
  double_t param_G;  // r * sqrt(sigma_G)

  // Functions associated to samplePre
  void samplePre_init(const mat_ZZ R, const double_t s);
  void samplePre_offline(const mat_ZZ_p A0, const mat_ZZ R);
  vec_ZZ samplePre(const mat_ZZ R, const vec_ZZ_p u);

  mat_ZZ_p gadget;  // matrix G with g = [1, 2, 4, ... 2^{k-1}]

  std::stack<int32_t*> G_precomputedPerturbations;
  std::stack<vec_ZZ> precomputed_p;
  std::stack<vec_ZZ_p> precomputed_w;
  std::stack<vec_ZZ_p> precomputed_w_bar;
};

int64_t Trapdoor::MP12::encapsulation::sampleZ(const double_t s,
                                               const double_t c) {
  return sampler->call(sampler, s, c);
}

void Trapdoor::MP12::encapsulation::sampleG_offline(const double_t s) {
  int32_t* perturbation{new int32_t[modulus_len]};
  double_t sigmas[modulus_len];
  double_t beta = 0.0;
  int32_t z[modulus_len];

  for (uint32_t i = 0; i < modulus_len; ++i) {
    sigmas[i] = s / G_l[i];
    double_t c = beta / G_l[i];
    int32_t cz = (int32_t)floor(c);
    z[i] = cz + sampleZ(sigmas[i], c - cz);
    beta = -z[i] * G_h[i];
  }

  perturbation[0] = z[0] + (z[0] << 2) + (z[1] << 1);  // b = 2
  for (uint32_t i = 1; i < modulus_len; ++i) {
    perturbation[i] = (z[i - 1] + (z[i] << 1) + z[i + 1]) << 1;
  }
  // add to the precomputed perturbations
  G_precomputedPerturbations.push(perturbation);
}

int32_t* Trapdoor::MP12::encapsulation::G_perturb(void) {
  int32_t* result = G_precomputedPerturbations.top();
  G_precomputedPerturbations.pop();
  return result;
}

vec_ZZ Trapdoor::MP12::encapsulation::G_sampleD(const double_t s,
                                                const double_t* c) {
  const double_t c_d = -c[modulus_len - 1] / G_d[modulus_len - 1];
  const int32_t c_dz = (int32_t)floor(c_d);
  vec_ZZ output;
  output.SetLength(modulus_len);
  output[modulus_len - 1] =
      c_dz + sampleZ(s / G_d[modulus_len - 1], c_d - c_dz);

  double_t cs[modulus_len];
  const int32_t zk_1 = conv<int32_t>(output[modulus_len - 1]);
  for (uint32_t i = 0; i < modulus_len; ++i) {
    cs[i] = zk_1 * G_d[i] - c[i];
  }

  for (uint32_t i = 0; i < modulus_len - 1; ++i) {
    const int32_t cs_iz = (int32_t)floor(cs[i]);
    output[i] = cs_iz + sampleZ(s, cs[i] - cs_iz);
  }

  return output;
}

vec_ZZ Trapdoor::MP12::encapsulation::G_sampleG(const double_t s,
                                                const ZZ_p coset_i) {
  vec_ZZ output;
  output.SetLength(modulus_len);

  // get bits of coset_i
  uint8_t u[modulus_len];
  getBits(u, conv<ZZ>(coset_i), modulus_len);

  // get perturbations
  int32_t* p = G_perturb();

  // compute cs
  double_t cs[modulus_len];
  cs[0] = (u[0] - p[0]) / 2;
  for (uint32_t i = 1; i < modulus_len; ++i) {
    cs[i] = (cs[i - 1] + u[i] - p[i]) / 2;
  }

  vec_ZZ z = G_sampleD(s, cs);

  // compute t
  output[0] = (z[0] << 1) + u[0];
  for (uint32_t i = 1; i < modulus_len - 1; ++i) {
    output[i] = (z[i] << 1) - z[i - 1] + modulus_bits[i] * z[i - 1] + u[i];
  }
  output[modulus_len - 1] = modulus_bits[modulus_len - 1] * z[modulus_len - 1] -
                            z[modulus_len - 2] + u[modulus_len - 1];

  delete[] p;
  return output;
}

vec_ZZ Trapdoor::MP12::encapsulation::sampleG(const double_t s,
                                              const vec_ZZ_p coset) {
  vec_ZZ output;
  for (uint32_t i = 0; i < dimension; ++i) {
    output.append(G_sampleG(s, coset[i]));
  }
  return output;
}

void Trapdoor::MP12::encapsulation::samplePre_init(const mat_ZZ R,
                                                   const double_t s) {
  // compute param_P:
  // param_P = sqrt(sigma_P - r^2)
  // sigma_P = s^2 * I - r^2 * sigma_G * [R; I] * [R.T | I]
  //         = s^2 * I - r^2 * sigma_G * [ R * R.T | R ]
  //                                     [ R.T     | I ]
  // sigma_G is fixed to 5

  const double_t s2 = s * s;

  mat_ZZ RT = transpose(R);
  mat_ZZ COV, cov1, cov2;
  mat_ZZ aux;
  mul(aux, R, RT);
  cov1 = mat_concat_horizontal(aux, R);
  cov2 = mat_concat_horizontal(RT, ident_mat_ZZ(omega));
  COV = mat_concat_vertical(cov1, cov2);

  mat_RR sigma_P;
  sigma_P = s2 * ident_mat_RR(m) - SIGMA_G * conv<mat_RR>(COV);
  param_P = r * cholesky(sigma_P);
}

void Trapdoor::MP12::encapsulation::samplePre_offline(const mat_ZZ_p A0,
                                                      const mat_ZZ R) {
  vec_RR d;
  vec_RR b;       // b = B_2 * D_1^m
  vec_ZZ p1, p2;  // p = [p1; p2]

  vec_ZZ w_bar, aux, aux1, aux2;
  vec_ZZ_p w;

  // sample p from D(Z^m, sigma_P)
  // - sample from continuous Gaussian over R^m centered on 0 and of std dev 1
  // - sample from discrete Gauss over Z_m centered on each entry w/ param r
  d.SetLength(m);
  std::normal_distribution<> normalDist(0, 1);
  for (uint32_t i = 0; i < m; i++) {
    d[i] = normalDist(*randomGenerator);
  }
  b = param_P * d;

  p1.SetLength(m_bar);
  p2.SetLength(omega);
  for (uint32_t i = 0; i < m_bar; i++) {
    p1[i] = sampleZ(r, conv<double_t>(b[i]));
  }
  for (uint32_t i = 0; i < omega; i++) {
    p2[i] = sampleZ(r, conv<double_t>(b[i + m_bar]));
  }

  // offline phase computations
  // w_bar = [I | A0] * (p1 - R * p2)
  //       = [I | A0] * [AUX1; AUX2]
  //       = [AUX1 + A0 * AUX2]
  aux = p1 - R * p2;
  aux1.SetLength(dimension);
  aux2.SetLength(dimension);
  for (uint32_t i = 0; i < dimension; i++) {
    aux1[i] = aux[i];
  }
  for (uint32_t i = 0; i < dimension; i++) {
    aux2[i] = aux[i + dimension];
  }

  w_bar = aux1 + conv<mat_ZZ>(A0) * aux2;
  w = gadget * conv<vec_ZZ_p>(p2);

  // add to the precomputed perturbations
  precomputed_w_bar.push(conv<vec_ZZ_p>(w_bar));
  precomputed_w.push(w);
  p1.append(p2);
  precomputed_p.push(p1);
}

// TODO(1): update this to accept input tag H
// thus all the trapdoor computations need to be updated accordingly
vec_ZZ Trapdoor::MP12::encapsulation::samplePre(const mat_ZZ R,
                                                const vec_ZZ_p u) {
  vec_ZZ_p w_bar_q;
  vec_ZZ_p w_q;
  vec_ZZ_p v;
  vec_ZZ p;
  vec_ZZ z;
  vec_ZZ aux;  // aux = [R; I] * z

  w_bar_q = precomputed_w_bar.top();
  precomputed_w_bar.pop();
  w_q = precomputed_w.top();
  precomputed_w.pop();
  v = u - w_bar_q - w_q;

  z = sampleG(param_G, v);
  aux = R * z;
  aux.append(z);
  p = precomputed_p.top();
  precomputed_p.pop();
  return p + aux;
}

Trapdoor::MP12::MP12(const uint32_t n,      // dimension
                     const uint32_t k,      // modulus_len
                     const uint32_t q,      // modulus
                     const double_t sigma,  // Gaussian parameter
                     const double_t r       // rounding parameter
                     )
    : impl(new encapsulation{.dimension = n,
                             .modulus = q,
                             .modulus_len = k,
                             .m_bar = 2 * n,
                             .omega = n * k,
                             .m = (k + 2) * n,
                             .r = r,
                             .sigma = sigma}) {
  impl->param_G =
      impl->r * sqrt(SIGMA_G);  // choose sigma_G = 5 if modulus is odd prime
  impl->sampler = dgs_rround_dp_init(impl->tau, DGS_RROUND_DEFAULT);

  // random generator for continuous Gauss
  std::random_device rd;
  impl->randomGenerator_t.reset(new std::mt19937_64(rd()));
  impl->randomGenerator = impl->randomGenerator_t.get();

  // calculate gadget = G
  vec_ZZ_p g;
  g.SetLength(k);
  for (uint32_t i = 0; i < k; i++) {
    g[i] = 1 << i;
  }

  impl->gadget.SetDims(n, n * k);
  for (uint32_t i = 0; i < n; i++) {
    for (uint32_t j = 0; j < k; j++) {
      impl->gadget[i][j + i * k] = g[j];
    }
  }

  // compute binary representation of modulus
  getBits(impl->modulus_bits, ZZ(q), k);

  // compute G_l, G_h and G_d
  impl->G_l[0] = sqrt(2 * (1 + 1 / k) + 1);
  impl->G_h[0] = 0;
  for (uint32_t i = 1; i < k; ++i) {
    impl->G_l[i] = sqrt(2 * (1 + 1 / (k - i)));
    impl->G_h[i] = sqrt(2 * (1 - 1 / (k - i + 1)));
  }
  impl->G_d[0] = impl->modulus_bits[0] / 2;
  for (uint32_t i = 1; i < k; ++i) {
    impl->G_d[i] = (impl->G_d[i - 1] + impl->modulus_bits[i]) / 2;
  }

  this->init_invert();
}

Trapdoor::MP12::~MP12(void) {
  dgs_rround_dp_clear(impl->sampler);
}

int64_t Trapdoor::MP12::sample(const double_t s, const double_t c) {
  return impl->sampleZ(s, c);
}

vec_ZZ Trapdoor::MP12::sampleD(const mat_ZZ R, const vec_ZZ_p u) {
  return impl->samplePre(R, u);
}

void Trapdoor::MP12::preCompute(const mat_ZZ R,
                                const mat_ZZ_p A0,
                                const mat_ZZ_p A1,
                                const double_t s,
                                const uint32_t iter) {
  impl->samplePre_init(R, s);
  for (uint32_t i = 0; i < iter; i++) {
    impl->sampleG_offline(s);
    impl->samplePre_offline(A0, R);
  }
}

// TODO: remove this function, and use generate_0 with ident matrix
// since we fixed to make invert() works with generate_0 only.
void Trapdoor::MP12::generate_1(mat_ZZ& R, mat_ZZ_p& A0, mat_ZZ_p& A1) {
  mat_ZZ R1, R2;

  // pick A0 at random
  A0.SetDims(impl->dimension, impl->dimension);
  generate_random_matrix(A0);

  // sample R from D(Z^{m_bar * omega}) = D(Z^{2n * nk})
  R1.SetDims(impl->dimension, impl->omega);
  R2.SetDims(impl->dimension, impl->omega);
  NTL_EXEC_RANGE(impl->dimension, first, last)
  for (uint32_t i = first; i < last; i++)
    for (uint32_t j = 0; j < impl->omega; j++) {
      R1[i][j] = impl->sampleZ(impl->sigma, 0);
      R2[i][j] = impl->sampleZ(impl->sigma, 0);
    }
  NTL_EXEC_RANGE_END

  // compute A
  conv(A1, R2);
  mul(A1, A0, A1);
  add(A1, A1, conv<mat_ZZ_p>(R1));
  sub(A1, impl->gadget, A1);

  R = mat_concat_vertical(R1, R2);
}

void Trapdoor::MP12::init_invert(void) {
  // compute the basis Sk of the dual lattice L*(g^T).
  // See [MP12, sect. 4.2]
  mat_ZZ Sk;
  vec_ZZ col;
  Sk.SetDims(impl->modulus_len, impl->modulus_len);
  col.SetLength(impl->modulus_len);
  col[0] = 2;
  col[1] = -1;
  for (uint32_t i = 0; i < impl->modulus_len - 1; i++) {
    Sk[i] = col;
    std::rotate(col.begin(), col.end() - 1, col.end());
  }
  Sk = transpose(Sk);
  for (uint32_t i = 0; i < impl->modulus_len; i++) {
    Sk[i][impl->modulus_len - 1] = impl->modulus_bits[i];
  }

  // compute the basis Vk of the lattice L(g)
  // Vk = q * (Sk^{-1})
  Vk = conv<mat_ZZ>(impl->modulus * inv(conv<mat_RR>(Sk)));

  // compute the Gram-Schmidth data
  gso_Vk_star.SetDims(impl->modulus_len, impl->modulus_len);
  mat_RR gso_mu;
  ComputeGS(Vk, gso_mu, gso_c);
  gso_Vk_star(1) = conv<vec_RR>(Vk(1));
  for (uint32_t i = 2; i <= impl->modulus_len; i++) {
    gso_Vk_star(i) = conv<vec_RR>(Vk(i));
    for (uint32_t j = 1; j < i; j++) {
      gso_Vk_star(i) -= gso_mu(i, j) * gso_Vk_star(j);
    }
  }
}

void Trapdoor::MP12::babai_nearest_plane(ZZ_p& s, const vec_ZZ_p t) {
  // inverse g_g := s*g + e = t with the basis Vk
  vec_ZZ small = conv<vec_ZZ>(t);
  for (int32_t i = impl->modulus_len - 1; i >= 0; i--) {
    ZZ c = RoundToZZ((conv<vec_RR>(small) * gso_Vk_star[i]) / gso_c[i]);
    small -= Vk[i] * c;
  }
  vec_ZZ_p lp = t - conv<vec_ZZ_p>(small);
  s = lp[0];
}

void Trapdoor::MP12::invert_G(vec_ZZ_p& s, const vec_ZZ_p gA) {
  s.SetLength(impl->dimension);
  NTL_EXEC_RANGE(impl->dimension, first, last)
  for (uint32_t i = first; i < last; i++) {
    ZZ_p small_s;
    vec_ZZ_p target;
    target.SetLength(impl->modulus_len);
    for (uint32_t j = 0; j < impl->modulus_len; j++) {
      target[j] = gA[i * impl->modulus_len + j];
    }
    this->babai_nearest_plane(small_s, target);
    s[i] = small_s;
  }
  NTL_EXEC_RANGE_END
}

void Trapdoor::MP12::invert(vec_ZZ_p& s,
                            vec_ZZ& e,
                            const mat_ZZ R,
                            const mat_ZZ_p A0,
                            const mat_ZZ_p A1,
                            const mat_ZZ_p H,
                            const vec_ZZ_p b) {
  // compute LWE instance of the gadgat matrix G
  vec_ZZ_p b1, b2;
  vec_ZZ_p b_hat;
  vec_ZZ_p s_hat;
  vec_ZZ e_hat;

  b1.SetLength(impl->m_bar);
  for (uint32_t i = 0; i < impl->m_bar; i++) {
    b1[i] = b[i];
  }
  b2.SetLength(impl->omega);
  for (uint32_t i = 0; i < impl->omega; i++) {
    b2[i] = b[impl->m_bar + i];
  }
  mul(b_hat, b1, conv<mat_ZZ_p>(R));
  add(b_hat, b_hat, b2);

  // call the oracle
  this->invert_G(s_hat, b_hat);

  mul(s, transpose(inv(H)), s_hat);

  // note: re-compute parity-check matrix A with given tag H
  mat_ZZ_p auxA;
  mul(auxA, H, impl->gadget);
  add(auxA, auxA, A1);

  vec_ZZ_p aux;
  e = conv<vec_ZZ>(s);
  mul(aux, s, A0);
  e.append(conv<vec_ZZ>(aux));
  mul(aux, s, auxA);
  e.append(conv<vec_ZZ>(aux));
  sub(e, conv<vec_ZZ>(b), e);
}

void Trapdoor::MP12::generate_0(mat_ZZ& R, mat_ZZ_p& A0, mat_ZZ_p& A1) {
  mat_ZZ R1, R2;

  // pick A0 at random
  A0.SetDims(impl->dimension, impl->dimension);
  generate_random_matrix(A0);

  // sample R from D(Z^{m_bar * omega}) = D(Z^{2n * nk})
  R1.SetDims(impl->dimension, impl->omega);
  R2.SetDims(impl->dimension, impl->omega);
  NTL_EXEC_RANGE(impl->dimension, first, last)
  for (uint32_t i = first; i < last; i++)
    for (uint32_t j = 0; j < impl->omega; j++) {
      R1[i][j] = impl->sampleZ(impl->sigma, 0);
      R2[i][j] = impl->sampleZ(impl->sigma, 0);
    }
  NTL_EXEC_RANGE_END

  // compute A1 = -(A0 * R2 + R1)
  mat_ZZ_p aux;
  conv(A1, R2);
  mul(A1, A0, A1);
  add(A1, A1, conv<mat_ZZ_p>(R1));
  aux.SetDims(A1.NumRows(), A1.NumCols());
  sub(A1, aux, A1);

  R = mat_concat_vertical(R1, R2);
}

mat_ZZ_p Trapdoor::MP12::gadget(void) {
  return impl->gadget;
}