#ifndef HYPERTRITON3STRUCTURES_H
#define HYPERTRITON3STRUCTURES_H

#include <Rtypes.h>

struct RHyperTriton {
  virtual ~RHyperTriton() = default;
  float centrality = -1.;
  float pt = -999.f;
  float phi = -999.f;
  float pz = -999.f;
  float ct = -1.f;
  float r = -1.f;
  float cosPA = -2.f;
  float m = -1;
  float cosPA_Lambda = -2.; 
  Double32_t cosTheta_LambdaD = 1.;   //[-1,1,8]
  Double32_t cosTheta_ProtonPiL = 1.; //[-1,1,8]
  Double32_t cosTheta_ProtonPiH = 1.; //[-1,1,8]
  Double32_t mppi_vert = -1.; //[1.077,1.203,8]
  Double32_t dca_lambda_hyper = -1.0; //[0.0,8.0,8]
  Double32_t dca_de = -1.0; //[0.0,8.0,8]
  Double32_t dca_pr = -1.0; //[0.0,8.0,8]
  Double32_t dca_pi = -1.0; //[0.0,8.0,8]
  Double32_t tpcNsig_de = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_pr = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_pi = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_de = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_pr = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_pi = -4.0; //[-4.0,4.0,8]
  Double32_t dca_de_pr = -4.0; //[0.0,8.0,8]
  Double32_t dca_de_pi = -4.0; //[0.0,8.0,8]
  Double32_t dca_pr_pi = -4.0; //[0.0,8.0,8]
  UChar_t tpcClus_de = 0u;
  UChar_t tpcClus_pr = 0u;
  UChar_t tpcClus_pi = 0u;
  UChar_t candidates = 0u;
  UChar_t trigger = 0u;
  bool hasTOF_de = false;
  bool hasTOF_pr = false;
  bool hasTOF_pi = false;
  bool positive = false;
  ClassDef(RHyperTriton, 1)
};

struct RHyperTriton3O2 : public RHyperTriton {
  RHyperTriton3O2() : RHyperTriton{} {}
  virtual ~RHyperTriton3O2() = default;
  Double32_t dca_de_sv = -4.0; //[0.0,8.0,8]
  Double32_t dca_pr_sv = -4.0; //[0.0,8.0,8]
  Double32_t dca_pi_sv = -4.0; //[0.0,8.0,8]
  Double32_t chi2 = -1.f;      //[0.0,16.,16]
  ClassDef(RHyperTriton3O2,1)
};

struct RHyperTriton3KF : public RHyperTriton {
  RHyperTriton3KF() : RHyperTriton{} {}
  virtual ~RHyperTriton3KF() = default;
  float chi2_deuprot = -1.f;
  float chi2_3prongs = -1.f;
  float chi2_topology = -1.f;
  ClassDef(RHyperTriton3KF,1)
};

template<class Hyper>
struct SHyperTriton : public Hyper {
  SHyperTriton() : Hyper{} {}
  SHyperTriton(const Hyper& other) : Hyper{other} {}
  virtual ~SHyperTriton() = default;
  float gPt = -999.f;
  float gPhi = -999.f;
  float gPz = -999.f;
  float gCt = -1.f;
  float gT = -1.f;
  bool  gPositive = false;
  bool  gReconstructed = false;
  ClassDef(SHyperTriton,1)
};

using SHyperTriton3KF = SHyperTriton<RHyperTriton3KF>;
using SHyperTriton3O2 = SHyperTriton<RHyperTriton3O2>;

#endif