/** @file 
    @brief Implementation of a Harmonic class */
#include "flow/AliFMDFlowHarmonic.h"
#include "flow/AliFMDFlowUtil.h"
// #include <cmath>

//====================================================================
void 
AliFMDFlowHarmonic::Add(Double_t phi, Double_t psi) 
{ 
  Double_t a       = NormalizeAngle(fOrder * (phi - psi));
  Double_t contrib = cos(a);
  AliFMDFlowStat::Add(contrib);
}
//____________________________________________________________________
Double_t 
AliFMDFlowHarmonic::Value(Double_t r, Double_t er2, Double_t& e2) const 
{
  // The corrected value is given by 
  //
  //          v_n^obs
  //    v_n = -------
  //             R
  // 
  // where 
  // 
  //              1
  //    v_n^obs = - \sum_i(cos(n(\phi_i - \Psi)))
  //              N 
  // 
  // and R is the resolution 
  // 
  // The error on the corrected value is given by 
  //
  //                dv_n                    dv_n
  //    d^2v_n = (--------)^2 d^2v_n^obs + (----)^2 d^2R
  //              dv_n^obs                   dR
  //
  //             d^2v_n^obs R^2 + d^2R v_n^obs^2
  //           = -------------------------------
  //                         R^4 
  // 
  Double_t a = fAverage;
  Double_t v = a / r;
  Double_t s = fSqVar / fN;
  e2       = (s * r * r + er2 * a * a) / pow(r, 4);
  return v;
}
//____________________________________________________________________
//
// EOF
//
