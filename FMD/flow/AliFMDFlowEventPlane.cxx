/** @file 
    @brief Implementation of an EventPlane class */
#include "flow/AliFMDFlowEventPlane.h"
#include "flow/AliFMDFlowUtil.h"
#include <TMath.h>
// #include <cmath>
#ifndef _GNU_SOURCE
extern "C" 
{
  /** Function to caculate @f$ \sin(a), \cos(a)@f$ in one go.  Note,
      that with GCC, this is a built-in, and we only have this to
      serve as a replacement function 
      @ingroup utils */
  inline void sincos(Double_t a, Double_t* sina, Double_t* cosa) {
    *sina = sin(a);
    *cosa = cos(a);
  }
}
#endif

//====================================================================
void 
AliFMDFlowEventPlane::Clear(Option_t*) 
{ 
  fSumSinMPhi = 0;
  fSumCosMPhi = 0;
  fCache      = -1;
}
//____________________________________________________________________
void 
AliFMDFlowEventPlane::Add(Double_t phi, Double_t weight) 
{ 
  Double_t a = NormalizeAngle(fOrder * phi);
  Double_t s, c;
  sincos(a, &s, &c);
  if (TMath::IsNaN(s) || !TMath::Finite(s) || 
      TMath::IsNaN(c) || !TMath::Finite(s)) return;
  fSumSinMPhi += weight * s;
  fSumCosMPhi += weight * c;
}
//____________________________________________________________________
Double_t 
AliFMDFlowEventPlane::Psi() const  
{ 
  if (fCache < 0) fCache = DoPsi(fSumSinMPhi, fSumCosMPhi);
  return fCache;
}
//____________________________________________________________________
Double_t 
AliFMDFlowEventPlane::Psi(Double_t phi, Double_t w) const  
{ 
  Double_t a = NormalizeAngle(fOrder * phi);
  Double_t s, c;
  sincos(a, &s, &c);
  if (TMath::IsNaN(s) || !TMath::Finite(s) || 
      TMath::IsNaN(c) || !TMath::Finite(s)) return Psi();
  Double_t psi = DoPsi(fSumSinMPhi - w * s, fSumCosMPhi - w * c);
  return psi;
}

//____________________________________________________________________
Double_t 
AliFMDFlowEventPlane::DoPsi(Double_t sumsin, Double_t sumcos) const
{
  Double_t psi = 0;
  // Make sure we get an angle everywhere 
  if      (sumcos != 0) psi =  atan2(sumsin, sumcos);
  else if (sumsin == 0) psi =  0;
  else if (sumsin >  0) psi =  M_PI / 2;
  else                  psi = -M_PI / 2;
  psi =  NormalizeAngle(psi);
  psi /= fOrder;
  return psi;
}

//____________________________________________________________________
//
// EOF
//

