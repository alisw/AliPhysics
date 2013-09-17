#ifndef ALILOCALRHOPARAMETER_H
#define ALILOCALRHOPARAMETER_H

// $Id$

#include <TMath.h>
#include <TF1.h>
#include <AliRhoParameter.h>

class AliLocalRhoParameter : public AliRhoParameter {
 public: 
  AliLocalRhoParameter();
  AliLocalRhoParameter(const char* name, Double_t val);
  void     SetLocalRho(TF1* f)             { fLocalRho = f;        }
  TF1*     GetLocalRho() const             { return fLocalRho;     }
  Double_t GetLocalVal(Double_t phi, Double_t r, Double_t n) const {
    if(!fLocalRho) return GetVal();
    Double_t denom(2*r*fLocalRho->GetParameter(0));
    return  (denom <= 0.) ? GetVal() : n*(fLocalRho->Integral(phi-r, phi+r)/denom); 
  }
  Double_t GetLocalVal(Double_t phi, Double_t r) const {
    return GetLocalVal(phi, r, GetVal());
  }
  Double_t GetLocalUncertainty(Double_t phi, Double_t r, Double_t n) const {
      if(!fLocalRho) return 999.;
      Double_t intError(fLocalRho->IntegralError(phi-r,phi+r));
      Double_t absConst(TMath::Abs(n/(2*r*fLocalRho->GetParameter(0))));
      return  intError*absConst;        // absolute error on local rho
  }
  Double_t GetLocalUncertainty(Double_t phi, Double_t r) const {
      return GetLocalUncertainty(phi, r, GetVal());
  }
 private:
  TF1*     fLocalRho;      // ! rho as function of phi

  AliLocalRhoParameter(const AliLocalRhoParameter&);             // not implemented
  AliLocalRhoParameter& operator=(const AliLocalRhoParameter&);  // not implemented

  ClassDef(AliLocalRhoParameter, 1); // Rho parameter for local (flow) variations
};
#endif
