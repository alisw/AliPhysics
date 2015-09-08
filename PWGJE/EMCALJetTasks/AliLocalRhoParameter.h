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
  Double_t GetLocalValInPhi(Double_t phi, Double_t r, Double_t n) const {
      // calls the old 'get local val', function added for transparency
      // integrates the local rho function only in phi 
      return GetLocalVal(phi, r, n); 
  }
  Double_t GetLocalValInEtaPhi(Double_t phi, Double_t r, Double_t n, Int_t gran = 20) const {
      // numerical approximation for integrating over a circular area in eta phi
      // instead of 'rectangle' around phi - gran determines step size of integral
      if(!fLocalRho) return GetVal();       // fallback value
      Double_t denom(2*r*fLocalRho->GetParameter(0));
      if (denom < 0) return GetVal();
      Double_t sum(0), sumWeight(0), weight(0);
      // make sure gran is an even number and then divide it by two 
      if(gran%2 != 0) gran++;
      gran/=2;
      gran = TMath::Nint(gran);
      for(Int_t i(0); i < gran; i++) {
          // since the function has no eta dependence, it is just weighted by the approximate
          // fraction of the circle that covers a given eta, phi region set by the granularity
          weight = TMath::Sqrt(r*r-((r/(double)gran)*(double)i)*(r/(double)gran)*(double)i);
          sumWeight+=weight;
          sum += TMath::Abs(weight * fLocalRho->Integral(phi + (r/(double)gran)*(double)(i), phi + (r/(double)gran)*(double)(i+1)));
          sum += TMath::Abs(weight * fLocalRho->Integral(phi - (r/(double)gran)*(double)(i), phi - (r/(double)gran)*(double)(i+1)));
      }
      if(sumWeight>0) sum/=sumWeight;
      else return GetVal();
      return (double)gran*n*(sum/denom);
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
