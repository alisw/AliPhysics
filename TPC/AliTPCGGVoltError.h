#ifndef ALI_TPC_GG_VOLT_ERROR_H
#define ALI_TPC_GG_VOLT_ERROR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCGGVoltError class                                                //
// The class calculates the electric field and space point distortions    //
// due a Gating Grid (GG) Error voltage. It uses the exact calculation    //
// technique based on bessel functions. (original code from STAR)         //
// The class allows "effective Omega Tau" corrections.                    // 
//                                                                        //
// date: 27/04/2010                                                       //
// Authors: Jim Thomas, Stefan Rossegger, Magnus Mager                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"

class AliTPCGGVoltError : public AliTPCCorrection {
public:
  AliTPCGGVoltError();
  virtual ~AliTPCGGVoltError();

  // common setters and getters for ExB
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
    const Double_t wt0=t2*omegaTau;
    fC0=1./(1.+wt0*wt0);
    const Double_t wt1=t1*omegaTau;
    fC1=wt1/(1.+wt1*wt1);
  };

  void SetC0C1(Double_t c0,Double_t c1) {fC0=c0;fC1=c1;} // CAUTION: USE WITH CARE
  Float_t GetC0() const {return fC0;}
  Float_t GetC1() const {return fC1;}

  // setters and getters for GG
  void SetDeltaVGGA(Double_t deltaVGGA) {fDeltaVGGA=deltaVGGA;}
  void SetDeltaVGGC(Double_t deltaVGGC) {fDeltaVGGC=deltaVGGC;}
  Double_t GetDeltaVGGA() const {return fDeltaVGGA;}
  Double_t GetDeltaVGGC() const {return fDeltaVGGC;}

  void InitGGVoltErrorDistortion();

  Float_t GetIntErOverEz(const Float_t x[],const Short_t roc);

  virtual void Print(Option_t* option="") const;

protected:
  virtual void GetCorrection(const Float_t x[],const Short_t roc, Float_t dx[]);
private:

  Float_t fC0; // coefficient C0                 (compare Jim Thomas's notes for definitions)
  Float_t fC1; // coefficient C1                 (compare Jim Thomas's notes for definitions)

  Double_t fDeltaVGGA;            // Missmatch of gating grid voltage on A-side [V]
  Double_t fDeltaVGGC;            // Missmatch of gating grid voltage on C-side [V]
  Double_t fGGVoltErrorER[kNZ][kNR]; // Array to store electric field for GGVoltError calculation

  ClassDef(AliTPCGGVoltError,1);
};

#endif
