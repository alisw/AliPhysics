#ifndef ALITPCEXBEFFECTIVESECTOR_H
#define ALITPCEXBEFFECTIVESECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCExBEffectiveSector class                                                   //
// date: 02/05/2010                                                       //
// Authors: Maarian Ivanov, Jim Thomas, Magnus Mager, Stefan Rossegger                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"
class TH3F;
class THnSparse;
class THnBase;

class AliTPCExBEffectiveSector : public AliTPCCorrection {
public:
  AliTPCExBEffectiveSector();
  virtual ~AliTPCExBEffectiveSector();
  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);
  // common setters and getters for ExB
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
    fT1=t1; fT2=t2;
    const Float_t wt1=t1*omegaTau;    fC1=wt1/(1.+wt1*wt1);
    const Float_t wt2=t2*omegaTau;    fC0=1/(1.+wt2*wt2);
  };
  Float_t GetC1() const {return fC1;}
  Float_t GetC0() const {return fC0;}
  void Print(const Option_t* option) const;
public:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);
public:
  Double_t fC0;                // coefficient C0 (compare Jim Thomas's notes for definitions)
  Double_t fC1;                // coefficient C1 (compare Jim Thomas's notes for definitions)
  TH3F *fCorrectionR;        // radial correction
  TH3F *fCorrectionRPhi;     // r-phi correction
  TH3F *fCorrectionZ;        // z correction
private:
  AliTPCExBEffectiveSector(const AliTPCExBEffectiveSector&);
  AliTPCExBEffectiveSector &operator=(const AliTPCExBEffectiveSector&);
  ClassDef(AliTPCExBEffectiveSector,2);
};

#endif
