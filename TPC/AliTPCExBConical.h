#ifndef ALITPCEXBCONICAL_H
#define ALITPCEXBCONICAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCExBConical class                                                   //
// date: 02/05/2010                                                       //
// Authors: Maarian Ivanov, Jim Thomas, Magnus Mager, Stefan Rossegger                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"

class AliTPCExBConical : public AliTPCCorrection {
public:
  AliTPCExBConical();
  virtual ~AliTPCExBConical();

  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);


  // common setters and getters for ExB
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
    fT1=t1; fT2=t2;
    const Float_t wt1=t1*omegaTau;
    fC1=wt1/(1.+wt1*wt1);
    const Float_t wt2=t2*omegaTau;
    fC2=wt2*wt2/(1.+wt2*wt2);
  };
  void SetC1C2(Float_t c1,Float_t c2) {fC1=c1;fC2=c2;} // CAUTION: USE WITH CARE
  Float_t GetC1() const {return fC1;}
  Float_t GetC2() const {return fC2;}

  // setters and getters for conical
  void SetConicalA(Float_t conicalA[3]);
  void SetConicalC(Float_t conicalC[3]);
  void SetConicalFactor(Float_t factor) { fConicalFactor=factor;}

  Float_t GetConicalA(Int_t i) const {return fConicalA[i];}
  Float_t GetConicalC(Int_t i) const {return fConicalC[i];}
  Float_t GetConicalFactor() const {return fConicalFactor;}

  virtual void Print(Option_t* option="") const;

protected:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

private:
  Float_t fC1; // coefficient C1                 (compare Jim Thomas's notes for definitions)
  Float_t fC2; // coefficient C2                 (compare Jim Thomas's notes for definitions)
  Float_t  fConicalFactor;                  // empirical factor - transform conical angle to delta
  Float_t  fConicalA[3];               // Conical shape parameterization A side
  Float_t  fConicalC[3];               // Conical shape parameterization C side

  ClassDef(AliTPCExBConical,1);
};

#endif
