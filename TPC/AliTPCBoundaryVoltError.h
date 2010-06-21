#ifndef ALITPCBOUNDARYVOLTERROR_H
#define ALITPCBOUNDARYVOLTERROR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCBoundaryVoltError class                                          //
// date: 01/06/2010                                                       //
// Authors: Jim Thomas, Stefan Rossegger                                  //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"


class AliTPCBoundaryVoltError : public AliTPCCorrection {
public:
  AliTPCBoundaryVoltError();
  virtual ~AliTPCBoundaryVoltError();

  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);


  // common setters and getters for tangled ExB effect
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
    fT1=t1; fT2=t2;
    const Double_t wt0=t2*omegaTau;     fC0=1./(1.+wt0*wt0);
    const Double_t wt1=t1*omegaTau;     fC1=wt1/(1.+wt1*wt1);
  };
  void SetC0C1(Float_t c0,Float_t c1) {fC0=c0;fC1=c1;} // CAUTION: USE WITH CARE
  Float_t GetC0() const {return fC0;}
  Float_t GetC1() const {return fC1;}

  // setters and getters for conical
  void SetBoundariesA(Float_t boundariesA[8]);
  void SetBoundariesC(Float_t boundariesC[6]); // CE settings from the A side

  Float_t GetBoundariesA(Int_t i) const {return fBoundariesA[i];}
  Float_t GetBoundariesC(Int_t i) const {return fBoundariesC[i];}

  void InitBoundaryVoltErrorDistortion();

  virtual void Print(const Option_t* option="") const;

protected:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

private:
  Float_t fC0; // coefficient C0                 (compare Jim Thomas's notes for definitions)
  Float_t fC1; // coefficient C1                 (compare Jim Thomas's notes for definitions)
  Float_t  fBoundariesA[8];            // Boundary values on the A side (see Setter function)
  Float_t  fBoundariesC[8];            // Boundary values on the C side (see Setter function)

  Bool_t fInitLookUp;                  // flag to check it the Look Up table was created

  Double_t fLookUpErOverEz[kNZ][kNR];  // Array to store electric field integral (int Er/Ez)

  // basic numbers for the poisson relaxation //can be set individually in each class
  enum {kRows   =257}; // grid size in r direction used in the poisson relaxation // ( 2**n + 1 ) eg. 65, 129, 257 etc.
  enum {kColumns=257}; // grid size in r direction used in the poisson relaxation // ( 2**m + 1 ) eg. 65, 129, 257 etc.
  enum {kIterations=100}; // Number of iterations within the poisson relaxation 

  ClassDef(AliTPCBoundaryVoltError,0); 
};

#endif
