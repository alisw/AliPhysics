#ifndef ALITPCEXBEFFECTIVE_H
#define ALITPCEXBEFFECTIVE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCExBEffective class                                                   //
// date: 02/05/2010                                                       //
// Authors: Maarian Ivanov, Jim Thomas, Magnus Mager, Stefan Rossegger                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"

class AliTPCExBEffective : public AliTPCCorrection {
public:
  AliTPCExBEffective();
  virtual ~AliTPCExBEffective();
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
  Double_t GetSum(const TMatrixD& pol, const TMatrixD&coef, Double_t r, Double_t drift, Double_t phi, Int_t coord=0) const;
  void SetPolynoms(const TMatrixD *polA, const TMatrixD *polC);
  void SetCoeficients(const TMatrixD *valA,const TMatrixD *valC);
  void Print(const Option_t* option) const;

public:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

private:
  Double_t fC0;                // coefficient C0 (compare Jim Thomas's notes for definitions)
  Double_t fC1;                // coefficient C1 (compare Jim Thomas's notes for definitions)
  TMatrixD *fPolynomA;         // correction polynoms A
  TMatrixD *fPolynomC;         // correction polynoms C
  TMatrixD *fPolynomValA;      // correction polynoms coefficient A
  TMatrixD *fPolynomValC;      //  correction polynoms coefficient C

  AliTPCExBEffective(const AliTPCExBEffective&);
  AliTPCExBEffective &operator=(const AliTPCExBEffective&);
  ClassDef(AliTPCExBEffective,1);
};

#endif
