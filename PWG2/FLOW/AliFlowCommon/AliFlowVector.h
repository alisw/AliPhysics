/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIFLOWVECTOR_H
#define ALIFLOWVECTOR_H

#include "TVector2.h"

//********************************************************************
// AliFlowVector:                                                    *
// Class to hold the flow vector and multiplicity for flow analysis. *
// Author: A. Bilandzic (anteb@nikhef.nl)                            *
//********************************************************************

class AliFlowVector: public TVector2 {
 public:
  AliFlowVector();
  AliFlowVector(const AliFlowVector& aVector);
  AliFlowVector(const TVector2 &p, const Double_t m);
  virtual ~AliFlowVector();

  AliFlowVector& operator=(const AliFlowVector& aVector);
  AliFlowVector& operator+=(const AliFlowVector& aVector);

  Bool_t  IsFolder() const {return kTRUE;};

  void SetMult(Double_t const mult) {this->fMult = mult;};
  Double_t GetMult() const {return this->fMult;};
        
 private:
  Double_t fMult;                 // multiplicity = sum of weights = w_1 + w_2 + ... + w_n
   
  ClassDef(AliFlowVector, 1) 
};

inline  AliFlowVector operator+(const AliFlowVector& aVector,const AliFlowVector& bVector) {
  AliFlowVector cVector;
  Double_t x = aVector.X() + bVector.X(); 
  Double_t y = aVector.Y() + bVector.Y(); 
  Double_t mult = aVector.GetMult() + bVector.GetMult();
  cVector.Set(x,y);
  cVector.SetMult(mult);
  
  return cVector;
}

#endif
