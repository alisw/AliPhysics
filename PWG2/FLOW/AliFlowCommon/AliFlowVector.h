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
  AliFlowVector(const TVector2 &p, const Double_t m, const Double_t sumPow2w, const Double_t sumPow3w, const Double_t sumPow4w, const Double_t sumPow5w, const Double_t sumPow6w, const Double_t sumPow7w, const Double_t sumPow8w);
  virtual ~AliFlowVector();

  AliFlowVector& operator=(const AliFlowVector& aVector);

  Bool_t  IsFolder() const {return kTRUE;};

  void SetMult(Double_t const mult) {this->fMult = mult;};
  Double_t GetMult() const {return this->fMult;};
  
  void SetSumOfWeightsToPower2(Double_t const p2w) {this->fSumOfWeightsToPower2 = p2w;};
  Double_t GetSumOfWeightsToPower2() const {return this->fSumOfWeightsToPower2;};

  void SetSumOfWeightsToPower3(Double_t const p3w) {this->fSumOfWeightsToPower3 = p3w;};
  Double_t GetSumOfWeightsToPower3() const {return this->fSumOfWeightsToPower3;};
  
  void SetSumOfWeightsToPower4(Double_t const p4w) {this->fSumOfWeightsToPower4 = p4w;};
  Double_t GetSumOfWeightsToPower4() const {return this->fSumOfWeightsToPower4;};
  
  void SetSumOfWeightsToPower5(Double_t const p5w) {this->fSumOfWeightsToPower5 = p5w;};
  Double_t GetSumOfWeightsToPower5() const {return this->fSumOfWeightsToPower5;};
  
  void SetSumOfWeightsToPower6(Double_t const p6w) {this->fSumOfWeightsToPower6 = p6w;};
  Double_t GetSumOfWeightsToPower6() const {return this->fSumOfWeightsToPower6;};
  
  void SetSumOfWeightsToPower7(Double_t const p7w) {this->fSumOfWeightsToPower7 = p7w;};
  Double_t GetSumOfWeightsToPower7() const {return this->fSumOfWeightsToPower7;};
  
  void SetSumOfWeightsToPower8(Double_t const p8w) {this->fSumOfWeightsToPower8 = p8w;};
  Double_t GetSumOfWeightsToPower8() const {return this->fSumOfWeightsToPower8;};
    
 private:
  Double_t fMult;                 // multiplicity = sum of weights = w_1 + w_2 + ... + w_n
  Double_t fSumOfWeightsToPower2; // pow(w_1,2) + pow(w_2,2) + ... + pow(w_n,2)
  Double_t fSumOfWeightsToPower3; // pow(w_1,3) + pow(w_2,3) + ... + pow(w_n,4)
  Double_t fSumOfWeightsToPower4; // pow(w_1,4) + pow(w_2,4) + ... + pow(w_n,4)
  Double_t fSumOfWeightsToPower5; // pow(w_1,5) + pow(w_2,5) + ... + pow(w_n,5)
  Double_t fSumOfWeightsToPower6; // pow(w_1,6) + pow(w_2,6) + ... + pow(w_n,6)
  Double_t fSumOfWeightsToPower7; // pow(w_1,7) + pow(w_2,7) + ... + pow(w_n,7)
  Double_t fSumOfWeightsToPower8; // pow(w_1,8) + pow(w_2,8) + ... + pow(w_n,8)
 
  ClassDef(AliFlowVector, 1) 
};
#endif


