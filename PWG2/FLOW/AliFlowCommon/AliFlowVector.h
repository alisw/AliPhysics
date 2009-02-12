/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliFlowVector_H
#define AliFlowVector_H

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
  AliFlowVector(const TVector2 &p, Double_t m);
  virtual ~AliFlowVector();

  AliFlowVector& operator=(const AliFlowVector& aVector);

  void SetMult(Double_t mult)       {this->fMult = mult;} ;
  Double_t GetMult() const          {return this -> fMult;} ;
  
 private:
  Double_t fMult;   //multiplicity
        
  ClassDef(AliFlowVector, 0) 
};
#endif


