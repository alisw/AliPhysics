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
  AliFlowVector(const TVector2 &p, Double_t m);                      // constructor: Add a weight to the TVector
  AliFlowVector(Double_t *y, Double_t m=1);                          // constructor: Analogue to TVector2(y) with multiplicity
  AliFlowVector(Double_t x, Double_t y, Double_t m=1);   // constructor: Sets the components individually
  virtual ~AliFlowVector();

  void SetMagPhi(Double_t size, Double_t angle, Double_t mult=1);          // Set vector and weighted multiplicity

  AliFlowVector& operator=(const AliFlowVector& aVector);                        // Assign to self
  AliFlowVector& operator+=(const AliFlowVector& aVector);                       // Add to self
  AliFlowVector& operator-=(const AliFlowVector& aVector);                       // Subtract from self
  AliFlowVector& operator*=(double w);                                     // Multiply by a weight
  AliFlowVector& operator/=(double w){ (*this)*=(1.0/w); return *this;};      // Divide by a weight
  const AliFlowVector operator+(const AliFlowVector&a) const { AliFlowVector v(*this); return v+=a; };   // Add and return by value
  const AliFlowVector operator-(const AliFlowVector&a) const { AliFlowVector v(*this); return v-=a; };   // Subtract and return by value
  const AliFlowVector operator*(double w) const { AliFlowVector v(*this); return v*=w; };          // Scale and return by value
  const AliFlowVector operator/(double w) const { AliFlowVector v(*this); return v/=w; };          // Scale and return by value

  Bool_t  IsFolder() const {return kTRUE;};

  void SetMult(Double_t mult) {fMult = mult;};           // Set sum of weights
  Double_t GetMult() const {return fMult;};                    // Get sum of weights
        
 private:
  Double_t fMult;                 // multiplicity = sum of weights = w_1 + w_2 + ... + w_n
   
  ClassDef(AliFlowVector, 1) 
};


/* Old, less efficient code
inline  AliFlowVector operator+(const AliFlowVector& aVector,const AliFlowVector& bVector) {
  AliFlowVector cVector;
  Double_t x = aVector.X() + bVector.X(); 
  Double_t y = aVector.Y() + bVector.Y(); 
  Double_t mult = aVector.GetMult() + bVector.GetMult();
  cVector.Set(x,y);
  cVector.SetMult(mult);
  
  return cVector;
}

inline  AliFlowVector operator-(const AliFlowVector& aVector,const AliFlowVector& bVector) 
{
   // Difference between two vectors
  AliFlowVector cVector;
  Double_t x = aVector.X() - bVector.X(); 
  Double_t y = aVector.Y() - bVector.Y(); 
  Double_t mult = aVector.GetMult() - bVector.GetMult();
  cVector.Set(x,y);
  cVector.SetMult(mult);
  
  return cVector;
}
*/

#endif
