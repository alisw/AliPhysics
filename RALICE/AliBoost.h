#ifndef ALIBOOST_H
#define ALIBOOST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <iostream.h>
#include <math.h>
 
#include "TObject.h"

#include "Ali4Vector.h" 

class AliBoost : public TObject
{
 public:
  AliBoost();                             // Default constructor
  ~AliBoost();                            // Default destructor
  void SetBeta(Ali3Vector b);             // Set boost parameters by beta 3-vector
  void Set4Momentum(Ali4Vector& p);       // Set boost parameters by 4-momentum
  Ali3Vector GetBetaVector();             // Provide the beta 3-vector
  Double_t GetBeta();                     // Provide norm of beta 3-vector
  Double_t GetGamma();                    // Provide gamma value
  void Info(TString f="car");             // Print boost parameter info in coord. frame f
  Ali4Vector Boost(Ali4Vector& v);        // Perform Lorentz boost on 4-vector v
  Ali4Vector Inverse(Ali4Vector& v);      // Perform inverse Lorentz boost on 4-vector v
  Double_t GetResultError();              // Provide error on scalar result
 
 protected:
  Ali3Vector fBeta;  // The beta 3-vector
  Double_t fGamma;   // The gamma factor
  Double_t fDgamma;  // Error on the gamma value
  Double_t fDresult; // Error on scalar result
 
 ClassDef(AliBoost,1) // Perform various Lorentz transformations.
};
#endif
