#ifndef ALISAMPLE_H
#define ALISAMPLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <math.h>
#include <iostream.h>

#include "Rtypes.h"
 
class AliSample
{
 public:
  AliSample();                                  // Default constructor
  virtual ~AliSample();                         // Default destructor
  void Reset();                                 // Reset complete statistics
  void Enter(Float_t x);                        // Enter value for 1-dim. sample
  void Remove(Float_t x);                       // Remove value from 1-dim. sample
  void Enter(Float_t x, Float_t y);             // Enter value for 2-dim. sample
  void Remove(Float_t x, Float_t y);            // Remove value from 2-dim. sample
  void Enter(Float_t x, Float_t y, Float_t z);  // Enter value for 3-dim. sample
  void Remove(Float_t x, Float_t y, Float_t z); // Remove value from 3-dim. sample
  Int_t GetDimension();                         // Provide dimension of the sample
  Int_t GetN();                                 // Provide the number of entries
  Float_t GetSum(Int_t i);                      // Provide sum for i-th variable
  Float_t GetMean(Int_t i);                     // Provide mean for i-th variable
  Float_t GetVar(Int_t i);                      // Provide variance for i-th variable
  Float_t GetSigma(Int_t i);                    // Standard deviation for i-th variable
  Float_t GetCov(Int_t i, Int_t j);             // Covariance for i-th and j-th variable
  Float_t GetCor(Int_t i, Int_t j);             // Correlation for i-th and j-th variable
  void Info();                                  // Stat. info for the complete sample
  void Info(Int_t i);                           // Stat. info for the i-th variable
  void Info(Int_t i, Int_t j);                  // Stat. info for i-th and j-th variable
 
 private:
  Int_t fDim;                      // Dimension of the sample
  Int_t fN;                        // Number of entries of the sample
  enum  {fMaxdim=3};               // Maximum supported dimension
  char  fNames[fMaxdim];           // Variable names i.e. X,Y,Z
  Float_t fSum[fMaxdim];           // Total sum for each variable
  Float_t fSum2[fMaxdim][fMaxdim]; // Total sum**2 for each variable
  void  Compute();                 // Compute the various quantities
  Float_t fMean[fMaxdim];          // Mean for each variable
  Float_t fVar[fMaxdim];           // Variation for each variable
  Float_t fSigma[fMaxdim];         // Standard deviation for each variable
  Float_t fCov[fMaxdim][fMaxdim];  // Covariances of pairs of variables
  Float_t fCor[fMaxdim][fMaxdim];  // Correlations of pairs of variables
};
#endif
