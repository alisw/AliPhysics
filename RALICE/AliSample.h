#ifndef ALISAMPLE_H
#define ALISAMPLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>

#include "Rtypes.h"
#include "TArrayF.h"
 
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
  Int_t GetDimension() const;                   // Provide dimension of the sample
  Int_t GetN() const;                           // Provide the number of entries
  Float_t GetSum(Int_t i) const;                // Provide sum for i-th variable
  Float_t GetMean(Int_t i) const;               // Provide mean for i-th variable
  Float_t GetVar(Int_t i) const;                // Provide variance for i-th variable
  Float_t GetSigma(Int_t i) const;              // Standard deviation for i-th variable
  Float_t GetCov(Int_t i, Int_t j) const;       // Covariance for i-th and j-th variable
  Float_t GetCor(Int_t i, Int_t j) const;       // Correlation for i-th and j-th variable
  Float_t GetMedian(Int_t i);                   // Provide median for i-th variable
  Float_t GetSpread(Int_t i);                   // Provide spread w.r.t. the median for i-th variable
  Float_t GetMinimum(Int_t i) const;            // Provide the minimum value for i-th variable
  Float_t GetMaximum(Int_t i) const;            // Provide the maximum value for i-th variable
  void Data();                                  // Stat. info for the complete sample
  void Data(Int_t i);                           // Stat. info for the i-th variable
  void Data(Int_t i, Int_t j) const;            // Stat. info for i-th and j-th variable
  void SetStoreMode(Int_t mode=1);              // Set mode for storage of all entered data
  Int_t GetStoreMode() const;                   // Provide storage mode of all entered data
 
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
  Float_t fMin[fMaxdim];           // Minimum value for each variable
  Float_t fMax[fMaxdim];           // Maximum value for each variable
  Int_t fRemove;                   // Flag to indicate that some entry has been removed
  Int_t fStore;                    // Flag to denote storage of all entered data 
  TArrayF* fX;                     // Storage array for the 1st variable (e.g. X)
  TArrayF* fY;                     // Storage array for the 2nd variable (e.g. Y)
  TArrayF* fZ;                     // Storage array for the 3rd variable (e.g. Z)
  TArrayF* fArr;                   // Temp. storage array for ordering

 ClassDef(AliSample,0) // Statistics tools for various multi-dimensional data samples.
};
#endif
