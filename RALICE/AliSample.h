#ifndef ALISAMPLE_H
#define ALISAMPLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// Class AliSample
// Perform statistics on various multi-dimensional data samples
// A data sample can be filled using the "Enter" and/or "Remove" functions,
// whereas the "Reset" function resets the complete sample to 'empty'.
// The info which can be extracted from a certain data sample are the
// sum, mean, variance, sigma, covariance and correlation.
// The "Info" function provides all statistics data for a certain sample.
// The variables for which these stat. parameters have to be calculated
// are indicated by the index of the variable which is passed as an
// argument to the various member functions.
// The index convention for a data point (x,y) is : x=1  y=2
//
// Example :
// ---------
// For an AliSample s a data point (x,y) can be entered as s.Enter(x,y) and
// the mean_x can be obtained as s.GetMean(1) whereas the mean_y is obtained
// via s.GetMean(2).
// The correlation between x and y is available via s.GetCor(1,2).
// The x-statistics are obtained via s.Info(1), y-statistics via s.Info(2),
// and the covariance and correlation between x and y via s.Info(1,2).
// All statistics of a sample are obtained via s.Info().
//
//--- NvE 30-mar-1996 CERN Geneva
///////////////////////////////////////////////////////////////////////////
 
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
