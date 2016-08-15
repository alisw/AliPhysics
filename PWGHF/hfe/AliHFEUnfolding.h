// This class unfolds using simple matrix inversion without Regularization
// martin.andreas.volkl@cern.ch

#ifndef ALIHFEUNFOLDING_H
#define ALIHFEUNFOLDING_H

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TVectorD.h"
#include <iostream>

class AliHFEUnfolding{
private:
  TH1D * fDataHistogram;  // The input Histogram
  TH1D * fUnfoldedHistogram; // The unfolded output
  TH2D * fResponseMatrixHistogram; // ResponseMatrix(true, reconstructed)
  Int_t fLowestBin; // First bin considered in the unfolding, same for input and output
  Int_t fHighestBin; // Last bin considered in the unfolding, same for input and output
  
  TMatrixD * fResponseMatrix; // Response Matrix as TMatrix
  TMatrixD * fCovarianceMatrix; // Covariance Matrix
  TVectorD * fDataVector; // Data Histogram as TVector
  
  
  Bool_t fDataIsSet; // Has the Data Histogram been read in
  Bool_t fResponseIsSet; // Has the Response Matrix been readin
  
  
public:
  AliHFEUnfolding();
  virtual ~AliHFEUnfolding();
  
  void SetData(TH1D * DataHistogram);
  void SetData(TH1D * DataHistogram, Int_t LowestBin, Int_t HighestBin);
  TH1D * GetData(void) const {return fDataHistogram;}
  void SetResponseMatrix(TH2D * ResponseMatrix);
  TH2D * GetResponseMatrix(void) const {return fResponseMatrixHistogram;}
  TMatrixD * GetCovarianceMatrix(void) const {return fCovarianceMatrix;}
  
  TH1D * Unfold(void);
  
  ClassDef(AliHFEUnfolding, 1);
};


#endif
