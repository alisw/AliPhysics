#ifndef ALIGAUSCORR_H
#define ALIGAUSCORR_H
/* Copyright(c) 2001-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TMatrixD.h>
class TArrayD;


class AliGausCorr : public TObject {
////////////////////////////////////////////////////////////////////////
// Class used to generate correlated gaussian numbers with mean
// zero and known covariance matrix.
// Adapted from the Fortran code in Cernlib V122 (corset, corgen)
// F. James, Monte Carlo theory and practice, 
// Rep. Prog. Phys. 43 (1980) 1145-1189. 
// M.Masera 14.03.2001 19:30 - last mod. 26.02.2002 17:45
////////////////////////////////////////////////////////////////////////
 public:
  //
  AliGausCorr();
  AliGausCorr(const TMatrixD & cov, Int_t size);
  AliGausCorr(const AliGausCorr & tgcorr);
  virtual ~AliGausCorr();
  void GetGaussN(TArrayD &vec) const;
  TMatrixD GetSqrtMatrix() const { return *fCv;}
  void PrintCv() const;
  AliGausCorr & operator=(const AliGausCorr & tgcorr);
  //
 private:
  //
  Int_t fSize;   // number of correlated gaussian random numbers
  TMatrixD *fCv; // 'square root' of the covariance matrix

  ClassDef(AliGausCorr,1)  
};


#endif



