#ifndef ALIGAUSCORR_H
#define ALIGAUSCORR_H
/* Copyright(c) 2001-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
// Class used to generate correlated gaussian numbers with mean
// zero and known covariance matrix.
// M.Masera 15.03.2001 9:30 - modified on 26.02.2002 17:40
////////////////////////////////////////////////////////////////////////

#include <TMatrixD.h>
class TArrayD;


class AliGausCorr : public TObject 
{
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



