#ifndef ALITRACKRESIDUALSFAST_H
#define ALITRACKRESIDUALSFAST_H

//************************************************************************
// AliTrackResidualsFast: derived class (from AliTrackResiduals) which   *
// implements a simple linear minimization of the track residuals chi2.  *
// The minimization relies on the fact that the alignment parameters     *
// (angles and translations) are small.                                  *
//                                                                       *
//                                                                       *
//************************************************************************

#include "TMatrixDSym.h"
#include "TMatrixD.h"

#include "AliAlignObj.h"
#include "AliTrackResiduals.h"

class AliTrackResidualsFast : public AliTrackResiduals {

 public:
  AliTrackResidualsFast();
  AliTrackResidualsFast(Int_t ntracks);
  AliTrackResidualsFast(const AliTrackResidualsFast &res);
  AliTrackResidualsFast& operator= (const AliTrackResidualsFast& res);
  virtual ~AliTrackResidualsFast() { }

  Bool_t Minimize();

 protected:

  void   AddPoints(AliTrackPoint &p, AliTrackPoint &pprime);
  Bool_t Update();

  Double_t fSum[27]; // Sums used during the chi2 minimization
  Double_t fSumR;    // Sum of r squared

  ClassDef(AliTrackResidualsFast,1)

};

#endif
