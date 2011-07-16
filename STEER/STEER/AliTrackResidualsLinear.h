#ifndef ALITRACKRESIDUALSLINEAR_H
#define ALITRACKRESIDUALSLINEAR_H

//************************************************************************
// AliTrackResidualsLinear: derived class (from AliTrackResiduals) which   *
// implements a simple linear minimization of the track residuals chi2.  *
// The minimization relies on the fact that the alignment parameters     *
// (angles and translations) are small.                                  *
//                                                                       *
//                                                                       *
//************************************************************************

#include "AliAlignObj.h"
#include "AliTrackResiduals.h"
class TLinearFitter;

class AliTrackResidualsLinear : public AliTrackResiduals {

 public:
  AliTrackResidualsLinear();
  AliTrackResidualsLinear(Int_t ntracks);
  AliTrackResidualsLinear(const AliTrackResidualsLinear &res);
  AliTrackResidualsLinear& operator= (const AliTrackResidualsLinear& res);
  virtual ~AliTrackResidualsLinear();
  Bool_t Minimize();
  void   SetRobust(Float_t fraction){fFraction=fraction;}
  const Double_t * GetParameters() const { return fParams;}
  const Double_t * GetCovariance() const { return fCovar;}
 protected:
  Bool_t Update();
  void   AddPoints(AliTrackPoint &p, AliTrackPoint &pprime);
  TLinearFitter *fFitter;           // !Linear fitter
  Float_t fFraction;                // fraction of points for robust fitter if less 0 - no robust fitter invoked
  Double_t fParams[6];               // resulting parameters 
  Double_t fCovar[36];               // covariance matrix 
  Double_t fChi2Orig;               // original chi2 
  ClassDef(AliTrackResidualsLinear,1)

};

#endif
