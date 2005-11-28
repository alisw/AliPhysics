#ifndef ALITRACKRESIDUALSCHI2_H
#define ALITRACKRESIDUALSCHI2_H

//************************************************************************
// AliTrackResidualsChi2: derived class (from AliTrackResiduals) which   *
// implements a minimization of the track residuals based on chi2        *
// approach.                                                             *
//                                                                       *
//************************************************************************

#include "AliAlignObj.h"
#include "AliTrackResiduals.h"

class AliTrackResidualsChi2 : public AliTrackResiduals {

 public:
  AliTrackResidualsChi2():AliTrackResiduals() { }
  AliTrackResidualsChi2(Int_t ntracks, AliAlignObj *alignobj):AliTrackResiduals(ntracks,alignobj) { }
  AliTrackResidualsChi2(const AliTrackResidualsChi2 &res):AliTrackResiduals(res) { }
  AliTrackResidualsChi2& operator= (const AliTrackResidualsChi2& res) { ((AliTrackResiduals *)this)->operator=(res); return *this; }
  virtual ~AliTrackResidualsChi2() { }

  Bool_t Minimize();

  void   Chi2(Int_t & /* npar */, Double_t * /* gin */, Double_t &f, Double_t *par, Int_t /* iflag */);

 protected:

  ClassDef(AliTrackResidualsChi2,1)

};

#endif
