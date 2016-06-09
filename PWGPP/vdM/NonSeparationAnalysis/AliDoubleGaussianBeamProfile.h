// -*- C++ -*-

#ifndef _ALI_DOUBLE_GAUSSIAN_BEAM_PROFILE_H_
#define _ALI_DOUBLE_GAUSSIAN_BEAM_PROFILE_H_

#include <TObject.h>
#include <TVectorD.h>

#include "AliLog.h"

class AliDoubleGaussianBeamProfile : public TObject {  
public:

  static Bool_t Eval(Double_t sepX, Double_t sepY, const TVectorD &par, TVectorD &profile);

  static Int_t GetNPar() { return 20; }
  static const char* GetParName(Int_t i);

  static Double_t EvalProfile0(Double_t *x, Double_t *par); // returns profile(0)

protected:
private:
  ClassDef(AliDoubleGaussianBeamProfile, 1);
} ;

#endif // _ALI_DOUBLE_GAUSSIAN_BEAM_PROFILE_H_
