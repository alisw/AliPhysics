// -*- C++ -*-

#include <TObject.h>
#include <TVectorD.h>
#include <TString.h>

#include "AliLog.h"

class AliDoubleGaussianBeamProfile : public TObject {  
public:

  static void Eval(Double_t sepX, Double_t sepY, const TVectorD &par, TVectorD &profile);

  static Int_t GetNPar() { return 20; }
  static TString GetParameterName(Int_t i);

protected:
private:
  ClassDef(AliDoubleGaussianBeamProfile, 1);
} ;
