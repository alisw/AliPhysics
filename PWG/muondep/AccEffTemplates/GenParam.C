#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TRandom.h"
#include "AliGenerator.h"
#include "AliGenParam.h"
#endif

AliGenerator* GenParam()
{
  AliGenMC* generator = new AliGenParam(1, VAR_GENPARAM_GENLIB_TYPE,VAR_GENPARAM_GENLIB_PARNAME);

  generator->SetMomentumRange(0,999);
  generator->SetPtRange(0,50.);
  generator->SetYRange(-4.1,-2.4);
  generator->SetPhiRange(0., 360.);
  generator->SetChildThetaRange(0.,180.);
  generator->SetForceDecay(kDiMuon);

  generator->SetCutOnChild(1);
  generator->SetChildPhiRange(0.,360.);
  generator->SetTrackingFlag(1);

  return generator;
}
