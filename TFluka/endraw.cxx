#include <Riostream.h>
#include "TVirtualMCApplication.h"

#ifndef WITH_ROOT
#include "TFluka.h"
#else
#include "TFlukaGeo.h"
#endif

#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#ifndef WIN32
# define endraw endraw_
#else
# define endraw ENDRAW
#endif
extern "C" {
void endraw(Int_t& icode, Int_t& mreg, Double_t& rull, Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
  TFluka* fluka = (TFluka*) gMC;
  fluka->SetCaller(3);
  fluka->SetIcode(icode);
  fluka->SetRull(rull);
  fluka->SetXsco(xsco);
  fluka->SetYsco(ysco);
  fluka->SetZsco(zsco);
  fluka->SetMreg(mreg);
  if (icode == 11) {
    cout << " For icode=" << icode << " Stepping is NOT called" << endl;
    return;
  }
  (TVirtualMCApplication::Instance())->Stepping();
  fluka->SetTrackIsNew(kFALSE);
} // end of endraw
} // end of extern "C"

