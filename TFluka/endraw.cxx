#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TFluka.h"
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
  ((TFluka*) gMC)->SetCaller(3);
  ((TFluka*) gMC)->SetIcode(icode);
  ((TFluka*) gMC)->SetRull(rull);
  ((TFluka*) gMC)->SetXsco(xsco);
  ((TFluka*) gMC)->SetYsco(ysco);
  ((TFluka*) gMC)->SetZsco(zsco);
  ((TFluka*) gMC)->SetMreg(mreg);
  (TVirtualMCApplication::Instance())->Stepping();
} // end of endraw
} // end of extern "C"

