#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TFluka.h"
#ifndef WIN32
# define bxdraw bxdraw_
#else
# define bxdraw BXDRAW
#endif
extern "C" {
void bxdraw(Int_t& icode, Int_t& mreg, Int_t& newreg,
            Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
  ((TFluka*) gMC)->SetIcode(icode);
  ((TFluka*) gMC)->SetMreg(mreg);
  ((TFluka*) gMC)->SetNewreg(newreg);
  ((TFluka*) gMC)->SetXsco(xsco);
  ((TFluka*) gMC)->SetYsco(ysco);
  ((TFluka*) gMC)->SetZsco(zsco);
//  cout << endl << " !!! I am in bxdraw - calling Stepping()" << mreg << endl;
  ((TFluka*) gMC)->FutoTest();
  (TVirtualMCApplication::Instance())->Stepping();
} // end of bxdraw
} // end of extern "C"

