#include <Riostream.h>
#include "TFluka.h"
#ifndef WIN32
# define usdraw usdraw_
#else
# define usdraw USDRAW
#endif
extern "C" {
void usdraw(Int_t& icode, Int_t& mreg, 
            Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
  ((TFluka*) gMC)->SetIcode(icode);
  ((TFluka*) gMC)->SetMreg(mreg);
  ((TFluka*) gMC)->SetXsco(xsco);
  ((TFluka*) gMC)->SetYsco(ysco);
  ((TFluka*) gMC)->SetZsco(zsco);
  cout << endl << " !!! I am in usdraw" << endl;
  ((TFluka*) gMC)->FutoTest();
} // end of usdraw
} // end of extern "C"

