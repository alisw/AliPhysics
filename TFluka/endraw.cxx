#include <Riostream.h>
#include "AliRun.h"
#include "TFluka.h"
#ifndef WIN32
# define endraw endraw_
#else
# define endraw ENDRAW
#endif
extern "C" {
void endraw(Int_t& icode, Int_t& mreg, Double_t& rull, Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
  ((TFluka*) gMC)->SetIcode(icode);
  ((TFluka*) gMC)->SetMreg(mreg);
  ((TFluka*) gMC)->SetRull(rull);
  ((TFluka*) gMC)->SetXsco(xsco);
  ((TFluka*) gMC)->SetYsco(ysco);
  ((TFluka*) gMC)->SetZsco(zsco);
  cout << endl << " !!! I am in endraw - calling gAlice->Stepping()" << endl;
  ((TFluka*) gMC)->FutoTest();
//  gAlice->Stepping();
} // end of endraw
} // end of extern "C"

