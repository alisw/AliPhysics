#include <Riostream.h>
#include "AliRun.h"
#include "TFluka.h"
#ifndef WIN32
# define mgdraw mgdraw_
#else
# define mgdraw MGDRAW
#endif

extern "C" {
void mgdraw(Int_t& icode, Int_t& mreg)
{
  ((TFluka*) gMC)->SetIcode(icode);
  ((TFluka*) gMC)->SetMreg(mreg);
  cout << endl << " !!! I am in mgdraw - calling gAlice->Stepping()" << endl;
  ((TFluka*) gMC)->FutoTest();
//  gAlice->Stepping();
} // end of mgdraw
} // end of extern "C"

