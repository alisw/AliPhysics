#include <Riostream.h>
#include "AliRun.h"
#include "TFluka.h"
#ifndef WIN32
# define eedraw eedraw_
#else
# define eedraw EEDRAW
#endif
extern "C" {
void eedraw(Int_t& icode)
{
  ((TFluka*) gMC)->SetIcode(icode);
  cout << endl << " !!! I am in eedraw - calling gAlice->Stepping()" << endl;
  ((TFluka*) gMC)->FutoTest();
//  gAlice->Stepping();
} // end of eedraw
} // end of extern "C"

