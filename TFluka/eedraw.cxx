#include <Riostream.h>
#include "TVirtualMCApplication.h"

#include "TFluka.h"

#ifndef WIN32
# define eedraw eedraw_
#else
# define eedraw EEDRAW
#endif
extern "C" {
void eedraw(Int_t& icode)
{
  ((TFluka*) gMC)->SetCaller(2);
  ((TFluka*) gMC)->SetIcode(icode);
} // end of eedraw
} // end of extern "C"

