#include <Riostream.h>
#include "TVirtualMCApplication.h"

#ifndef WITH_ROOT
#include "TFluka.h"
#else
#include "TFlukaGeo.h"
#endif

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

