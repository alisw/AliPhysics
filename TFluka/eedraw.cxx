#include <Riostream.h>
#include "TVirtualMCApplication.h"

#include "TFluka.h"
#include "TFlukaCodes.h"

#ifndef WIN32
# define eedraw eedraw_
#else
# define eedraw EEDRAW
#endif
extern "C" {
void eedraw(Int_t& icode)
{
  ((TFluka*) gMC)->SetCaller(kEEDRAW);
  ((TFluka*) gMC)->SetIcode((FlukaProcessCode_t) icode);
} // end of eedraw
} // end of extern "C"

