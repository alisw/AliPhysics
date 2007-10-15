#include <Riostream.h>
#include "TVirtualMCApplication.h"

#include "TFluka.h"
#include "TFlukaCodes.h"

#ifndef WIN32
# define sodraw sodraw_
#else
# define sodraw SODRAW
#endif
extern "C" {
void sodraw()
{
  ((TFluka*) gMC)->SetCaller(kSODRAW);
  ((TFluka*) gMC)->SetIcode((FlukaProcessCode_t)0);
  (TVirtualMCApplication::Instance())->Stepping();
} // end of sodraw
} // end of extern "C"

