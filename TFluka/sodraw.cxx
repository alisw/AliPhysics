#include <Riostream.h>
#include "TVirtualMCApplication.h"

#include "TFluka.h"

#ifndef WIN32
# define sodraw sodraw_
#else
# define sodraw SODRAW
#endif
extern "C" {
void sodraw()
{
  ((TFluka*) gMC)->SetCaller(5);
  ((TFluka*) gMC)->SetIcode(0);
} // end of sodraw
} // end of extern "C"

