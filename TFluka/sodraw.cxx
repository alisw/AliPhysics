#include <Riostream.h>
#include "TFluka.h"
#ifndef WIN32
# define sodraw sodraw_
#else
# define sodraw SODRAW
#endif
extern "C" {
void sodraw()
{
  ((TFluka*) gMC)->SetIcode(0);
  cout << endl << " !!! I am in sodraw" << endl;
  ((TFluka*) gMC)->FutoTest();
} // end of sodraw
} // end of extern "C"

