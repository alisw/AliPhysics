#include <Riostream.h>
#include "AliRun.h"
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
  cout << endl << " !!! I am in sodraw - calling gAlice->Stepping()" << endl;
  ((TFluka*) gMC)->FutoTest();
//  gAlice->Stepping();
} // end of sodraw
} // end of extern "C"

