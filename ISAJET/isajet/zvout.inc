#ifndef CERNLIB_ISAJET_ZVOUT_INC
#define CERNLIB_ISAJET_ZVOUT_INC
*
*
* zvout.inc
*
      COMMON/ZVOUT/ZVOUT(512)
      SAVE /ZVOUT/
#if defined(CERNLIB_LEVEL2)
#include "isajet/l2zout.inc"
* Ignoring t=pass
#endif
      REAL      ZVOUT

#endif
