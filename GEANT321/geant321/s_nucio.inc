*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:44  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_S_NUCIO_INC
#define CERNLIB_GEANT321_S_NUCIO_INC
*
*
* s_nucio.inc
*
      COMMON/NUCIN /TECM,AMASS(18),NPG,KGENEV
      COMMON/NUCOUT/PCM(5,18),WGT
#if !defined(CERNLIB_SINGLE)
#include "geant321/nuciod.inc"
#endif
C

#endif
