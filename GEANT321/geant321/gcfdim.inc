*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:36  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCFDIM_INC
#define CERNLIB_GEANT321_GCFDIM_INC
*
*
* gcfdim.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCFDIM/NPROC,NMPTOT,TSEQTO,TOTMBY,TSEQ,
     +TLAT,TNET,X1MIN,X1MAX,Y1MIN,Y1MAX,Z1MIN,Z1MAX,
     +XCUT,YCUT,ZCUT
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtfdim.inc"
#endif
C
#endif

#endif
