*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:32  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCPARM_INC
#define CERNLIB_GEANT321_GCPARM_INC
*
*
* gcparm.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCPARM/IPARAM,PCUTGA,PCUTEL,PCUTNE,PCUTHA,PCUTMU
     +             ,NSPARA,MPSTAK,NPGENE
      REAL PACUTS(5)
      EQUIVALENCE (PACUTS(1),PCUTGA)
#if defined(CERNLIB_TYPE)
#include "geant321/gtparm.inc"
#endif
C
#endif
*
* gcparm.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcparm
#include "geant321/incend.inc"
#endif

#endif
