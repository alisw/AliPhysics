*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:44  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CRANDM_INC
#define CERNLIB_GEANT321_CRANDM_INC
*
*
* crandm.inc
*
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION RANDC
C
#endif
#if defined(CERNLIB_SINGLE)
      REAL RANDC
C
#endif

#endif
