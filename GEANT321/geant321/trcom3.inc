*
* $Id$
*
* $Log$
* Revision 1.2  1996/07/25 15:01:26  ravndal
* CMZ to CVS conversion artefact removed
*
* Revision 1.1  1996/03/06 15:32:35  mclareni
* Add geane321 includes
*
*
#ifndef CERNLIB_ERPRCDES_TRCOM3_INC
#define CERNLIB_ERPRCDES_TRCOM3_INC
*
*      Sequences of the EMC routines
*
* trcom3.inc
*
      DOUBLE PRECISION A, B, S, TN, T, COSL, SINL, COSP, SINP, COSL1
      INTEGER          NEW
      COMMON /TRCOM3/ A(5,5), B(5,5), S(15),
     +                TN(3), T(5), COSL, SINL, COSP, SINP, COSL1, NEW
*
* trcom3.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
     +    TRCOM3
#include "geant321/incend.inc"
#endif

#endif
