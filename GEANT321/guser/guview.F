*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUVIEW(IDS,IVS,ICS,IVIEW)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine for interactive version                     *
C.    *                                                                *
C.    *    ==>Called by : <GXINT>, GINC1                               *
C.    *                                                                *
C.    ******************************************************************
C.
C.
C.    ------------------------------------------------------------------
C.
      CHARACTER*(*) ICS
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
#if !defined(CERNLIB_USRJMP)
      CALL GUSTEP
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUSTEP)
#endif
      END
