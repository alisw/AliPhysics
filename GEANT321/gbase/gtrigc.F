*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GTRIGC
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Clear event partition                                    *
C.    *                                                                *
C.    *    ==>Called by : GRUN   ,<GXINT> GINC4                        *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcstak.inc"
#include "geant321/gcflag.inc"
C.
C.    ------------------------------------------------------------------
C.
*
      IF (NJTMAX.LT.0) NJTMAX = -NJTMAX
*
      CALL GRNDMQ(IQ(JRUNG+19),IQ(JRUNG+20),0,'G')
      IQ(JRUNG+29) = IEVENT
*
      IF (JSTAK.NE.0) THEN
         IQ(JSTAK+1) = 0
         IQ(JSTAK+3) = 0
      ENDIF
*
      CALL MZWIPE (IXDIV)
*                                                             END GTRIGC
      END
 
