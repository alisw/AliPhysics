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
      SUBROUTINE GUSWIM(CHARGE,STEP,VECT,VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to control tracking of one track            *
C.    *       in a magnetic field                                      *
C.    *                                                                *
C.    *    ==>Called by : GTELEC,GTHADR,GTMUON                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gctmed.inc"
      DIMENSION VECT(6),VOUT(6)
C.
C.    ------------------------------------------------------------------
C.
         IF (IFIELD.EQ.3)THEN
            CALL GHELX3(FIELDM*CHARGE,STEP,VECT,VOUT)
         ELSEIF(IFIELD.EQ.2)THEN
            CALL GHELIX(CHARGE,STEP,VECT,VOUT)
         ELSE
            CALL GRKUTA(CHARGE,STEP,VECT,VOUT)
         ENDIF
C
      END
