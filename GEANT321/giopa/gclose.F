*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:16  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GCLOSE (LUN,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to close I/O units                               *
C.    *                                                                *
C.    *       LUN      Logical unit number                             *
C.    *       IER      error flag                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGLAST                               *
C.    *       Authors    R.Brun, F.Carena  *********                   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
C.
C.    ------------------------------------------------------------------
C.
      IER=0
      IF(LUN.EQ.0)THEN
         CALL FZENDO(LUN,'T')
         CALL FZENDI(LUN,'T')
      ELSE
         DO 10 I=1,NUNITS
            IF(LUN.EQ.ABS(LUNITS(I)))THEN
               IF(LUNITS(I).LT.0)THEN
                  CALL FZENDI(LUN,'T')
               ELSE
                  CALL FZENDO(LUN,'T')
               ENDIF
               LUNITS(I)=0
               GOTO 999
            ENDIF
   10    CONTINUE
         IER=1
      ENDIF
  999 END
