*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GRLEAS(JBANK)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Release unused space in data structure pointed by JBANK. *
C.    *       JBANK can be either JHITS or JDIGI                       *
C.    *                                                                *
C.    *    ==>Called by : GSAVE                                        *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF(JBANK.LE.0)GO TO 99
      NS=IQ(JBANK-2)
      IF(NS.LE.0)GO TO 99
      DO 20 IS=1,NS
         JB=LQ(JBANK-IS)
         IF(JB.GT.0)THEN
            ND=IQ(JB-1)
            IF(ND.GT.0)THEN
               DO 10 ID=1,ND
                  JBD=LQ(JB-ID)
                  IF(JBD.GT.0)THEN
                     NPUSH=IQ(JBD-1)-IQ(JB+ID)
                     IF(NPUSH.GT.0)CALL MZPUSH(IXDIV,JBD,0,-NPUSH,'I')
                  ENDIF
  10           CONTINUE
            ENDIF
         ENDIF
  20  CONTINUE
C
  99  RETURN
      END
