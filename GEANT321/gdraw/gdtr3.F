*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR3(JM)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Collect multiplicity for one level of JVOLUM structure   *
C.    *                                                                *
C.    *       JM = mother node (input)                                 *
C.    *                                                                *
C.    *    ==>Called by : GDTR4                                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
      JX=IQ(JXON+JM)
      IF (JX.NE.0) THEN
   10    CONTINUE
         JXX=IQ(JXON+JM)
         IF (IQ(JSCA+JX).EQ.0) THEN
   20       CONTINUE
            IF (JXX.EQ.0) GO TO 30
            IF (IQ(JNAM+JXX).EQ.IQ(JNAM+JX)) THEN
***             IQ(JDUP+JX)=IQ(JDUP+JX)+1
               IQ(JDUP+JX)=IQ(JDUP+JX)+IQ(JPSM+JXX)
               IQ(JSCA+JXX)=-1
            ENDIF
            JXX=IQ(JBRO+JXX)
            GO TO 20
   30       CONTINUE
         ENDIF
         JX=IQ(JBRO+JX)
         IF (JX.NE.0) GO TO 10
      ENDIF
      IQ(JSCA+JM)=1
      END
