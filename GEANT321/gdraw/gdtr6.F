*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR6(JMAX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Build all optimized tree                                 *
C.    *                                                                *
C.    *       JMAX = total number of nodes (output)                    *
C.    *                                                                *
C.    *    ==>Called by : GDTR0                                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
C.
      DO 10 J=1,NUMNOD
         IQ(JSCA+J)=0
   10 CONTINUE
      JNEW=0
      DO 40 J=1,NUMNOD
         IF (IQ(JDUP+J).EQ.0) THEN
            JX=IQ(JXON+J)
            IF (JX.NE.0) THEN
               IQ(JDUP+JX)=0
               JB=IQ(JBRO+JX)
   20          IF (JB.EQ.0) GO TO 30
               IQ(JDUP+JB)=0
               JB=IQ(JBRO+JB)
               GO TO 20
   30          CONTINUE
            ENDIF
         ELSE
            JNEW=JNEW+1
            IQ(JSCA+J)=JNEW
         ENDIF
   40 CONTINUE
      JMAX=JNEW
      IQ(JNAM1+1)=IQ(JNAM+1)
      IQ(JMOT1+1)=0
      IQ(JBRO1+1)=0
      DO 50 J=1,NUMNOD
         IF (IQ(JSCA+J).NE.0) CALL GDTR5(J)
   50 CONTINUE
      END
