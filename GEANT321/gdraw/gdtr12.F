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
      SUBROUTINE GDTR12(JJ)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Plot a node of the tree                                  *
C.    *                                                                *
C.    *       JJ = node pointer; negative when tree is being           *
C.    *            scanned from bottom to up; positive when            *
C.    *            it is scanned top to bottom or left to right        *
C.    *            (input)                                             *
C.    *                                                                *
C.    *    ==>Called by : GDTR10, GDTR11                               *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
C.
      J=IABS(JJ)
      IQ(JSCA1+J)=1
      IF (LEVVER.EQ.MAXV.AND.IQ(JXON1+J).NE.0) IQ(JSCA1+J)=-1
      IF (JJ.GT.0) THEN
         LEVHOR=LEVHOR+1
         Q(JULEV+J)=LEVHOR
      ELSE
         JX=IQ(JXON1+J)
         UL=Q(JULEV+JX)
   10    JB=JX
         JX=IQ(JBRO1+JB)
         IF (JX.NE.0) GO TO 10
         UR=Q(JULEV+JB)
         Q(JULEV+J)=(UL+UR)/2.
      ENDIF
      Q(JVLEV+J)=LEVVER
      END
