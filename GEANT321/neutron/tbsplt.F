*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE TBSPLT(A,E,NP,Y)
C       THIS ROUTINE DETERMINES A CROSS SECTION AT A GIVEN
C       ENERGY FROM A CROSS SECTION VERSUS ENERGY TABLE USING
C       A TABLE SPLITTING METHOD
      DIMENSION A(1)
      SAVE
      IPP=1
      IF(E.LE.A(1))GO TO 40
      IF(E.GE.A(2*NP-1))GO TO 50
      INDXH=NP
      INDXL=0
   10 IF(INDXL+1.EQ.INDXH)GO TO 30
      J=(INDXH+INDXL)/2
      N=2*J-1
      IF(E.LE.A(N))GO TO 20
      INDXL=J
      GO TO 10
   20 INDXH=J
      GO TO 10
   30 N=2*INDXH-1
      Y=A(N-1)+(E-A(N-2))*(A(N+1)-A(N-1))/(A(N)-A(N-2))
      RETURN
   40 Y=A(IPP+1)
      RETURN
   50 Y=A(2*NP)
      RETURN
      END
