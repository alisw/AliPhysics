*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDMAN(U0,V0)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw a 2D-man at position (U0,V0)                        *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDSPEC                      *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
      DIMENSION U(24),V(24),LEN(4),UN(24,4),VN(24,4)
      DIMENSION UN1(24),VN1(24),UN2(8),VN2(8)
      DIMENSION UN3(4),VN3(4),UN4(12),VN4(12)
C
      EQUIVALENCE (UN1(1),UN(1,1)),(VN1(1),VN(1,1))
      EQUIVALENCE (UN2(1),UN(1,2)),(VN2(1),VN(1,2))
      EQUIVALENCE (UN3(1),UN(1,3)),(VN3(1),VN(1,3))
      EQUIVALENCE (UN4(1),UN(1,4)),(VN4(1),VN(1,4))
C
      SAVE LEN,UN1,VN1,UN2,VN2,UN3,VN3,UN4,VN4
      DATA LEN/24,8,4,12/
C
      DATA UN1/1.3,-0.4,-0.7,0.3,-1.2,-1.6,-1.7,-2.8,-3.6,-1.7,
     +         -1.0,-0.1,0.7,1.9,2.4,2.5,1.4,2.9,3.3,3.1,
     +         1.9,0.9,1.0,1.3/
      DATA VN1/6.2,6.3,4.7,4.6,3.7,2.8,-1.9,-9.5,-9.8,-9.8,
     +         -5.5,-3.6,-5.5,-10.0,-10.4,-9.8,-3.6,-1.2,0.2,1.7,
     +         3.8,4.4,4.7,6.2/
C
      DATA UN2/0.0,0.2,1.3,0.9,0.2,-0.3,-0.8,0.2/
      DATA VN2/0.0,2.6,4.1,4.4,2.6,4.3,3.9,2.6/
C
      DATA UN3/1.8,2.4,1.9,1.8/
      DATA VN3/-1.1,0.4,1.9,-1.1/
C
      DATA UN4/-1.2,-1.7,-2.3,-2.1,-2.2,-2.5,-2.6,-1.7,-1.7,
     +         -2.2,-2.1,-1.6/
      DATA VN4/3.7,2.9,-1.3,-1.4,-1.8,-1.7,-2.8,-3.0,-1.9,
     +         -1.8,-1.4,-1.5/
C.
C.    ------------------------------------------------------------------
C.
      DO 20 I=1,4
         JMAX=LEN(I)
         DO 10 J=1,JMAX
            U(J)=UN(J,I)*10.0*GSCU*GZUA+U0
            V(J)=VN(J,I)*10.0*GSCV*GZVA+V0
C
C             Normalize to PLTRNX,PLTRNY
C
            U(J)=U(J)*PLTRNX/20.
            V(J)=V(J)*PLTRNY/20.
   10    CONTINUE
         CALL GDRAWV(U,V,JMAX)
   20 CONTINUE
C
      END
