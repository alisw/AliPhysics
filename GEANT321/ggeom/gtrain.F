*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GTRAIN(PARU,P)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *              SUBROUTINE GTRAIN(PARU,P)                         *
C.    *    Routine to convert from user supplied parameters PARU to    *
C.    *    internal parameters P for the twisted trapezoidal shape     *
C.    *    GTRA.                                                       *
C.    *          Called by : GSPOSP, GSVOLU                            *
C.    *          A.C.McPherson    22nd April 1985.                     *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION PARU(12),P(30),CORN(3,8)
C.
C.                   -------------------------------------
C.
C                First copy the 12 user parameters into the internal
C                parameter array.
C
      DO 10 I=1,12
      P(I)=PARU(I)
   10 CONTINUE
C
C                Now some trignometrical quantities including the first
C                two extra parameters.
C
      P(13) = TAN(P(2)*DEGRAD)*COS(P(3)*DEGRAD)
      P(14) = TAN(P(2)*DEGRAD)*SIN(P(3)*DEGRAD)
      STWIS = SIN(P(4)*DEGRAD)
      CTWIS = COS(P(4)*DEGRAD)
      DXDY1 = TAN(P(8)*DEGRAD)
      DXDY2 = TAN(P(12)*DEGRAD)
C
C                Next the edges outside the constant z faces (derived
C                from the eight corner points).
C
      P(15) = (CTWIS*(P(5)*DXDY1+P(7)+P(9)*DXDY2+P(11)) +
     +STWIS*(P(5)-P(9)))*0.5
      P(16) = (CTWIS*(P(5)+P(9)) -
     +STWIS*(P(5)*DXDY1+P(7)-P(9)*DXDY2-P(11)))*0.5
      P(17) = P(13) - (CTWIS*(P(5)*DXDY1+P(7)-P(9)*DXDY2-P(11)) +
     +STWIS*(P(5)+P(9)))*0.5/P(1)
      P(18) = P(14) - (CTWIS*(P(5)-P(9)) -
     +STWIS*(P(5)*DXDY1+P(7)+P(9)*DXDY2+P(11)))*0.5/P(1)
C
      P(19) = (CTWIS*(P(5)*DXDY1-P(7)+P(9)*DXDY2-P(11)) +
     +STWIS*(P(5)-P(9)))*0.5
      P(20) = (CTWIS*(P(5)+P(9)) -
     +STWIS*(P(5)*DXDY1-P(7)-P(9)*DXDY2+P(11)))*0.5
      P(21) = P(13) - (CTWIS*(P(5)*DXDY1-P(7)-P(9)*DXDY2+P(11)) +
     +STWIS*(P(5)+P(9)))*0.5/P(1)
      P(22) = P(14) - (CTWIS*(P(5)-P(9)) -
     +STWIS*(P(5)*DXDY1-P(7)+P(9)*DXDY2-P(11)))*0.5/P(1)
C
      P(23) = -(CTWIS*(P(5)*DXDY1+P(6)+P(9)*DXDY2+P(10)) +
     +STWIS*(P(5)-P(9)))*0.5
      P(24) = -(CTWIS*(P(5)+P(9)) -
     +STWIS*(P(5)*DXDY1+P(6)-P(9)*DXDY2-P(10)))*0.5
      P(25) = P(13) + (CTWIS*(P(5)*DXDY1+P(6)-P(9)*DXDY2-P(10)) +
     +STWIS*(P(5)+P(9)))*0.5/P(1)
      P(26) = P(14) + (CTWIS*(P(5)-P(9)) -
     +STWIS*(P(5)*DXDY1+P(6)+P(9)*DXDY2+P(10)))*0.5/P(1)
C
      P(27) = -(CTWIS*(P(5)*DXDY1-P(6)+P(9)*DXDY2-P(10)) +
     +STWIS*(P(5)-P(9)))*0.5
      P(28) = -(CTWIS*(P(5)+P(9)) -
     +STWIS*(P(5)*DXDY1-P(6)-P(9)*DXDY2+P(10)))*0.5
      P(29) = P(13) + (CTWIS*(P(5)*DXDY1-P(6)-P(9)*DXDY2+P(10)) +
     +STWIS*(P(5)+P(9)))*0.5/P(1)
      P(30) = P(14) + (CTWIS*(P(5)-P(9)) -
     +STWIS*(P(5)*DXDY1-P(6)+P(9)*DXDY2-P(10)))*0.5/P(1)
C
C               All done.
C
      DO 970 IL=1,4
      I0=IL*4+11
      CORN(3,IL)=-P(1)
      CORN(1,IL)=P(I0)+P(I0+2)*CORN(3,IL)
      CORN(2,IL)=P(I0+1)+P(I0+3)*CORN(3,IL)
      CORN(3,IL+4)=P(1)
      CORN(1,IL+4)=P(I0)+P(I0+2)*CORN(3,IL+4)
      CORN(2,IL+4)=P(I0+1)+P(I0+3)*CORN(3,IL+4)
  970 CONTINUE
C
#if defined(CERNLIB_DEBUG)
C      WRITE(CHMAIL,1000) P,CORN
C      CALL GMAIL(0,0)
C 1000 FORMAT(' P =',10F10.4,/,4X,10F10.4,/,4X,10F10.4,/,' CORN =',
C     +3F10.4,7X,3F10.4,3(/,7X,3F10.4,7X,3F10.4))
#endif
C
      RETURN
      END
