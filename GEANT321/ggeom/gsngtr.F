*
* $Id$
*
* $Log$
* Revision 1.2  1996/02/27 09:54:45  ravndal
* Double precision for some locals, thanks to G.Poulard
*
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSNGTR(X,P,IACT,SNEXT,SNXT,SAFE,INSIDE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *    Routine to determine the shortest distance from the point   *
C.    *    X(1-3) to the boundary of the shape of type GTRA defined    *
C.    *    by the parameters P along the vector X(4-6). If INSIDE is   *
C.    *    1 then the point is inside the shape and the distance is    *
C.    *    returned as SNEXT. If INSIDE is 0 then the point is         *
C.    *    outside the shape and if the line hits the shape then       *
C.    *    if the new distance is less than the                        *
C.    *    old value of SNEXT the new distance is returned as SNEXT.   *
C.    *                                                                *
C.    *          Called by : GNEXT, GTNEXT                             *
C.    *          A.C.McPherson   22nd April 1985.                      *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
C
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION X0,Y0,DXDZ,DYDZ,A,B,C,DISC,X1,X2,X3,SN,CP,SMALL
#endif
      PARAMETER (SMALL=1E-10)
C
      DIMENSION X(6),P(30),SN(2,5),IOUT(5),X0(4),Y0(4),DXDZ(4),DYDZ(4)
C.
C.                ---------------------------------------------
C.
C
C               Compute Safety distance
C
      DOUBLE PRECISION SL,SL1,SM,SM1
      IF(IACT.LT.3) CALL GSAGTR(X,P,SAFE,INSIDE)
      SNXT=BIG
      IF (IACT .EQ. 0) GO TO 999
      IF (IACT .EQ. 1) THEN
        IF (SNEXT .LT. SAFE) GO TO 999
      ENDIF
C
C               First compute the distance along the line to the
C               boundaries.
C
C               The distance to the planes defined by z=+/-P(1).
C
      IF(X(6).EQ.0.0) THEN
          SN(1,1)=BIG
          SN(2,1)=BIG
          GOTO 10
      ENDIF
      SN(1,1)=(-P(1)-X(3))/X(6)
      SN(2,1)=(P(1)-X(3))/X(6)
      IF(X(6).GT.0.0) GO TO 10
      ST=SN(2,1)
      SN(2,1)=SN(1,1)
      SN(1,1)=ST
   10 CONTINUE
C
C               The distance to the remaining four surfaces.
C
      DO 20 I=1,4
      X0(I)=P(I*4+11)
      Y0(I)=P(I*4+12)
      DXDZ(I)=P(I*4+13)
      DYDZ(I)=P(I*4+14)
   20 CONTINUE
C
      DO 65 I=1,4
      J=I+1
      IF(J.EQ.5) J=1
C
      A=(X(4)-DXDZ(I)*X(6))*(DYDZ(J)-DYDZ(I))*X(6) -
     +(X(5)-DYDZ(I)*X(6))*(DXDZ(J)-DXDZ(I))*X(6)
C
      B=(X(4)-DXDZ(I)*X(6))*(Y0(J)-Y0(I)+(DYDZ(J)-DYDZ(I))*X(3)) +
     +(X(1)-X0(I)-DXDZ(I)*X(3))*(DYDZ(J)-DYDZ(I))*X(6) -
     +(X(5)-DYDZ(I)*X(6))*(X0(J)-X0(I)+(DXDZ(J)-DXDZ(I))*X(3)) -
     +(X(2)-Y0(I)-DYDZ(I)*X(3))*(DXDZ(J)-DXDZ(I))*X(6)
C
      C=(X(1)-X0(I)-DXDZ(I)*X(3))*(Y0(J)-Y0(I)+(DYDZ(J)-DYDZ(I))*X(3))
     + - (X(2)-Y0(I)-DYDZ(I)*X(3))*(X0(J)-X0(I)+(DXDZ(J)-DXDZ(I))*X(3))
C
      IOUT(I+1)=0
      IF(C.GT.0.0) IOUT(I+1)=1
C
C             The solutions are in the normal form:
C             s = (-B+/-SQRT(B*B-4.0*A*C))*0.5/A
C
      SN(1,I+1)=BIG
      SN(2,I+1)=BIG
      IF(ABS(A).GT.1.0E-10) GO TO 30
C
C             A = 0 only one solution.
C
      IF(ABS(B).LT.1.0E-10) GO TO 60
C
      SN(1,I+1)=-C/B
      GO TO 60
C
   30 CONTINUE
      IF(ABS(C).GT.1.0E-10) GO TO 40
      SN(1,I+1)=0.0
      SN(2,I+1)=0.0
      IF(ABS(B).LT.1.0E-10) GO TO 60
      SN(1,I+1)=-C/B
      IF(C.EQ.0.0) SN(1,I+1)=SIGN(SMALL,B)
      SN(2,I+1)=-B/A
      GO TO 50
C
   40 CONTINUE
      DISC=B*B-A*C*4.0
      IF(DISC.LT.0.0) GO TO 60
      IF(DISC.GT.0.0) DISC=SQRT(DISC)
      SN(1,I+1)=(-B-DISC)*0.5/A
      SN(2,I+1)=(-B+DISC)*0.5/A
C
   50 CONTINUE
      IF(SN(2,I+1).GT.SN(1,I+1)) GO TO 60
      ST=SN(2,I+1)
      SN(2,I+1)=SN(1,I+1)
      SN(1,I+1)=ST
C
   60 CONTINUE
C
      DO 65 K=1,2
      IF(ABS(SN(K,I+1)).GT.1.0E+05.OR.ABS(SN(K,I+1)).LT.1.0E-05)
     +GO TO 65
C
      X1=X(1)+SN(K,I+1)*X(4)
      X2=X(2)+SN(K,I+1)*X(5)
      X3=X(3)+SN(K,I+1)*X(6)
      CP=(X1-X0(I)-DXDZ(I)*X3)*(Y0(J)-Y0(I)+(DYDZ(J)-DYDZ(I))*X3)
     + - (X2-Y0(I)-DYDZ(I)*X3)*(X0(J)-X0(I)+(DXDZ(J)-DXDZ(I))*X3)
      CP=CP/SQRT((X0(J)-X0(I)+(DXDZ(J)-DXDZ(I))*X3)**2+
     +   (Y0(J)-Y0(I)+(DYDZ(J)-DYDZ(I))*X3)**2)
C
      IF(ABS(CP).LT.0.0001) GO TO 65
      IF(ABS(CP/SN(K,I+1)).LT.1.0E-06) GO TO 65
      WRITE(CHMAIL,1020) I,K,SN(K,I+1)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1021) X
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1022) X1,X2,X3,CP
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1023) A,B,C,DISC
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1024) INSIDE
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1025) X0
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1026) Y0
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1027) DXDZ
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1028) DYDZ
      CALL GMAIL(0,0)
 1020 FORMAT('0 GSNGTR ERROR - I,K =',2I2,' SN =',E13.5)
 1021 FORMAT(' X =',6F15.6)
 1022 FORMAT(' X1,X2,X3 =',3F15.6,' CP =',E15.6)
 1023 FORMAT(' A =',E15.6,' B =',E15.6,' C =',E15.6,' DISC =',E15.6)
 1024 FORMAT(' INSIDE =',I3)
 1025 FORMAT('   X0 =',4E15.6)
 1026 FORMAT('   Y0 =',4E15.6)
 1027 FORMAT(' DXDZ =',4E15.6)
 1028 FORMAT(' DYDZ =',4E15.6)
C
   65 CONTINUE
C
C
C             Have computed the two distances for the z planes and
C             the four surfaces. Combine them accordingly as to
C             whether the point is inside or outside the shape.
C
      IF(INSIDE.EQ.0) GO TO 80
C
C             Point is inside shape.
C
      DO 70 I=1,5
      DO 70 J=1,2
      IF(SN(J,I).GT.0.0.AND.SN(J,I).LT.SNXT) SNXT=SN(J,I)
   70 CONTINUE
      GO TO 999
C
   80 CONTINUE
C
C             Point is outside shape.
C
      IOUT(1)=0
      IF(ABS(X(3)).GT.P(1)) IOUT(1)=1
C
C             For each of five sets of SN and IOUT, IOUT(I) equal to 1
C             indicates that the point is outside the shape by the Ith
C             test, SN(1,I) is the distance to the first change in the
C             test and SN(2,I) is the distance to the second change.
C             The remaining logic just attempts to find a distance when
C             the line is inside by all five tests, bearing in mind that
C             for some tests the line can start inside, leave and return
C             inside.
C
      SL=-1.0
      SM=BIG
      SM1=BIG
      DO 100 I=1,5
      IF(IOUT(I).EQ.0) GO TO 90
      IF(SN(2,I).LT.0.0) GO TO 999
      IF(SN(1,I).LT.0.0.AND.SN(2,I).GT.SL) SL=SN(2,I)
      IF(SN(1,I).GT.SL) SL=SN(1,I)
      IF(SN(1,I).GE.0.0.AND.SN(2,I).LT.SM) SM=SN(2,I)
      GO TO 100
   90 CONTINUE
      IF(SN(1,I).LT.0.0.AND.SN(2,I).GE.0.0.AND.SN(2,I).LT.SM) SM=SN(2,I)
      IF(SN(1,I).LT.0.0.OR.SN(1,I).GT.SM1) GO TO 100
      IF(SN(1,I).GE.SN(2,I)) GO TO 100
      SM1=SN(1,I)
      SL1=SN(2,I)
  100 CONTINUE
C
C             SL is the largest of the five distances to the first
C             time the line is inside. SM is the smallest to the
C             last time the point is inside. SM1 is the smallest
C             distance to when the line is temporarily outside
C             one of the tests.
C
      IF(SM.LE.SL) GO TO 999
      IF(SM1.GT.SL) GO TO 130
C
  110 CONTINUE
C
C             In this loop SL is updated by the return after SM1
C             if SM1 is less than SL.
C
      SL=SL1
      IF(SM.LE.SL) GO TO 999
      SM1=SM
C
      DO 120 I=1,5
      IF(IOUT(I).EQ.1) GO TO 120
      IF(SN(2,I).LE.SL.OR.SN(1,I).GT.SM1) GO TO 120
      IF(SN(1,I).GE.SN(2,I)) GO TO 120
      SM1=SN(1,I)
      SL1=SN(2,I)
  120 CONTINUE
C
      IF(SM1.GT.SL) GO TO 130
C
      GO TO 110
  130 CONTINUE
C
      IF(SL.LT.SNXT) SNXT=SL
C
  999 CONTINUE
      END
