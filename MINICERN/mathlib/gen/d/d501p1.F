*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:19  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE D501P1(K,M,NC,X,NX,Y,SY,MODE,EPS0,EPS,MAXIT,IPRT,
     +                  N,A,AL,AU,NERROR,VERS)

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE CHECKS THE VALUES OF INPUT PARAMETERS OF THE
*   SUBROUTINES  DSUMSQ, DFUNFT, DMAXLK  DEPENDING ON THE VALUE OF
*   THE PARAMETER  VERS.
*   IF  IPRT < 0  ALL VALUES OF THE INPUT PARAMETERS ARE PRINTED.
*
*************************************************************************

#include "gen/imp64.inc"
      CHARACTER VERS*6,TIT1(3)*35,TIT2(0:1)*35
      DIMENSION X(*),Y(*),SY(*),A(*),AL(*),AU(*)
      PARAMETER (Z1 = 1, R10 = Z1/10)

      DATA TIT1(1) /'MINIMIZATION OF A SUM OF SQUARES'/
      DATA TIT1(2) /'LEAST-SQUARES DATA FITTING'/
      DATA TIT1(3) /'MAXIMUM LIKELIHOOD ESTIMATION'/

      DATA TIT2(0) /'APPROXIMATE DERIVATIVES (MODE = 0)'/
      DATA TIT2(1) /'ANALYTICAL DERIVATIVES (MODE = 1)'/

      IF(IPRT .NE. 0) THEN
       IF(VERS .EQ. 'DSUMSQ') IV=1
       IF(VERS .EQ. 'DFUNFT') IV=2
       IF(VERS .EQ. 'DMAXLK') IV=3

       WRITE(6,1000) VERS,TIT1(IV),TIT2(MODE)
      ENDIF

************************************************************************
*   PRINT INPUT PARAMETERS (IF IPRT .LT. 0)
************************************************************************

      IF(IPRT .LT. 0) THEN
       WRITE(6,1010) VERS,M,N
       IF(VERS .NE. 'DMAXLK') WRITE(6,1020) NC
       IF(VERS .NE. 'DSUMSQ') WRITE(6,1030) K,NX
       WRITE(6,1040) MAXIT,MODE,IPRT,EPS
      ENDIF

************************************************************************
*   CHECK VALUES OF INPUT PARAMETERS, AND PRINT THEM (IF IPRT .LT. 0)
************************************************************************

      NERROR=0

      IF(     MAXIT .LT. 1
     1   .OR. K     .LT. 1
     2   .OR. N     .LT. 1
     3   .OR. M     .LT. N
     4   .OR. NC    .LT. N
     5   .OR. NX    .LT. K ) THEN
       NERROR=1
       RETURN
      ENDIF

      IF(IPRT .LT. 0) THEN
       WRITE(6,1050) (AL(I), I=1,N)
       WRITE(6,1060) (AU(I), I=1,N)
       WRITE(6,1070) (A(I),  I=1,N)
      IF(VERS .NE. 'DSUMSQ')WRITE(6,1080)((X(I),I=L,M*NX,NX),L=1,K)
       IF(VERS .EQ. 'DFUNFT') THEN
        WRITE(6,1090) (Y(I), I=1,M)
        WRITE(6,1100) (SY(I),I=1,M)
       ENDIF
      ENDIF

      DO 10 I=1,N
      IF(AL(I) .GT. AU(I)) THEN
       NERROR=1
       RETURN
      ENDIF
   10 CONTINUE

************************************************************************
*   IF VALUES OF THE PARAMETERS A, SY, MODE OR EPS ARE NOT PRACTICABLE
*   SET RECOMMENDED VALUES FOR THIS PARAMETERS
************************************************************************

      DO 20 I=1,N
      IF(A(I) .GT. AU(I)) A(I)=AU(I)
      IF(A(I) .LT. AL(I)) A(I)=AL(I)
   20 CONTINUE

      IF (VERS .EQ. 'DFUNFT') THEN
       DO 30 I = 1,M
       IF(SY(I) .LE. 0) THEN
        NERROR=1
        RETURN
       ENDIF
   30  CONTINUE
      ENDIF

CC    IF(STEP .LE. 0) STEP=1
      IF(MODE .NE. 1) MODE=0
      IF (EPS .LT. EPS0  .OR.  EPS .GT. R10) EPS=10*EPS0

      RETURN

 1000 FORMAT(7(/),30X,'MATHLIB PACKAGE   D501   VERSION 15.03.93'//
     1       30X,'PACKAGE LEAMAX  ****  ROUTINE ',A6,' ****'///
     2       15X,A35,A35//)
 1010 FORMAT(' INPUT  OF  ',A6,' :'//'  M :',I5,6X,'N :',I5)
 1020 FORMAT('  NC:',I5)
 1030 FORMAT('  K :',I5,6X,'NX:',I5)
 1040 FORMAT('  MAXIT :',I5,8X,'MODE :',I5,8X,'IPRT :',I5/
     +       '  EPS   :',1PD11.1/)
 1050 FORMAT(/'  AL :',/(5(1PD15.5)))
 1060 FORMAT( '  AU :',/(5(1PD15.5)))
 1070 FORMAT( '  A :', /(5(1PD15.5)))
 1080 FORMAT( '  X :', /(5(1PD15.5)))
 1090 FORMAT( '  Y :', /(5(1PD15.5)))
 1100 FORMAT( '  SY :',/(5(1PD15.5)))

      END



