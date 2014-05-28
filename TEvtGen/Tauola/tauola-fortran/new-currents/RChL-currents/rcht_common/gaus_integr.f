* THIS IS ITERATIVE INTEGRATION PROCEDURE
* ORIGINATES  PROBABLY FROM CERN LIBRARY
* IT SUBDIVIDES INEGRATION RANGE UNTIL REQUIRED PRECISION IS REACHED
* PRECISION IS A DIFFERENCE FROM 8 AND 16 POINT GAUSS ITEGR. RESULT
* EEPS POSITIVE TREATED AS ABSOLUTE PRECISION
* EEPS NEGATIVE TREATED AS RELATIVE PRECISION
      DOUBLE PRECISION FUNCTION GAUS(F,A,B,EEPS)
*     *************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      DIMENSION W(12),X(12)
      DATA CONST /1.0D-19/
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/

C-----------------------------------------------------------------------------
C NEW INTEGRATION METHOD - constant step size (provided by fourth parameter)
C comment out following 5 lines to return to old solutions of advantages too
C-----------------------------------------------------------------------------
      RESULT=0
      EPS=0
c      CALL STEPGAUSS(F,A,B,6,RESULT,EPS)
      CALL CHANGEGAUSS(F,A,B,2,RESULT,EPS)
      GAUS=RESULT
      RETURN
C-----------------------------------------------------------------------------
C PREVIOUS INTEGRATION METHOD - step size depends on EPSSQ
C-----------------------------------------------------------------------------
      EPS=dABS(EEPS)
      DELTA=CONST*dABS(A-B)
      GAUS=0.d0
      AA=A
    5 Y=B-AA
      IF(dABS(Y) .LE. DELTA) RETURN
    2 BB=AA+Y
      C1=0.5d0*(AA+BB)
      C2=C1-AA
      S8=0.d0
      S16=0.d0
      DO 1 I=1,4
      U=X(I)*C2
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 3 I=5,12
      U=X(I)*C2
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(EEPS.LT.0D0) THEN
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4
      ELSE
        IF(dABS(S16-S8) .GT. EPS) GO TO 4
      ENDIF
      GAUS=GAUS+S16
      AA=BB
      GO TO 5
    4 Y=0.5d0*Y
      IF(dABS(Y) .GT. DELTA) GOTO 2
      PRINT 7
      GAUS=0.d0
      RETURN
    7 FORMAT(1X,36HGAUS  ... TOO HIGH ACCURACY REQUIRED)
      END


      DOUBLE PRECISION FUNCTION GAUS2(F,A,B,EEPS)
*     *************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      DIMENSION W(12),X(12)
      DATA CONST /1.0D-19/
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/

C-----------------------------------------------------------------------------
C NEW INTEGRATION METHOD - constant step size (provided by fourth parameter)
C comment out following 5 lines to return to old solutions of advantages too
C-----------------------------------------------------------------------------
      RESULT=0
      EPS=0
c      CALL STEPGAUSS2(F,A,B,6,RESULT,EPS)
      CALL CHANGEGAUSS2(F,A,B,2,RESULT,EPS)
      GAUS2=RESULT
      RETURN
C-----------------------------------------------------------------------------
C PREVIOUS INTEGRATION METHOD - step size depends on EPSSQ
C----------------------------------------------------------------------------- 
      EPS=dABS(EEPS)
      DELTA=CONST*dABS(A-B)
      GAUS2=0.D0
      AA=A
    5 Y=B-AA
      IF(dABS(Y) .LE. DELTA) RETURN
    2 BB=AA+Y
      C1=0.5d0*(AA+BB)
      C2=C1-AA
      S8=0.d0
      S16=0.d0
      DO 1 I=1,4
      U=X(I)*C2
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 3 I=5,12
      U=X(I)*C2
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(EEPS.LT.0D0) THEN
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4
      ELSE
        IF(dABS(S16-S8) .GT. EPS) GO TO 4
      ENDIF
      GAUS2=GAUS2+S16
      AA=BB
      GO TO 5
    4 Y=0.5d0*Y
      IF(dABS(Y) .GT. DELTA) GOTO 2
      PRINT 7
      GAUS2=0.D0
      RETURN
    7 FORMAT(1X,36HGAUS2 ... TOO HIGH ACCURACY REQUIRED)
      END

      DOUBLE PRECISION FUNCTION GAUS3(F,A,B,EEPS)
*     *************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      DIMENSION W(12),X(12)
      DATA CONST /1.0D-19/
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/

C-----------------------------------------------------------------------------
C NEW INTEGRATION METHOD - constant step size (provided by fourth parameter)
C comment out following 5 lines to return to old solutions of advantages too
C-----------------------------------------------------------------------------
      RESULT=0
      EPS=0
c      CALL STEPGAUSS3(F,A,B,6,RESULT,EPS)
      CALL CHANGEGAUSS3(F,A,B,2,RESULT,EPS)
      GAUS3=RESULT
      RETURN
C-----------------------------------------------------------------------------
C PREVIOUS INTEGRATION METHOD - step size depends on EPSSQ
C-----------------------------------------------------------------------------
      EPS=dABS(EEPS)
      DELTA=CONST*dABS(A-B)
      GAUS3=0.D0
      AA=A
    5 Y=B-AA
      IF(dABS(Y) .LE. DELTA) RETURN
    2 BB=AA+Y
      C1=0.5d0*(AA+BB)
      C2=C1-AA
      S8=0.d0
      S16=0.d0
      DO 1 I=1,4
      U=X(I)*C2
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 3 I=5,12
      U=X(I)*C2
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(EEPS.LT.0D0) THEN
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4
      ELSE
        IF(dABS(S16-S8) .GT. EPS) GO TO 4
      ENDIF
      GAUS3=GAUS3+S16
      AA=BB
      GO TO 5
    4 Y=0.5d0*Y
      IF(dABS(Y) .GT. DELTA) GOTO 2
      PRINT 7
      GAUS3=0.D0
      RETURN
    7 FORMAT(1X,36HGAUS3 ... TOO HIGH ACCURACY REQUIRED)
      END

C--------------------------------------------------------------------------------------------------
C--------------------------------------------------------------------------------------------------
C--------------------------------------------------------------------------------------------------
C OBSOLETE 1.13      901108  3.40 CERN PROGRAM LIBRARY OBSOLETE PAM
C          CERN OBSOLETE PROGRAMS AND SUBROUTINES PAM-FILE.
C          MANY OF THE ROUTINES ARE FOR ONE FORTRAN DIALECT ONLY.
C          THE FLAG PATCH F77 SHOULD BE USED FOR THE FORTRAN 77 VERSIONS
C          SHOULD THEY EXIST - IF NOT FORTRAN 4 WILL BE OBTAINED.
C          THE CERN LIBRARIES GENLIB AND KERNLIB ARE USUALLY REQUIRED.
C          OBSOLETE ROUTINES OF THE PROGRAM LIBRARY ARE OFTEN REQUIRED.
C
      SUBROUTINE STEPGAUSS(F,A,B,LAMBDA,RESULT,EPS)
      implicit real*8 (a-h,o-z)
      DIMENSION W(12),X(12)
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      RESULT=0.
      EPS=0.
      AU=A
      C=(B-A)/FLOAT(LAMBDA)
      DO 1 I = 1,LAMBDA
      FI=I
      AO=A+FI*C
      C1=0.5*(AU+AO)
      C2=C1-AU
      S8=0.
      S16=0.
      DO 2 J = 1,4
      U=X(J)*C2
    2 S8=S8+W(J)*(F(C1+U)+F(C1-U))
      DO 3 J = 5,12
      U=X(J)*C2
    3 S16=S16+W(J)*(F(C1+U)+F(C1-U))
      RESULT=RESULT+C2*S16
      EPS=EPS+ABS(C2*(S16-S8))
    1 AU=AO
      RETURN
      END

      SUBROUTINE STEPGAUSS2(F,A,B,LAMBDA,RESULT,EPS)
      implicit real*8 (a-h,o-z)
      DIMENSION W(12),X(12)
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      RESULT=0.
      EPS=0.
      AU=A
      C=(B-A)/FLOAT(LAMBDA)
      DO 1 I = 1,LAMBDA
      FI=I
      AO=A+FI*C
      C1=0.5*(AU+AO)
      C2=C1-AU
      S8=0.
      S16=0.
      DO 2 J = 1,4
      U=X(J)*C2
    2 S8=S8+W(J)*(F(C1+U)+F(C1-U))
      DO 3 J = 5,12
      U=X(J)*C2
    3 S16=S16+W(J)*(F(C1+U)+F(C1-U))
      RESULT=RESULT+C2*S16
      EPS=EPS+ABS(C2*(S16-S8))
    1 AU=AO
      RETURN
      END

      SUBROUTINE STEPGAUSS3(F,A,B,LAMBDA,RESULT,EPS)
      implicit real*8 (a-h,o-z)
      DIMENSION W(12),X(12)
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      RESULT=0.
      EPS=0.
      AU=A
      C=(B-A)/FLOAT(LAMBDA)
      DO 1 I = 1,LAMBDA
      FI=I
      AO=A+FI*C
      C1=0.5*(AU+AO)
      C2=C1-AU
      S8=0.
      S16=0.
      DO 2 J = 1,4
      U=X(J)*C2
    2 S8=S8+W(J)*(F(C1+U)+F(C1-U))
      DO 3 J = 5,12
      U=X(J)*C2
    3 S16=S16+W(J)*(F(C1+U)+F(C1-U))
      RESULT=RESULT+C2*S16
      EPS=EPS+ABS(C2*(S16-S8))
    1 AU=AO
      RETURN
      END

*******************************************************************************
* FUNCTIONS FOR CHANGE OF VARIABLES *******************************************
*******************************************************************************
      DOUBLE PRECISION FUNCTION F_CHANGE(X,F,AMS1,AMS2)
      IMPLICIT NONE
      EXTERNAL F
      DOUBLE PRECISION F
      DOUBLE PRECISION X, NEW_X
      DOUBLE PRECISION AMS1,AMS2
      DOUBLE PRECISION ALP1,ALP2,ALP
      DOUBLE PRECISION AMRX,GAMRX
      DATA AMRX/0.77/,GAMRX/1.8/

      ALP1  = ATAN((AMS1-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP2  = ATAN((AMS2-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP   = ALP1 + X*(ALP2-ALP1)            ! change of variables

      NEW_X = AMRX**2+AMRX*GAMRX*TAN(ALP)     ! second change of variables

      F_CHANGE = F(NEW_X)

! Jacobian for change of variables
      F_CHANGE = F_CHANGE * (ALP2-ALP1) * ((NEW_X-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)

      RETURN
      END

      DOUBLE PRECISION FUNCTION F_CHANGE2(X,F,AMS1,AMS2)
      IMPLICIT NONE
      EXTERNAL F
      DOUBLE PRECISION F
      DOUBLE PRECISION X, NEW_X
      DOUBLE PRECISION AMS1,AMS2
      DOUBLE PRECISION ALP1,ALP2,ALP
      DOUBLE PRECISION AMRX,GAMRX
      DATA AMRX/0.77/,GAMRX/1.8/

      ALP1  = ATAN((AMS1-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP2  = ATAN((AMS2-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP   = ALP1 + X*(ALP2-ALP1)            ! change of variables

      NEW_X = AMRX**2+AMRX*GAMRX*TAN(ALP)     ! second change of variables

      F_CHANGE2 = F(NEW_X)

! Jacobian for change of variables
      F_CHANGE2 = F_CHANGE2 * (ALP2-ALP1) * ((NEW_X-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)

      RETURN
      END

      DOUBLE PRECISION FUNCTION F_CHANGE3(X,F,AMS1,AMS2)
      IMPLICIT NONE
      EXTERNAL F
      DOUBLE PRECISION F
      DOUBLE PRECISION X, NEW_X
      DOUBLE PRECISION AMS1,AMS2
      DOUBLE PRECISION ALP1,ALP2,ALP
      DOUBLE PRECISION AMRX,GAMRX
      DATA AMRX/0.77/,GAMRX/1.8/

      ALP1  = ATAN((AMS1-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP2  = ATAN((AMS2-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP   = ALP1 + X*(ALP2-ALP1)            ! change of variables

      NEW_X = AMRX**2+AMRX*GAMRX*TAN(ALP)     ! second change of variables

      F_CHANGE3 = F(NEW_X)

! Jacobian for change of variables
      F_CHANGE3 = F_CHANGE3 * (ALP2-ALP1) * ((NEW_X-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)

      RETURN
      END

*******************************************************************************
* INTEGRATION WITH CHANGED VARIABLES ******************************************
*******************************************************************************
      SUBROUTINE CHANGEGAUSS(F,AMS1,AMS2,LAMBDA,RESULT,EPS)
      implicit real*8 (a-h,o-z)
      EXTERNAL F,F2
      DIMENSION W(12),X(12)
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      RESULT=0.
      EPS=0.

! Range is changed from [AMS1,AMS2] to [0,1]
! F, AMS1 and AMS2 are parameters of F_CHANGE
      A=0
      B=1

      AU=A
      C=(B-A)/FLOAT(LAMBDA)
      DO 1 I = 1,LAMBDA
      FI=I
      AO=A+FI*C
      C1=0.5*(AU+AO)
      C2=C1-AU
      S8=0.
      S16=0.
      DO 2 J = 1,4
      U=X(J)*C2
    2 S8=S8+W(J)*(F_CHANGE(C1+U,F,AMS1,AMS2)+F_CHANGE(C1-U,F,AMS1,AMS2))
      DO 3 J = 5,12
      U=X(J)*C2
    3 S16=S16+W(J)*(F_CHANGE(C1+U,F,AMS1,AMS2)+F_CHANGE(C1-U,F,AMS1,AMS2))
      RESULT=RESULT+C2*S16
      EPS=EPS+ABS(C2*(S16-S8))
    1 AU=AO
      RETURN
      END

      SUBROUTINE CHANGEGAUSS2(F,AMS1,AMS2,LAMBDA,RESULT,EPS)
      implicit real*8 (a-h,o-z)
      EXTERNAL F,F2
      DIMENSION W(12),X(12)
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      RESULT=0.
      EPS=0.

! Range is changed from [AMS1,AMS2] to [0,1]
! F, AMS1 and AMS2 are parameters of F_CHANGE
      A=0
      B=1

      AU=A
      C=(B-A)/FLOAT(LAMBDA)
      DO 1 I = 1,LAMBDA
      FI=I
      AO=A+FI*C
      C1=0.5*(AU+AO)
      C2=C1-AU
      S8=0.
      S16=0.
      DO 2 J = 1,4
      U=X(J)*C2
    2 S8=S8+W(J)*(F_CHANGE2(C1+U,F,AMS1,AMS2)+F_CHANGE2(C1-U,F,AMS1,AMS2))
      DO 3 J = 5,12
      U=X(J)*C2
    3 S16=S16+W(J)*(F_CHANGE2(C1+U,F,AMS1,AMS2)+F_CHANGE2(C1-U,F,AMS1,AMS2))
      RESULT=RESULT+C2*S16
      EPS=EPS+ABS(C2*(S16-S8))
    1 AU=AO
      RETURN
      END

      SUBROUTINE CHANGEGAUSS3(F,AMS1,AMS2,LAMBDA,RESULT,EPS)
      implicit real*8 (a-h,o-z)
      EXTERNAL F,F2
      DIMENSION W(12),X(12)
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      RESULT=0.
      EPS=0.

! Range is changed from [AMS1,AMS2] to [0,1]
! F, AMS1 and AMS2 are parameters of F_CHANGE
      A=0
      B=1

      AU=A
      C=(B-A)/FLOAT(LAMBDA)
      DO 1 I = 1,LAMBDA
      FI=I
      AO=A+FI*C
      C1=0.5*(AU+AO)
      C2=C1-AU
      S8=0.
      S16=0.
      DO 2 J = 1,4
      U=X(J)*C2
    2 S8=S8+W(J)*(F_CHANGE3(C1+U,F,AMS1,AMS2)+F_CHANGE3(C1-U,F,AMS1,AMS2))
      DO 3 J = 5,12
      U=X(J)*C2
    3 S16=S16+W(J)*(F_CHANGE3(C1+U,F,AMS1,AMS2)+F_CHANGE3(C1-U,F,AMS1,AMS2))
      RESULT=RESULT+C2*S16
      EPS=EPS+ABS(C2*(S16-S8))
    1 AU=AO
      RETURN
      END
