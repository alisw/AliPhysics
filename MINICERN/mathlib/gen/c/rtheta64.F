*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:12  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DTHETA( K, X, Q )
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RTHETA/DTHETA')
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RTHETA( K, X, Q )
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RTHETA')
#endif
      PARAMETER ( TWO=2, HALF=1/TWO, FOURTH=HALF**2 )
      PARAMETER ( PI = 3.14159 26535 89793 24 D0, PISQ=PI**2 )
      PARAMETER ( IADIM = 20 )
      DIMENSION A1(0:IADIM), A3(0:IADIM)
C  Set NBITS = number of mantissa bits.
#if !defined(CERNLIB_DOUBLE)
      PARAMETER ( NBITS = 48 )
C  Set QZ to the value specified in the comments.
      PARAMETER ( QZ = 0.7504 )
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER ( NBITS = 56 )
C  Set QZ to the value specified in the comments.
      PARAMETER ( QZ = 0.7809 )
#endif
C
C  The following statements set various constants (see comments).
      PARAMETER ( ELOW = -0.693*NBITS )
      PARAMETER ( QZSQ = QZ**2 )
      PARAMETER ( C10 = QZSQ, C11 = 2*QZ*(1-QZ) )
      PARAMETER ( C30 = QZSQ**2, C31 = 2*QZSQ*(1-QZSQ) )
C
      CHARACTER*130 ERRTXT

      DATA  QPRV1A, QPRV1B, QPRV3A, QPRV3B / 4*-9999. /
      SAVE A1, A3, QPRV1A, QPRV1B, QPRV3A, QPRV3B, B3, R, COEFF
      SAVE N1, N3
C
C  *********************************************************************
C
C  function is set equal to the value of theta(K,pi*X,Q), where
C  theta(k,u,q) is the theta function with subscript k (=1,2,3,4)
C  as defined in "Handbook of Mathematical Functions", M.Abramowitz
C  & I.A. Stegun (eds.), Washington, 1964.  K=0 is accepted as equiv-
C  alent to K=4.
C
C  The constant QZ in the parameter statement should be set as follows:
C           QZ = exp( (pi**2)/log(eps/3) ),
C  where eps=2**(-NBITS).
C
C  The constants Cj0 and Cj1 (j=1,3) then have values such that the two
C  inequalities
C
C           0 .le. t .le. 1/2,
C           0 .le. q .le. cj0 + cj1*t,
C
C  define a region within which the theta function can be approximated
C  with relative error less than eps by using either one (k=1,2) or two
C  (k=3,4) terms of the transformed series.
C
C  The constant ELOW has a value such that 1-EXP(u) is equal to 1 to
C  machine accuracy if u.lt.ELOW.
C
C  NOTE.  The value of NBITS and the corresponding QZ do not have to be
C         set separately for each computer.  All that is required is
C         that the number of mantissa bits shall not be greater than
C         NBITS.  Too large a value of NBITS merely causes the program
C         to run somewhat slower for certain values of X and Q.
C
C  (Version 11.02.1992)
C
C  *********************************************************************
C
C--Start.  Check Q and select basic function theta1 or theta3.
      IF ( Q.LT.0 .OR. Q.GT.1 ) THEN
         WRITE(ERRTXT,101) K,X,Q
         CALL MTLPRT(NAME,'C349.1',ERRTXT)
         RESULT = 0
         GO TO 50
      ENDIF
      IF ( K.EQ.1 ) THEN
         T = X
         GO TO 10
      ELSEIF ( K.EQ.2 ) THEN
         T = HALF - ABS(X)
         GO TO 10
      ELSEIF ( K.EQ.3 ) THEN
         T = X
         GO TO 30
      ELSEIF ( K.EQ.4 .OR. K.EQ.0 ) THEN
         T = HALF - ABS(X)
         GO TO 30
      ELSE
         WRITE(ERRTXT,102) K,X,Q
         CALL MTLPRT(NAME,'C349.2',ERRTXT)
         RESULT = 0
         GO TO 50
      ENDIF
C----------------------------------------------------------------------
C
C--Basic function THETA1.
C
C  Reduce argument T to the interval (0,1/2).
   10 FSIGN = +1
      IF ( T.LT.0 ) FSIGN = -FSIGN
      T = ABS(T)
      TINT = AINT(T)
      IF ( MOD(TINT,TWO) .NE. 0 ) FSIGN = -FSIGN
      T = T - TINT
      IF ( T.GE.HALF ) THEN
         T = 1 - T
      ELSEIF( T.EQ.0 ) THEN
         RESULT = 0
         GO TO 50
      ENDIF
      IF ( Q .GT. C10+C11*T ) GO TO 14
C
C  Classical series.
C  (If Q has changed value, compute N1 and A1(j), j=0,1,...,N1)
      IF ( Q.NE.QPRV1A ) THEN
         AZERO = 2*SQRT(SQRT(Q))
         QSQ = Q**2
         P = -1
         ATEMP = 1
         DO 11 J = 0,IADIM
            A1(J) = AZERO*ATEMP
            N1 = J
            P = QSQ*P
            ATEMP = P*ATEMP
            IF ( 1-ABS(ATEMP) .EQ. 1 ) GO TO 12
   11    CONTINUE
   12    QPRV1A = Q
      ENDIF
C  (Compute the series using forwards recurrence for sines)
      S = SIN(PI*T)
      ALPHA = 2 - 4*(S**2)
      U = -S
      SUM = A1(0)*S
      DO 13 J = 1,N1
         V = U
         U = S
         S = ALPHA*U - V
         SUM = SUM + A1(J)*S
   13 CONTINUE
      RESULT = FSIGN*SUM
      GO TO 50
C
C  One-term sinh approximation.
   14 IF ( Q.EQ.1 ) THEN
         RESULT = 0
      IF ( T.EQ.HALF ) THEN
         WRITE(ERRTXT,103) K,X,Q
         CALL MTLPRT(NAME,'C349.3',ERRTXT)
      ENDIF
         GO TO 50
      ENDIF
C  (If Q has changed value, compute R and COEFF)
      IF ( Q.NE.QPRV1B ) THEN
         ABSLOG = ABS(LOG(Q))
         R = PISQ/ABSLOG
         COEFF = SQRT(PI/ABSLOG)
         QPRV1B = Q
      ENDIF
C  (Compute function)
      U = R*T
      IF ( U.LT.HALF ) THEN
         E = -R*( T**2 + FOURTH )
         ETERM = EXP(E)
         IF ( ETERM.NE.0 ) THEN
            RESULT = 2*COEFF*ETERM*SINH(U)
         ELSE
            RESULT = EXP( E + LOG(2*COEFF) )*SINH(U)
         ENDIF
      ELSE
         E = -R*(T-HALF)**2
         ETERM = EXP(E)
         IF ( ETERM.NE.0 ) THEN
            RESULT = COEFF*ETERM
         ELSE
            RESULT = EXP( E + LOG(COEFF) )
         ENDIF
         E = -2*U
         IF ( E.GT.ELOW ) RESULT = RESULT*( 1 - EXP(E) )
      ENDIF
      RESULT = FSIGN*RESULT
      GO TO 50
C----------------------------------------------------------------------
C
C--Basic function THETA3.
C
C  Reduce argument T to the interval (0,1/2).
   30 T = ABS(T) - AINT(ABS(T))
      IF ( T.GE.HALF ) T = 1 - T
      IF ( Q .GT. C30+C31*T ) GO TO 34
C
C  Classical series.
C  (If Q has changed value, compute N3 and A3(j), j=0,1,...,N3)
      IF ( Q.NE.QPRV3A ) THEN
         QSQ = Q**2
         P = Q
         ATEMP = 1
         DO 31 J = 0,IADIM
            A3(J) = ATEMP
            N3 = J
            ATEMP = P*ATEMP
            P = QSQ*P
            IF ( 1-ATEMP .EQ. 1 ) GO TO 32
   31    CONTINUE
   32    QPRV3A= Q
      ENDIF
C  (Compute the series using Clenshaw recurrence)
      ALPHA = 2*COS(2*PI*T)
      W = 0
      V = 0
      U = A3(N3)
      DO 33 J = N3-1,0,-1
         W = V
         V = U
         U = A3(J) + ALPHA*V - W
   33 CONTINUE
      RESULT = U - W
      GO TO 50
C
C  Two-term cosh approximation.
   34    IF ( Q.EQ.1 ) THEN
            IF ( T.EQ.0 ) THEN
               WRITE(ERRTXT,103) K,X,Q
               CALL MTLPRT(NAME,'C349.3',ERRTXT)
            ENDIF
            RESULT = 0
            GO TO 50
         ENDIF
C
C  (If Q has changed value, compute R, B3 and COEFF)
         IF ( Q.NE.QPRV3B ) THEN
            ABSLOG = ABS(LOG(Q))
            R = PISQ/ABSLOG
            IF ( -R.GE.ELOW ) B3 = 2*EXP(-R)
            COEFF = SQRT(PI/ABSLOG)
            QPRV3B = Q
         ENDIF
C  (Compute function)
      E = - R*(T**2)
      RESULT = COEFF*EXP(E)
      IF ( RESULT.EQ.0 ) THEN
         RESULT = EXP( E + LOG(COEFF) )
      ENDIF
      E = -R*(1-2*T)
      IF ( E.GT.ELOW ) THEN
         ETEST = -R*(1+2*T)
         IF ( ETEST.GE.ELOW ) THEN
            RESULT = RESULT*( 1 + B3*COSH(2*R*T) )
         ELSE
            RESULT = RESULT*( 1 + EXP(E) )
         ENDIF
      ENDIF
C
C
C--Terminate.
#if !defined(CERNLIB_DOUBLE)
   50 RTHETA = RESULT
#endif
#if defined(CERNLIB_DOUBLE)
   50 DTHETA = RESULT
#endif
      RETURN
C----------------------------------------------------------------------
C
C--Formats for error messages.
  101 FORMAT('Q < 0 or Q > 1       ',
     +       ' K=',I2,5X,'X=',1P,E13.5,5X,'Q=',E13.5)
  102 FORMAT('Impermissible K value',
     $       ' K=',I2,5X,'X=',1P,E13.5,5X,'Q=',E13.5)
  103 FORMAT('Function is infinite ',
     $       ' K=',I2,5X,'X=',1P,E13.5,5X,'Q=',E13.5)
      END
