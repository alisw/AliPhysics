*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION DSPPS1(K,M,NDER,X,T,C,NERR)

#include "gen/imp64.inc"
      DIMENSION T(*),C(*),B(27)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSPPS1')

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   DSPPS1 COMPUTES FUNCTION VALUES, VALUES OF DERIVATIVES, AND THE
*   VALUE OF INTEGRAL, RESPECTIVELY, OF A POLYNOMIAL SPLINE  S(X)  IN
*   B-SPLINE REPRESENTATION
*
*          S(X) =  SUMME(I=1,...,M-K-1)  C(I) * B(I,K)(X) .
*
*   THE FUNCTIONS  B(I,K)(X)  ARE NORMALIZED B-SPLINES OF DEGREE  K
*   (0<= K <= 25)  WITH INDEX  I  (1 <= I <= M-K-1) OVER A SET OF
*   SPLINE-KNOTS
*              T(1),T(2), ... ,T(M)    ( M >= 2*K+2 )
*   (KNOTS IN ASCENDING ORDER, WITH MULTIPLICITIES NOT GREATER
*   THAN  K+1).
*
*   THE FUNCTION VALUE OF THE NORMALIZED B-SPLINE  B(I,K)(X)  IS
*   IDENTICALLY ZERO OUTSIDE THE INTERVAL  T(I) <= X < T(I+K+1).
*
*   THE NORMALIZATION OF  N(X) = B(I,K)(X)  IS SUCH THAT THE INTGRAL OF
*   N(X)  OVER THE WHOLE X-RANGE EQUALS
*                  ( T(I+K+1) - T(I) ) / (K+1)  .
*
*   C(1),...,C(M-K-NDER-1)  MUST CONTAIN THE COEFFICIENTS OF THE
*   POLYNOMIAL SPLINE  S(X)  OR ITS DERIVATIVE, REPRESENTED BY
*   NORMALIZED B-SPLINES OF DEGREE  K-NDER .
*
*   FOR TRANSFORMATION THE COEFFICIENTS OF THE POLYNOMIAL SPLINE  S(X)
*   TO THE CORRESPONDING COEFFICIENTS OF THE NDER-TH DERIVATIVE OF
*   S(X)  THE ROUTINE  DSPCD1  MAY BE USED.
*
*   ESPECIALLY FOR COMPUTING THE COEFFICIENTS  C(1),...,C(M-K-1)  OF
*   THE POLYNOMIAL VARIATION DIMINISHING SPLINE APPROXIMATION OF A USER
*   SUPPLIED FUNCTION  F(X)  THE ROUTINE  DSPVD1  MAY BE USED.
*
*   PARAMETERS:
*
*   K     (INTEGER) DEGREE (= ORDER - 1) OF B-SPLINES.
*   M     (INTEGER) NUMBER OF KNOTS FOR THE B-SPLINES.
*   NDER  (INTEGER) ON ENTRY, NDER MUST CONTAIN AN INTEGER VALUE .GE. -1
*         = -1: DSPPS1 COMPUTES THE INTEGRAL OF  S(TAU)  OVER THE
*               RANGE  TAU <= X.
*         =  0: DSPPS1 COMPUTES THE FUNCTION VALUE OF THE POLYNOMIAL
*               SPLINE  S(X)  AT  X.
*         >= 1: DSPPS1 COMPUTES THE VALUE OF THE NDER-TH DERIVATIVE OF
*               S(X)  AT X  (IF  NDER > K  ZERO RETURNS).
*   X     (DOUBLE PRECISION) ON ENTRY, X MUST CONTAIN THE VALUE OF THE
*         INDEPENDENT VARIABLE X OF  S(X).
*   T     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER M CONTAINING THE
*         KNOTS, ON ENTRY.
*   C     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER  M-K-NDER-1,
*         ON ENTRY C(J) CONTAINS THE J-TH COEFFICIENT OF THE B-SPLINE
*         REPRESENTATION OF S(X) OR OF ITS DERIVATIVE.
*   NERR  (INTEGER) ERROR INDICATOR. ON EXIT:
*         = 0: NO ERROR DETECTED
*         = 1: AT LEAST ONE OF THE CONSTANTS K , M , NDER IS ILLEGAL
*
*   USAGE:
*
*   THE FUNCTION-CALL
*         Y = DSPPS1(K,M,NDER,X,T,C,NERR)
*   RETURNS
*       - THE VALUE OF THE INTEGRAL           (NDER = -1) OR
*       - THE FUNCTION VALUE                  (NDER = 0 ) OR
*       - THE VALUE OF THE NDER-TH DERIVATIVE (NDER > 0 )
*   OF THE POLYNOMIAL SPLINE (IN B-SPLINE REPRESENTATION)
*          S(X) = SUMME(I=1,...,M-K-1)  C(I)*B(I,K)(X)
*   OF DEGREE K WITH THE SET OF KNOTS T(1),...,T(M)  AT X.
*
*   ERROR MESSAGES:
*
*   IF ONE OF THE FOLLOWING RELATION IS SATISFIED BY THE CHOSEN INPUT-
*   PARAMETERS THE PROGRAM RETURNS, AND AN ERROR MESSAGE IS PRINTED:
*     K < 0      OR    K > 25    OR
*     M < 2*K+2  OR
*     NDER < -1  .
*
************************************************************************

      PARAMETER (Z0 = 0 , Z1 = 1)

      NERR=1
      IF(K .LT. 0 .OR. K .GT. 25) THEN
       WRITE(ERRTXT,101) 'K',K
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(M .LT. 2*K+2) THEN
       WRITE(ERRTXT,101) 'M',M
       CALL MTLPRT(NAME,'E210.2',ERRTXT)
      ELSEIF(NDER .LT. -1) THEN
       WRITE(ERRTXT,101) 'NDER',NDER
       CALL MTLPRT(NAME,'E210.5',ERRTXT)
      ELSE

       NERR=0
       DSPPS1=Z0
       IF(     X .LT. T(K+1)
     +    .OR. X .GT. T(M-K) .AND. NDER .GE. 0
     +    .OR. K .LT. NDER                      ) RETURN

       IF(NDER .EQ. -1) THEN
        IF(X .GE. T(M-K)) THEN
         R=0
         DO 10 I=1,M-K-1
   10    R=R+C(I)*(T(I+K+1)-T(I))
         R=R/(K+1)
        ELSE
         KK=LKKSPL(X,T(K+1),M-2*K-1)+K
         R=0
         DO 20 I=1,KK-K-2
   20    R=R+C(I)*(T(I+K+1)-T(I))
         R=R/(K+1)
         IF(K .EQ. 0) THEN
          K1=MAX(1,KK-1)
          R=R+(X-T(K1))*C(K1)
         ELSE
          DO 50 I=MAX(1,KK-K-1),KK-1
          CALL DVSET(K+1,Z0,B(1),B(2))
          B(KK-I)=1/(T(KK)-T(KK-1))

          DO 30 L=1,K
          DO 30 J=MAX(1,KK-I-L),MIN(K+1-L,KK-I)
          DIF=T(I+J+L)-T(I+J-1)
          B0=Z0
          IF(DIF .NE. 0) B0=((X-T(I+J-1))*B(J)+(T(I+J+L)-X)*B(J+1))/DIF
   30     B(J)=B0
          S=Z0
          DO 40 L=1,KK-I
   40     S=S+(X-T(I+L-1))*B(L)
          S=S*(T(I+K+1)-T(I))/(K+1)
   50     R=R+C(I)*S
         ENDIF
        ENDIF
       ELSE

        KK=LKKSPL(X,T(K+1),M-2*K-1)+K
        E1=X-T(KK-1)
        B(1)=Z1
        DO 60 J=2,K-NDER+1
   60   B(J)=E1*B(J-1)/(T(KK-2+J)-T(KK-1))
        IF(KK .NE. 0  .OR.  NDER .NE. 0) THEN
         E2=T(KK)-X
         DO 70 J=1,K-NDER
         E3=X-T(KK-1-J)
         B(1)=E2*B(1)/(T(KK)-T(KK-J))
         DO 70 L=2,K-NDER+1-J
   70    B(L)=E3*B(L-1)/(T(KK-2+L)-T(KK-1-J))+
     +                  (T(KK-1+L)-X)*B(L)/(T(KK-1+L)-T(KK-J))
        ENDIF
        R=DVMPY(K-NDER+1,C(KK-1-K),C(KK-K),B(1),B(2))
       ENDIF
       DSPPS1=R
      ENDIF
      RETURN

  101 FORMAT(1X,A5,' =',I6,'   NOT IN RANGE')
      END



