*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION DSPNB1(K,M,I,NDER,X,T,NERR)

#include "gen/imp64.inc"
      DIMENSION T(*),B(27),D(27)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSPNB1')

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   DSPNB1 COMPUTES FUNCTION VALUES, VALUES OF DERIVATIVES, AND THE
*   VALUE OF INTEGRAL, RESPECTIVELY, OF NORMALIZED B-SPLINES
*                      B(I,K)(X)
*   OF DEGREE  K  ( 0<= K <= 25 )  WITH INDEX  I  ( 1 <= I <= M-K-1 )
*   OVER A SET OF SPLINE-KNOTS
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
*   PARAMETERS:
*
*   K     (INTEGER) DEGREE (= ORDER - 1) OF B-SPLINES.
*   M     (INTEGER) NUMBER OF KNOTS FOR THE B-SPLINES.
*   NDER  (INTEGER) ON ENTRY, NDER MUST CONTAIN AN INTEGER VALUE .GE. -1
*         = -1: DSPNB1 COMPUTES THE INTEGRAL OF  B(I,K)(TAU)  OVER THE
*               RANGE  TAU <= X.
*         =  0: DSPNB1 COMPUTES THE FUNCTION VALUE  B(I,K)(X)  FOR
*               FOR THE SPECIFIED VALUES OF I, K, AND X.
*         >= 1: DSPNB1 COMPUTES THE VALUE OF THE NDER-TH DERIVATIVE OF
*               B(I,K)(X)  FOR THE SPECIFIED VALUES OF I, K, AND X
*               (IF  NDER > K  ZERO RETURNS).
*   X     (DOUBLE PRECISION) ON ENTRY, X MUST CONTAIN THE VALUE OF THE
*         INDEPENDENT VARIABLE X OF  B(I,K)(X)
*   T     (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER M CONTAINING THE
*         KNOTS, ON ENTRY.
*   I     (INTEGER) INDEX OF THE B-SPLINE  B(I,K)(X)
*   NERR  (INTEGER) ERROR INDICATOR. ON EXIT:
*         = 0: NO ERROR DETECTED
*         = 1: AT LEAST ONE OF THE CONSTANTS K , M , I , NDER IS ILLEGAL
*
*   USAGE:
*
*   THE FUNCTION-CALL
*         A = DSPNB1(K,M,I,NDER,X,T,NERR)
*   RETURNS
*       - THE VALUE OF THE INTEGRAL           (NDER = -1) OR
*       - THE FUNCTION VALUE                  (NDER = 0 ) OR
*       - THE VALUE OF THE NDER-TH DERIVATIVE (NDER > 0 )
*   OF THE NORMALIZED B-SPLINE  B(I,K)(X) OF DEGREE K WITH INDEX I AT X.
*
*   ERROR MESSAGES:
*
*   IF ONE OF THE FOLLOWING RELATION IS SATISFIED BY THE CHOSEN INPUT-
*   PARAMETERS THE PROGRAM RETURNS, AND AN ERROR MESSAGE IS PRINTED:
*     K < 0      OR    K > 25    OR
*     M < 2*K+2  OR
*     NDER < -1  OR
*     I < 1      OR    I > M-K-1
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
      ELSEIF(I . LT. 1 .OR. I .GT. M-K-1) THEN
       WRITE(ERRTXT,101) 'I',I
       CALL MTLPRT(NAME,'E210.3',ERRTXT)
      ELSEIF(NDER .LT. -1) THEN
       WRITE(ERRTXT,101) 'NDER',NDER
       CALL MTLPRT(NAME,'E210.5',ERRTXT)
      ELSE

       NERR=0
       DSPNB1=Z0
       IF(     X .LT. T(I)
     +    .OR. X .GT. T(I+K+1) .AND. NDER .GE. 0
     +    .OR. K .LT. NDER                       ) RETURN

       IF(NDER .EQ. -1) THEN
        IF(X .GE. T(I+K+1)) THEN
         R=(T(I+K+1)-T(I))/(K+1)
        ELSE
         KK=LKKSPL(X,T(I),MIN(2*(K+1),M-K-I))+I-1
         IF(K .EQ. 0) THEN
          R=X-T(KK-1)
         ELSE
          CALL DVSET(K+1,Z0,B(1),B(2))
          B(KK-I)=1/(T(KK)-T(KK-1))
          DO 10 L=1,K
          DO 10 J=MAX(1,KK-I-L),MIN(K+1-L,KK-I)
          DIF=T(I+J+L)-T(I+J-1)
          B0=Z0
          IF(DIF .NE. 0) B0=((X-T(I+J-1))*B(J)+(T(I+J+L)-X)*B(J+1))/DIF
   10     B(J)=B0
          S=Z0
          DO 20 L=1,KK-I
   20     S=S+(X-T(I+L-1))*B(L)
          R=S*(T(I+K+1)-T(I))/(K+1)
         ENDIF
        ENDIF
        DSPNB1=R
        RETURN
       ENDIF

       IF(K .EQ. 0) THEN
        R=Z1
       ELSE
        KK=LKKSPL(X,T(I),MIN(2*(K+1),M-K-I))+I-1
        I0=I+K+2-KK
        IF(I0 .EQ. 0) THEN
         R=Z0
        ELSE
         E1=X-T(KK-1)
         B(1)=Z1
         DO 30 J=2,K-NDER+1
   30    B(J)=E1*B(J-1)/(T(KK-2+J)-T(KK-1))
         IF(KK .NE. I+1 .OR. NDER .NE. 0) THEN
          E2=T(KK)-X
          DO 40 J=1,K-NDER
          E3=X-T(KK-1-J)
          B(1)=E2*B(1)/(T(KK)-T(KK-J))
          DO 40 L=2,K-NDER+1-J
   40     B(L)=E3*B(L-1)/(T(KK-2+L)-T(KK-1-J))+
     +         (T(KK-1+L)-X)*B(L)/(T(KK-1+L)-T(KK-J))
         ENDIF
         IF(NDER .EQ. 0) THEN
          R=B(I0)
         ELSE
          CALL DVSET(K+2,Z0,D(1),D(2))
          D(I+K+2-KK)=1
          DO 50 J=1,NDER
          A=K-J+1
          DO 50 L=1,K-J+2
          DIF=T(L+KK-1)-T(L+KK-K-2+J)
          D0=Z0
          IF(DIF .NE. 0) D0=A*(D(L+1)-D(L))/DIF
   50     D(L)=D0
          R=DVMPY(K-NDER+1,B(1),B(2),D(1),D(2))
         ENDIF
        ENDIF
       ENDIF
       DSPNB1=R
      ENDIF
      RETURN

  101 FORMAT(1X,A5,' =',I6,'   NOT IN RANGE')
      END


