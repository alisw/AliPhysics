*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION DSPNB2(KX,KY,MX,MY,I,J,NDERX,NDERY,X,Y,TX,TY,NERR)

#include "gen/imp64.inc"
      DIMENSION TX(*),TY(*),B(27),D(27)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSPNB2')

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   DSPNB2 COMPUTES FUNCTION VALUES, VALUES OF DERIVATIVES, AND THE
*   VALUE OF INTEGRAL, RESPECTIVELY, OF NORMALIZED TWO-DIMENSIONAL
*   B-SPLINES IN THE FORM OF A PRODUCT OF TWO ONE-DIMENSIONAL NORMALIZED
*   B-SPLINES
*              B(I,J)(X,Y) = BX(I,KX)(X) * BY(J,KY)(Y)
*   OF DEGREE  KX AND KY ( 0 <= KX , KY <= 25 )   WITH INDICES  I , J
*   ( 1 <= I <= MX-KX-1 , 1 <= J <= MY-KY-1 ) OVER TWO SETS OF SPLINE-
*   KNOTS
*       TX(1),TX(2),...,TX(MX)   ( MX >= 2*KX+2 )
*       TY(1),TY(2),...,TY(MY)   ( MY >= 2*KY+2 ) ,
*   RESPECTIVELY.
*   FOR FURTHER DETAILS TO THE ONE-DIMENSIONAL NORMALIZED B-SPLINES SEE
*   THE COMMENTS TO  DSPNB1.
*
*   PARAMETERS:
*
*   KX    (INTEGER) DEGREE OF ONE-DIMENSIONAL B-SPLINES IN X-DIRECTION
*         OVER THE SET OF KNOTS  TX.
*   KY    (INTEGER) DEGREE OF ONE-DIMENSIONAL B-SPLINES IN Y-DIRECTION
*         OVER THE SET OF KNOTS  TY.
*   MX    (INTEGER) NUMBER OF KNOTS FOR THE B-SPLINES IN X-DIRECTION.
*   MY    (INTEGER) NUMBER OF KNOTS FOR THE B-SPLINES IN Y-DIRECTION.
*   NDERX (INTEGER) ON ENTRY, NDERX MUST CONTAIN AN INTEGER VALUE >= -1.
*         = -1: DSPNB2 COMPUTES THE INTEGRAL OF  B(I,J)(TAU,Y)  OVER THE
*               RANGE  TAU <= X.
*         =  0: DSPNB2 COMPUTES THE FUNCTION VALUE  B(I,J)(X,Y)  FOR
*               FOR THE SPECIFIED VALUES OF I, J, X, AND Y.
*         >= 1: DSPNB2 COMPUTES THE VALUE OF THE NDERX-TH PARTIAL
*               DERIVATIVE OF  B(I,J)(X,Y)  WITH RESPECT TO X
*               FOR THE SPECIFIED VALUES OF I, J, X, AND Y.
*               (IF  NDERX > KX  ZERO RETURNS).
*   NDERY (INTEGER) ON ENTRY, NDERY MUST CONTAIN AN INTEGER VALUE >= -1.
*         THE MEANING OF  NDERY  IS THE SAME AS THAT OF THE PARAMETER
*         NDERX WITH RESPECT TO Y-DIRECTION INSTEAD OF X-DIRECTION.
*   X     (DOUBLE PRECISION) ON ENTRY, X MUST CONTAIN THE VALUE OF THE
*         INDEPENDENT VARIABLE X OF  B(I,J)(X,Y)
*   Y     (DOUBLE PRECISION) ON ENTRY, Y MUST CONTAIN THE VALUE OF THE
*         INDEPENDENT VARIABLE Y OF  B(I,J)(X,Y)
*   TX    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER MX CONTAINING THE
*         KNOTS IN X-DIRECTION, ON ENTRY.
*   TY    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER MY CONTAINING THE
*         KNOTS IN Y-DIRECTION, ON ENTRY.
*   I     (INTEGER) INDEX OF THE B-SPLINE  B(I,J)(X,Y)
*   J     (INTEGER) INDEX OF THE B-SPLINE  B(I,J)(X,Y)
*   NERR  (INTEGER) ERROR INDICATOR. ON EXIT:
*         = 0: NO ERROR DETECTED
*         = 1: AT LEAST ONE OF THE CONSTANTS KX , KY , MX , MY , I , J ,
*              NDERX , NDERY  IS ILLEGAL
*
*   USAGE:
*
*   THE FUNCTION-CALL
*         A = DSPNB2(KX,KY,MX,MY,I,J,NDERX,NDERY,X,Y,TX,TY,NERR)
*   RETURNS
*       - THE VALUE OF THE INTEGRAL OR
*       - THE FUNCTION VALUE        OR
*       - THE VALUE OF THE NDERX-TH AND NDERY-TH PARTIAL DERIVATIVE
*   OF THE NORMALIZED B-SPLINE  B(I,J)(X,Y)  WITH INDICES I,J AT (X,Y).
*
*   ERROR MESSAGES:
*
*   IF ONE OF THE FOLLOWING RELATION IS SATISFIED BY THE CHOSEN INPUT-
*   PARAMETERS THE PROGRAM RETURNS, AND AN ERROR MESSAGE IS PRINTED:
*     KX < 0       OR    KX > 25      OR    KY < 0   OR   KY > 25   OR
*     MX < 2*KX+2  OR    MY < 2*KY+2  OR
*     NDERX < -1   OR    NDERY < -1   OR
*     I < 1        OR    I > MX-KX-1  OR    J < 1    OR   J > MY-KY-1 .
*
************************************************************************

      PARAMETER (Z0 = 0 , Z1 = 1)

      NERR=1
      IF(KX .LT. 0 .OR. KX .GT. 25) THEN
       WRITE(ERRTXT,101) 'KX',KX
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(KY .LT. 0 .OR. KY .GT. 25) THEN
       WRITE(ERRTXT,101) 'KY',KY
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(MX .LT. 2*KX+2) THEN
       WRITE(ERRTXT,101) 'MX',MX
       CALL MTLPRT(NAME,'E210.2',ERRTXT)
      ELSEIF(MY .LT. 2*KY+2) THEN
       WRITE(ERRTXT,101) 'MY',MY
       CALL MTLPRT(NAME,'E210.2',ERRTXT)
      ELSEIF(I .LT. 1 .OR. I .GT. MX-KX-1) THEN
       WRITE(ERRTXT,101) 'I',I
       CALL MTLPRT(NAME,'E210.3',ERRTXT)
      ELSEIF(J .LT. 1 .OR. J .GT. MY-KY-1) THEN
       WRITE(ERRTXT,101) 'J',J
       CALL MTLPRT(NAME,'E210.3',ERRTXT)
      ELSEIF(NDERX .LT. -1) THEN
       WRITE(ERRTXT,101) 'NDERX',NDERX
       CALL MTLPRT(NAME,'E210.5',ERRTXT)
      ELSEIF(NDERY .LT. -1) THEN
       WRITE(ERRTXT,101) 'NDERY',NDERY
       CALL MTLPRT(NAME,'E210.5',ERRTXT)
      ELSEIF(NDERX .LT. 0 .AND. NDERY .GT. 0  .OR.
     +       NDERY .LT. 0 .AND. NDERX .GT. 0  ) THEN
       WRITE(ERRTXT,102) 'NDERX',NDERX,'NDERY',NDERY
       CALL MTLPRT(NAME,'E210.6',ERRTXT)
      ELSE

       NERR=0
       DSPNB2=Z0

       IF(     X  .LT. TX(I)      .OR.  Y  .LT. TY(J)
     +    .OR. X  .GT. TX(I+KX+1) .AND. NDERX .GE. 0
     +    .OR. Y  .GT. TY(J+KY+1) .AND. NDERY .GE. 0
     +    .OR. KX .LT. NDERX      .OR.  KY .LT. NDERY  ) RETURN
*
*   COMPUTING  SPLNBX
*
       IF(NDERX .EQ. -1) THEN
        IF(X .GE. TX(I+KX+1)) THEN
         SPLNBX=(TX(I+KX+1)-TX(I))/(KX+1)
        ELSE
         KK=LKKSPL(X,TX(I),MIN(2*(KX+1),MX-KX-I))+I-1
         IF(KX .EQ. 0) THEN
          SPLNBX=X-TX(KK-1)
         ELSE
          CALL DVSET(KX+1,Z0,B(1),B(2))
          B(KK-I)=1/(TX(KK)-TX(KK-1))
          DO 10 L=1,KX
          DO 10 K=MAX(1,KK-I-L),MIN(KX+1-L,KK-I)
          DIF=TX(I+K+L)-TX(I+K-1)
          B0=Z0
          IF(DIF .NE. 0)
     +     B0=((X-TX(I+K-1))*B(K)+(TX(I+K+L)-X)*B(K+1))/DIF
   10     B(K)=B0
          S=Z0
          DO 20 L=1,KK-I
   20     S=S+(X-TX(I+L-1))*B(L)
          SPLNBX=S*(TX(I+KX+1)-TX(I))/(KX+1)
         ENDIF
        ENDIF
        GO TO 60
       ENDIF

       IF(KX .EQ. 0) THEN
        SPLNBX=Z1
       ELSE
        KK=LKKSPL(X,TX(I),MIN(2*(KX+1),MX-KX-I))+I-1
        I0=I+KX+2-KK
        IF(I0 .EQ. 0) THEN
         DSPNB2=Z0
         RETURN
        ELSE
         E1=X-TX(KK-1)
         B(1)=Z1
         DO 30 K=2,KX-NDERX+1
   30    B(K)=E1*B(K-1)/(TX(KK-2+K)-TX(KK-1))
         IF(KK .NE. I+1 .OR. NDERX .NE. 0) THEN
          E2=TX(KK)-X
          DO 40 K=1,KX-NDERX
          E3=X-TX(KK-1-K)
          B(1)=E2*B(1)/(TX(KK)-TX(KK-K))
          DO 40 L=2,KX-NDERX+1-K
   40     B(L)=E3*B(L-1)/(TX(KK-2+L)-TX(KK-1-K))+
     +         (TX(KK-1+L)-X)*B(L)/(TX(KK-1+L)-TX(KK-K))
         ENDIF
         IF(NDERX .EQ. 0) THEN
          SPLNBX=B(I0)
         ELSE
          CALL DVSET(KX+2,Z0,D(1),D(2))
          D(I+KX+2-KK)=1
          DO 50 K=1,NDERX
          A=KX-K+1
          DO 50 L=1,KX-K+2
          DIF=TX(L+KK-1)-TX(L+KK-KX-2+K)
          D0=Z0
          IF(DIF .NE. 0) D0=A*(D(L+1)-D(L))/DIF
   50     D(L)=D0
          SPLNBX=DVMPY(KX-NDERX+1,B(1),B(2),D(1),D(2))
         ENDIF
        ENDIF
       ENDIF
*
*   COMPUTING  SPLNBY
*
   60  IF(NDERY .EQ. -1) THEN
        IF(Y .GE. TY(J+KY+1)) THEN
         SPLNBY=(TY(J+KY+1)-TY(J))/(KY+1)
        ELSE
         KK=LKKSPL(Y,TY(J),MIN(2*(KY+1),MY-KY-J))+J-1
         IF(KY .EQ. 0) THEN
          SPLNBY=Y-TY(KK-1)
         ELSE
          CALL DVSET(KY+1,Z0,B(1),B(2))
          B(KK-J)=1/(TY(KK)-TY(KK-1))
          DO 70 L=1,KY
          DO 70 K=MAX(1,KK-J-L),MIN(KY+1-L,KK-J)
          DIF=TY(J+K+L)-TY(J+K-1)
          B0=Z0
          IF(DIF .NE. 0)
     +     B0=((Y-TY(J+K-1))*B(K)+(TY(J+K+L)-Y)*B(K+1))/DIF
   70     B(K)=B0
          S=Z0
          DO 80 L=1,KK-J
   80     S=S+(Y-TY(J+L-1))*B(L)
          SPLNBY=S*(TY(J+KY+1)-TY(J))/(KY+1)
         ENDIF
        ENDIF
        GO TO 120
       ENDIF

       IF(KY .EQ. 0) THEN
        SPLNBY=Z1
       ELSE
        KK=LKKSPL(Y,TY(J),MIN(2*(KY+1),MY-KY-J))+J-1
        I0=J+KY+2-KK
        IF(I0 .EQ. 0) THEN
         DSPNB2=Z0
         RETURN
        ELSE
         E1=Y-TY(KK-1)
         B(1)=Z1
         DO 90 K=2,KY-NDERY+1
   90    B(K)=E1*B(K-1)/(TY(KK-2+K)-TY(KK-1))
         IF(KK .NE. J+1 .OR. NDERY .NE. 0) THEN
          E2=TY(KK)-Y
          DO 100 K=1,KY-NDERY
          E3=Y-TY(KK-1-K)
          B(1)=E2*B(1)/(TY(KK)-TY(KK-K))
          DO 100 L=2,KY-NDERY+1-K
  100     B(L)=E3*B(L-1)/(TY(KK-2+L)-TY(KK-1-K))+
     +         (TY(KK-1+L)-Y)*B(L)/(TY(KK-1+L)-TY(KK-K))
         ENDIF
         IF(NDERY .EQ. 0) THEN
          SPLNBY=B(I0)
         ELSE
          CALL DVSET(KY+2,Z0,D(1),D(2))
          D(J+KY+2-KK)=1
          DO 110 K=1,NDERY
          A=KY-K+1
          DO 110 L=1,KY-K+2
          DIF=TY(L+KK-1)-TY(L+KK-KY-2+K)
          D0=Z0
          IF(DIF .NE. 0) D0=A*(D(L+1)-D(L))/DIF
  110     D(L)=D0
          SPLNBY=DVMPY(KY-NDERY+1,B(1),B(2),D(1),D(2))
         ENDIF
        ENDIF
       ENDIF

  120  DSPNB2=SPLNBX*SPLNBY

      ENDIF
      RETURN

  101 FORMAT(1X,A5,' =',I6,'   NOT IN RANGE')
  102 FORMAT(1X,A5,' =',I6,A7,' =',I6,'   INCONSISTENT')
      END


