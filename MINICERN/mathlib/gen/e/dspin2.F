*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:25  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DSPIN2(KX,KY,NX,NY,XI,YI,ZI,NDIMZ,KNOT,TX,TY,C,NDIMC,
     +                  W,IW,NERR)

#include "gen/imp64.inc"
      DIMENSION XI(*),YI(*),ZI(NDIMZ,*),TX(*),TY(*)
      DIMENSION W(*),IW(*),C(NDIMC,*)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DSPIN2')


************************************************************************
*   NORBAS, VERSION: 10.02.1993
************************************************************************
*
*   DSPIN2 COMPUTES THE COEFFICIENTS
*          C(I,J)   (I=1,...,NX , J=1,...,NY)
*   OF A TWO-DIMENSIONAL POLYNOMIAL INTERPOLATION SPLINE  Z = S(X,Y)  IN
*   REPRESENTATION OF NORMALIZED TWO-DIMENSIONAL B-SPLINES  B(I,J)(X,Y)
*
*     S(X,Y) = SUMME(I=1,...,NX)
*              SUMME(J=1,...,NY)  C(I,J) * B(I,J)(X,Y)
*
*   TO A USER SUPPLIED DATA SET
*
*     (XI(I),YI(J),ZI(I,J))       (I=1,...,NX , J=1,...,NY)
*
*   OF A FUNCTION  Z = F(X,Y) , I.E.
*
*     S(XI(I),YI(J)) = Z(I,J)     (I=1,...,NX , J=1,...,NY) .
*
*   THE TWO-DIMENSIONAL B-SPLINES  B(I,J)(X,Y)  ARE THE PRODUCT OF TWO
*   ONE-DIMENSIONAL B-SPLINES  BX , BY
*          B(I,J)(X,Y) = BX(I,KX)(X) * BY(J,KY)(Y)
*   OF DEGREE  KX AND KY ( 0 <= KX , KY <= 25 )   WITH INDICES  I , J
*   ( 1 <= I <= NX , 1 <= J <= NY ) OVER TWO SETS OF SPLINE-KNOTS
*       TX(1),TX(2),...,TX(MX)     ( MX = NX+KX+1 )
*       TY(1),TY(2),...,TY(MY)     ( MY = NY+KY+1 ) ,
*   RESPECTIVELY.
*   FOR FURTHER DETAILS TO THE ONE- AND TWO-DIMENSIONAL NORMALIZED
*   B-SPLINES SEE THE COMMENTS TO  DSPNB1  AND  DSPNB2.
*
*   PARAMETERS:
*
*   NX    (INTEGER) NUMBER OF INTERPOLATION POINTS IN X-DIRECTION :
*         XI(I) ,  I=1,...,NX .
*   NY    (INTEGER) NUMBER OF INTERPOLATION POINTS IN Y-DIRECTION :
*         YI(J) ,  J=1,...,NY .
*   KX    (INTEGER) DEGREE OF ONE-DIMENSIONAL B-SPLINES IN X-DIRECTION
*         OVER THE SET OF KNOTS  TX,
*         WITH  KX <= NX-1  AND  0 <= KX <= 25 .
*   KY    (INTEGER) DEGREE OF ONE-DIMENSIONAL B-SPLINES IN Y-DIRECTION
*         OVER THE SET OF KNOTS  TY,
*         WITH  KY <= NY-1  AND  0 <= KY <= 25 .
*   NDIMC (INTEGER) DECLARED FIRST DIMENSION OF ARRAY  C  IN THE
*         CALLING PROGRAM, WITH  NDIMC >= NX .
*   NDIMZ (INTEGER) DECLARED FIRST DIMENSION OF ARRAY  ZI  IN THE
*         CALLING PROGRAM, WITH  NDIMZ >= NX .
*   XI    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER NX .
*         XI MUST CONTAIN THE INTERPOLATION POINTS IN X-DIRECTION IN
*         ASCENDING ORDER, ON ENTRY.
*   YI    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER NY .
*         YI MUST CONTAIN THE INTERPOLATION POINTS IN Y-DIRECTION IN
*         ASCENDING ORDER, ON ENTRY.
*   ZI    (DOUBLE PRECISION) ARRAY OF ORDER  (NDIMZ , >= NY) .
*         ON ENTRY ZI MUST CONTAIN THE GIVEN FUNCTION VALUES
*           Z(I,J)  AT THE INTERPOLATION POINTS  (X(I),Y(J))
*                            ( I=1,...,NX , J=1,...,NY ).
*   KNOT  (INTEGER) PARAMETER FOR STEERING THE CHOICE OF KNOTS.
*         ON ENTRY:
*         = 1 : KNOTS ARE COMPUTED BY  DSPIN2  IN THE FOLLOWING WAY:
*               TX(J) = XI(1) ,                J = 1,...,KX+1
*               TX(J) = XI(1)+(J-KX-1)*(XI(NX)-XI(1)) ,
*                                              J = KX+2,...,NX
*               TX(N+J) = XI(NX) ,             J = 1,...,KX+1
*         = 2 : KNOTS ARE COMPUTED BY  DSPIN2  IN THE FOLLOWING WAY:
*               TX(J) = XI(1) ,                J = 1,...,KX+1
*               TX(J) = (XI(J-KX-1)+XI(J))/2 , J = KX+2,...,NX
*               TX(N+J) = XI(NX) ,             J = 1,...,KX+1
*         OTHERWISE KNOTS ARE USER SUPPLIED. RECOMMENDED CHOICE :
*               T(1) <= ... <= T(K+1) <= XI(1)
*               XI(1) < T(K+2) < ... < T(N) < XI(N)
*               XI(N) <= T(N+1) <= ... <= T(N+K+1)
*         OTHERWISE KNOTS ARE USER SUPPLIED. RECOMMENDED CHOICE :
*               TX(1) <= ... <= TX(KX+1) <= XI(1)
*               XI(1) < TX(KX+2) < ... < TX(NX) < XI(NX)
*               XI(NX) <= TX(NX+1) <= ... <= TX(NX+KX+1)
*         IN ALL CASES THE SAME CHOICE IS USED FOR KNOTS  TY  IN
*         Y-DIRECTION .
*   TX    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER MX.
*         IF THE INPUT VALUE OF THE PARAMETER KNOT IS  1 OR 2  THE
*         KNOTS ARE COMPUTED BY  DSPIN2  AND THEY ARE GIVEN IN THE
*         ARRAY TX, ON EXIT.
*         IN THE OTHER CASES THE ARRAY  TX  MUST CONTAIN THE USER
*         SUPPLIED KNOTS IN X-DIRECTION, ON ENTRY.
*   TY    (DOUBLE PRECISION) ARRAY OF AT LEAST ORDER MY.
*         IF THE INPUT VALUE OF THE PARAMETER KNOT IS  1 OR 2  THE
*         KNOTS ARE COMPUTED BY  DSPIN2  AND THEY ARE GIVEN IN THE
*         ARRAY TY, ON EXIT.
*         IN THE OTHER CASES THE ARRAY  TY  MUST CONTAIN THE USER
*         SUPPLIED KNOTS IN Y-DIRECTION, ON ENTRY.
*   W     (DOUBLE PRECISION) WORKING ARRAY OF AT LEAST ORDER (L+1)*NINT,
*         WITH  NINT=NX*NY ,  K=KX*NY ,  L=3*K+1 .
*   IW    (INTEGER) WORKING ARRAY OF AT LEAST ORDER  NX*NY.
*   C     (DOUBLE PRECISION) ARRAY OF ORDER (NDIMC, >= NY).
*         ON EXIT  C(I,J)  CONTAINS THE (I,J)-TH COEFFICIENT OF THE
*         TWO-DIMENSIONAL B-SPLINE REPRESENTATION OF  S(X,Y) .
*   NERR  (INTEGER) ERROR INDICATOR. ON EXIT:
*         = 0: NO ERROR DETECTED
*         = 1: AT LEAST ONE OF THE CONSTANTS  KX , KY , NX , NY
*              IS ILLEGAL
*         = 2: THE LAPACK ROUTINES  DGBTRF , DGBTRS  COULD NOT SOLVE
*              THE LINEAR SYSTEM OF EQUATIONS WITH BAND-MATRIX FOR
*              COMPUTING C(1,1),...,C(NX,NY) . IT INDICATES THAT
*              A SOLUTION OF THE INTERPOLATION PROBLEM DOES NOT EXIST.
*              (ESPECIALLY, THE EXISTENCE OF A SOLUTION DEPENDS ON THE
*              SET OF KNOTS!)
*
*   ERROR MESSAGES:
*
*   IF ONE OF THE FOLLOWING RELATION IS SATISFIED BY THE CHOSEN INPUT-
*   PARAMETERS THE PROGRAM RETURNS, AND AN ERROR MESSAGE IS PRINTED:
*     KX < 0      OR    KX > 25    OR
*     KY < 0      OR    KY > 25    OR
*     NX < KX+1   OR    NY < KY+1 .
*
************************************************************************


      PARAMETER (Z1 = 1 , Z2 = 2 , HALF = Z1/Z2)

      NERR=1
      IF(KX .LT. 0 .OR. KX .GT. 25) THEN
       WRITE(ERRTXT,101) 'KX',KX
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(KY .LT. 0 .OR. KY .GT. 25) THEN
       WRITE(ERRTXT,101) 'KY',KY
       CALL MTLPRT(NAME,'E210.1',ERRTXT)
      ELSEIF(NX .LT. KX+1) THEN
       WRITE(ERRTXT,101) 'NX',NX
       CALL MTLPRT(NAME,'E210.4',ERRTXT)
      ELSEIF(NY .LT. KY+1) THEN
       WRITE(ERRTXT,101) 'NY',NY
       CALL MTLPRT(NAME,'E210.4',ERRTXT)
      ELSE

       K=KX*NY
       L=3*K+1
       NINT=NX*NY
       MX=NX+KX+1
       MY=NY+KY+1
*
*   COMPUTE KNOTS FROM INTERPOLATION POINTS ( IF KNOT=1 OR 2 )
*
       IF (KNOT .EQ. 1) THEN
        CALL DSPKN2(KX,KY,MX,MY,XI(1),XI(NX),YI(1),YI(NY),TX,TY,NERR)
       ELSEIF (KNOT .EQ. 2) THEN
        DO 10 I=1,KX+1
        TX(I)=XI(1)
   10   TX(NX+I)=XI(NX)
        DO 20 I=KX+2,NX
   20   TX(I)=HALF*(XI(I-KX-1)+XI(I))
        DO 30 J=1,KY+1
        TY(J)=YI(1)
   30   TY(NY+J)=YI(NY)
        DO 40 J=KY+2,NY
   40   TY(J)=HALF*(YI(J-KY-1)+YI(J))
       ENDIF
*
*   COMPUTE BAND-MATRIX  W
*
*   NUMBER OF SUB-, SUPER-DIAGONALS: KL=KU=(KX-1)*NY+KY-1<=KX*NY
*
       DO 90 IC=1,NX
       X=XI(IC)
       DO 90 JC=1,NY
       Y=YI(JC)
       IJ=(IC-1)*NX+JC
       DO 90 I=1,NX
       DO 90 J=1,NY
       JZ=(I-1)*NX+J
       IZ=IJ+2*K+1-JZ
       IF (K+1 .LE. IZ  .AND.  IZ .LE. 3*K+1)
     +    W((JZ-1)*L+IZ)=DSPNB2(KX,KY,MX,MY,I,J,0,0,X,Y,TX,TY,NERR)
   90  CONTINUE
*
*   SOLVE SYSTEM OF EQUATIONS FOR COMPUTING  C
*
       DO 110 J=1,NINT
  110  IW(J)=J
       NERR=2
       CALL DGBTRF(NINT,NINT,K,K,W,L,IW,INFO)
       IF(INFO .NE. 0) RETURN
       LW=L*NINT
       DO 120 I=1,NX
       DO 120 J=1,NY
  120  W(LW+(I-1)*NY+J)=ZI(I,J)
       CALL DGBTRS('N',NINT,K,K,1,W,L,IW,W(LW+1),NINT,INFO)
       IF(INFO .NE. 0) RETURN
       DO 130 I=1,NX
       DO 130 J=1,NY
  130  C(I,J)=W(LW+(I-1)*NY+J)
       NERR=0
      ENDIF

      RETURN

  101 FORMAT(1X,A5,' =',I6,'   NOT IN RANGE')
      END



