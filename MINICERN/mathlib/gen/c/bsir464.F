*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:03  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION  BSIR4(X,NU)
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION DBSIR4(X,NU)
#include "gen/imp64.inc"
#endif
      CHARACTER*80 ERRTXT
      CHARACTER NAMEI*(*),NAMEK*(*),NAMEIE*(*),NAMEKE*(*)
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAMEI = 'BSIR4/DBSIR4', NAMEIE = 'EBSIR4/DEBIR4')
      PARAMETER (NAMEK = 'BSKR4/DBSKR4', NAMEKE = 'EBSKR4/DEBKR4')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAMEI = 'BSIR4/DBSIR4', NAMEIE = 'EBSIR4/DEBIR4')
      PARAMETER (NAMEK = 'BSKR4/DBSKR4', NAMEKE = 'EBSKR4/DEBKR4')
#endif
      LOGICAL LEX

      DIMENSION BC(0:23,2),CC(0:15,2),PP(-3:3),GG(-3:3)

      PARAMETER (EPS = 1D-14)
      PARAMETER (Z1 = 1, HF =Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (W2 = 1.41421 35623 73095 05D0)
      PARAMETER (G1 = 3.62560 99082 21908 31D0)
      PARAMETER (G3 = 1.22541 67024 65177 65D0)
      PARAMETER (PIH = PI/2, RPIH = 2/PI, RPI = 1/PI, RW2 = 1/W2)
      PARAMETER (C1 = PI/(2*W2))
      PARAMETER (GM = 2*(1/G3-4/G1), GP = (4/G1+1/G3)/2)

      DATA GG(-3) /0.27581 56628 30209 31D0/, PP(-3) /-0.75D0/
      DATA GG(-1) /0.81604 89390 98262 98D0/, PP(-1) /-0.25D0/
      DATA GG( 1) /1.10326 26513 20837 26D0/, PP( 1) / 0.25D0/
      DATA GG( 3) /1.08806 52521 31017 31D0/, PP( 3) / 0.75D0/

      DATA BC( 0,1) / 1.00619 92270 14122 57D0/
      DATA BC( 1,1) / 0.00631 99620 31140 72D0/
      DATA BC( 2,1) / 0.00012 56131 27965 64D0/
      DATA BC( 3,1) / 0.00000 52052 40761 57D0/
      DATA BC( 4,1) / 0.00000 03591 84411 39D0/
      DATA BC( 5,1) / 0.00000 00355 85362 89D0/
      DATA BC( 6,1) / 0.00000 00036 05011 66D0/
      DATA BC( 7,1) /-0.00000 00001 26294 10D0/
      DATA BC( 8,1) /-0.00000 00002 96595 12D0/
      DATA BC( 9,1) /-0.00000 00001 18337 70D0/
      DATA BC(10,1) /-0.00000 00000 21655 68D0/
      DATA BC(11,1) / 0.00000 00000 03032 04D0/
      DATA BC(12,1) / 0.00000 00000 03041 10D0/
      DATA BC(13,1) / 0.00000 00000 00530 77D0/
      DATA BC(14,1) /-0.00000 00000 00204 53D0/
      DATA BC(15,1) /-0.00000 00000 00105 49D0/
      DATA BC(16,1) / 0.00000 00000 00005 50D0/
      DATA BC(17,1) / 0.00000 00000 00014 36D0/
      DATA BC(18,1) / 0.00000 00000 00001 14D0/
      DATA BC(19,1) /-0.00000 00000 00001 87D0/
      DATA BC(20,1) /-0.00000 00000 00000 32D0/
      DATA BC(21,1) / 0.00000 00000 00000 26D0/
      DATA BC(22,1) / 0.00000 00000 00000 06D0/
      DATA BC(23,1) /-0.00000 00000 00000 04D0/

      DATA BC( 0,2) / 0.98980 19115 24008 91D0/
      DATA BC( 1,2) /-0.01035 09365 14827 02D0/
      DATA BC( 2,2) /-0.00015 85263 84973 08D0/
      DATA BC( 3,2) /-0.00000 60527 21962 69D0/
      DATA BC( 4,2) /-0.00000 04158 38597 31D0/
      DATA BC( 5,2) /-0.00000 00487 99346 57D0/
      DATA BC( 6,2) /-0.00000 00089 86835 44D0/
      DATA BC( 7,2) /-0.00000 00019 83283 58D0/
      DATA BC( 8,2) /-0.00000 00003 58969 60D0/
      DATA BC( 9,2) /-0.00000 00000 08766 62D0/
      DATA BC(10,2) / 0.00000 00000 25819 45D0/
      DATA BC(11,2) / 0.00000 00000 09780 24D0/
      DATA BC(12,2) / 0.00000 00000 00565 05D0/
      DATA BC(13,2) /-0.00000 00000 00851 66D0/
      DATA BC(14,2) /-0.00000 00000 00270 25D0/
      DATA BC(15,2) / 0.00000 00000 00040 96D0/
      DATA BC(16,2) / 0.00000 00000 00040 50D0/
      DATA BC(17,2) / 0.00000 00000 00001 11D0/
      DATA BC(18,2) /-0.00000 00000 00005 25D0/
      DATA BC(19,2) /-0.00000 00000 00000 70D0/
      DATA BC(20,2) / 0.00000 00000 00000 70D0/
      DATA BC(21,2) / 0.00000 00000 00000 14D0/
      DATA BC(22,2) /-0.00000 00000 00000 10D0/
      DATA BC(23,2) /-0.00000 00000 00000 02D0/

      DATA CC( 0,1) / 0.99128 81656 75147 07D0/
      DATA CC( 1,1) /-0.00850 62567 20022 24D0/
      DATA CC( 2,1) / 0.00019 70491 57408 35D0/
      DATA CC( 3,1) /-0.00000 80377 10166 54D0/
      DATA CC( 4,1) / 0.00000 04554 01498 43D0/
      DATA CC( 5,1) /-0.00000 00323 27352 82D0/
      DATA CC( 6,1) / 0.00000 00027 16130 28D0/
      DATA CC( 7,1) /-0.00000 00002 60644 07D0/
      DATA CC( 8,1) / 0.00000 00000 27882 69D0/
      DATA CC( 9,1) /-0.00000 00000 03267 69D0/
      DATA CC(10,1) / 0.00000 00000 00414 09D0/
      DATA CC(11,1) /-0.00000 00000 00056 17D0/
      DATA CC(12,1) / 0.00000 00000 00008 09D0/
      DATA CC(13,1) /-0.00000 00000 00001 23D0/
      DATA CC(14,1) / 0.00000 00000 00000 20D0/
      DATA CC(15,1) /-0.00000 00000 00000 03D0/

      DATA CC( 0,2) / 1.01476 24350 64637 87D0/
      DATA CC( 1,2) / 0.01449 34617 87809 66D0/
      DATA CC( 2,2) /-0.00025 87162 07241 80D0/
      DATA CC( 3,2) / 0.00000 96912 18911 49D0/
      DATA CC( 4,2) /-0.00000 05261 29313 99D0/
      DATA CC( 5,2) / 0.00000 00363 96854 29D0/
      DATA CC( 6,2) /-0.00000 00030 05472 76D0/
      DATA CC( 7,2) / 0.00000 00002 84827 80D0/
      DATA CC( 8,2) /-0.00000 00000 30182 91D0/
      DATA CC( 9,2) / 0.00000 00000 03511 10D0/
      DATA CC(10,2) /-0.00000 00000 00442 27D0/
      DATA CC(11,2) / 0.00000 00000 00059 70D0/
      DATA CC(12,2) /-0.00000 00000 00008 56D0/
      DATA CC(13,2) / 0.00000 00000 00001 30D0/
      DATA CC(14,2) /-0.00000 00000 00000 21D0/
      DATA CC(15,2) / 0.00000 00000 00000 03D0/

      LEX=.FALSE.
      GO TO 8

#if !defined(CERNLIB_DOUBLE)
      ENTRY EBSIR4(X,NU)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBIR4(X,NU)
#endif
      LEX=.TRUE.

    8 MU=ABS(NU)
      IF(MU .NE. 1 .AND. MU .NE. 2 .AND. MU .NE. 3 .OR.
     1   NU .LT. 0 .AND. X .LE. 0 .OR. NU .GT. 0 .AND. X .LT. 0) THEN
       S=0
       WRITE(ERRTXT,101) X,NU
       IF(.NOT.LEX) CALL MTLPRT(NAMEI ,'C327.1',ERRTXT)
       IF(     LEX) CALL MTLPRT(NAMEIE,'C327.1',ERRTXT)
      ELSEIF(X .EQ. 0) THEN
       S=0
      ELSEIF(NU .EQ. -2) THEN
       IF(LEX) THEN
        S=HF*(1+EXP(-X-X))/SQRT(PIH*X)
       ELSE
        S=COSH(X)/SQRT(PIH*X)
       ENDIF
      ELSEIF(NU .EQ. 2) THEN
       IF(LEX) THEN
        IF(X .LT. HF) THEN
         S=SINH(X)*EXP(-X)/SQRT(PIH*X)
        ELSE
         S=HF*(1-EXP(-X-X))/SQRT(PIH*X)
        ENDIF
       ELSE
        S=SINH(X)/SQRT(PIH*X)
       ENDIF
      ELSEIF(X .LT. 8) THEN
       Y=(HF*X)**2
       XN=PP(NU)
       XL=XN+2
       A0=1
       A1=1+2*Y/((XL+1)*(XN+1))
       A2=1+Y*(4+3*Y/((XL+2)*(XN+2)))/((XL+3)*(XN+1))
       B0=1
       B1=1-Y/(XL+1)
       B2=1-Y*(1-Y/(2*(XL+2)))/(XL+3)
       T1=3+XL
       V1=3-XL
       V3=XL-1
       V2=V3+V3
       C=0
       DO 33 N = 3,30
       C0=C
       T1=T1+2
       T2=T1-1
       T3=T2-1
       T4=T3-1
       T5=T4-1
       T6=T5-1
       V1=V1+1
       V2=V2+1
       V3=V3+1
       U1=N*T4
       E=V3/(U1*T3)
       U2=E*Y
       F1=1+Y*V1/(U1*T1)
       F2=(1+Y*V2/(V3*T2*T5))*U2
       F3=-Y*Y*U2/(T4*T5*T5*T6)
       A=F1*A2+F2*A1+F3*A0
       B=F1*B2+F2*B1+F3*B0
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 34
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
   33  CONTINUE
   34  S=GG(NU)*(HF*X)**PP(NU)*C
       IF(LEX) S=EXP(-X)*S
      ELSE
       K=(MU+1)/2
       R=1/X
       W=SQRT(RPI*R)
       H=16*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = 23,0,-1
       B0=BC(I,K)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       S=RW2*W*(B0-H*B2)
       IF(.NOT.LEX) S=EXP(X)*S
       T=0
       IF(NU .LT. 0) THEN
        H=10*R-1
        ALFA=H+H
        B1=0
        B2=0
        DO 3 I = 15,0,-1
        B0=CC(I,K)+ALFA*B1-B2
        B2=B1
    3   B1=B0
        R=EXP(-X)
        T=W*R*(B0-H*B2)
        IF(LEX) T=R*T
       ENDIF
       S=S+T
      ENDIF
      GO TO 99

#if !defined(CERNLIB_DOUBLE)
      ENTRY  BSKR4(X,NU)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DBSKR4(X,NU)
#endif
      LEX=.FALSE.
      GO TO 9

#if !defined(CERNLIB_DOUBLE)
      ENTRY EBSKR4(X,NU)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBKR4(X,NU)
#endif
      LEX=.TRUE.

    9 MU=ABS(NU)
      IF(MU .NE. 1 .AND. MU .NE. 2 .AND. MU .NE. 3 .OR. X .LE. 0) THEN
       S=0
       WRITE(ERRTXT,101) X,NU
       IF(.NOT.LEX) CALL MTLPRT(NAMEK ,'C327.1',ERRTXT)
       IF(     LEX) CALL MTLPRT(NAMEKE,'C327.1',ERRTXT)
      ELSEIF(MU .EQ. 2) THEN
       S=SQRT(PIH/X)
       IF(.NOT.LEX) S=EXP(-X)*S
      ELSEIF(X .LE. 1) THEN
       A0=PP(-1)
       B=HF*X
       D=-LOG(B)
       F=A0*D
       E=EXP(F)
       G=(GM*A0+GP)*E
       BK=C1*(HF*GM*(E+1/E)+GP*D*SINH(F)/F)
       F=BK
       E=A0**2
       P=HF*C1*G
       Q=HF/G
       C=1
       D=B**2
       BK1=P
       DO 11 N = 1,15
       FN=N
       F=(FN*F+P+Q)/(FN**2-E)
       C=C*D/FN
       P=P/(FN-A0)
       Q=Q/(FN+A0)
       G=C*(P-FN*F)
       H=C*F
       BK=BK+H
       BK1=BK1+G
       IF(H*BK1+ABS(G)*BK .LE. EPS*BK*BK1) GO TO 12
   11  CONTINUE
   12  S=BK
       IF(MU .EQ. 3) S=BK1/B
       IF(LEX) S=EXP(X)*S
      ELSEIF(X .LE. 5) THEN
       XN=4*PP(MU)**2
       A=9-XN
       B=25-XN
       C=768*X**2
       C0=48*X
       A0=1
       A1=(16*X+7+XN)/A
       A2=(C+C0*(XN+23)+XN*(XN+62)+129)/(A*B)
       B0=1
       B1=(16*X+9-XN)/A
       B2=(C+C0*B)/(A*B)+1
       C=0
       DO 24 N = 3,30
       C0=C
       FN=N
       FN2=FN+FN
       FN1=FN2-1
       FN3=FN1/(FN2-3)
       FN4=12*FN**2-(1-XN)
       FN5=16*FN1*X
       RAN=1/((FN2+1)**2-XN)
       F1=FN3*(FN4-20*FN)+FN5
       F2=28*FN-FN4-8+FN5
       F3=FN3*((FN2-5)**2-XN)
       A=(F1*A2+F2*A1+F3*A0)*RAN
       B=(F1*B2+F2*B1+F3*B0)*RAN
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 25
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
   24  CONTINUE
   25  S=C/SQRT(RPIH*X)
       IF(.NOT.LEX) S=EXP(-X)*S
      ELSE
       K=(MU+1)/2
       R=1/X
       Y=5*R
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 13 I = 15,0,-1
       B0=CC(I,K)+ALFA*B1-B2
       B2=B1
   13  B1=B0
       S=SQRT(PIH*R)*(B0-H*B2)
       IF(.NOT.LEX) S=EXP(-X)*S
      ENDIF
#if !defined(CERNLIB_DOUBLE)
   99  BSIR4=S
#endif
#if defined(CERNLIB_DOUBLE)
   99 DBSIR4=S
#endif
      RETURN
  101 FORMAT('INCORRECT ARGUMENT OR INDEX, X = ',1P,E15.6,5X,'NU = ',I5)
      END
