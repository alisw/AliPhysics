*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:07  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION  BSIR3(X,NU)
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION DBSIR3(X,NU)
#include "gen/imp64.inc"
#endif
      CHARACTER*80 ERRTXT
      CHARACTER NAMEI*(*),NAMEK*(*),NAMEIE*(*),NAMEKE*(*)
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAMEI = 'BSIR3', NAMEIE = 'EBSIR3')
      PARAMETER (NAMEK = 'BSKR3', NAMEKE = 'EBSKR3')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAMEI = 'BSIR3/DBSIR3', NAMEIE = 'EBSIR3/DEBIR3')
      PARAMETER (NAMEK = 'BSKR3/DBSKR3', NAMEKE = 'EBSKR3/DEBKR3')
#endif
      LOGICAL LEX

      DIMENSION BC(0:23,2),CC(0:15,2),PP(-2:2),GG(-2:2)

      PARAMETER (EPS = 1D-15)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (W3 = 1.73205 08075 68877 29D0)
      PARAMETER (G1 = 2.67893 85347 07747 63D0)
      PARAMETER (G2 = 1.35411 79394 26400 42D0)
      PARAMETER (PIH = PI/2, RPIH = 2/PI, RPI2 = 1/(2*PI))
      PARAMETER (C1 = 2*PI/(3*W3))
      PARAMETER (GM = 3*(1/G2-3/G1)/2, GP = (3/G1+1/G2)/2)

      DATA PP(-2) /-0.66666 66666 66666 67D0/
      DATA PP(-1) /-0.33333 33333 33333 33D0/
      DATA PP( 1) / 0.33333 33333 33333 33D0/
      DATA PP( 2) / 0.66666 66666 66666 67D0/

      DATA GG(-2) / 0.37328 21739 07395 23D0/
      DATA GG(-1) / 0.73848 81116 21648 31D0/
      DATA GG( 1) / 1.11984 65217 22185 68D0/
      DATA GG( 2) / 1.10773 21674 32472 47D0/

      DATA BC( 0,1) / 1.00458 61710 93207 35D0/
      DATA BC( 1,1) / 0.00467 34791 99873 60D0/
      DATA BC( 2,1) / 0.00009 08034 04815 04D0/
      DATA BC( 3,1) / 0.00000 37262 16110 59D0/
      DATA BC( 4,1) / 0.00000 02520 73237 90D0/
      DATA BC( 5,1) / 0.00000 00227 82110 77D0/
      DATA BC( 6,1) / 0.00000 00012 91332 28D0/
      DATA BC( 7,1) /-0.00000 00006 11915 16D0/
      DATA BC( 8,1) /-0.00000 00003 75616 85D0/
      DATA BC( 9,1) /-0.00000 00001 16415 46D0/
      DATA BC(10,1) /-0.00000 00000 14443 25D0/
      DATA BC(11,1) / 0.00000 00000 05373 69D0/
      DATA BC(12,1) / 0.00000 00000 03074 27D0/
      DATA BC(13,1) / 0.00000 00000 00297 66D0/
      DATA BC(14,1) /-0.00000 00000 00265 20D0/
      DATA BC(15,1) /-0.00000 00000 00091 37D0/
      DATA BC(16,1) / 0.00000 00000 00015 52D0/
      DATA BC(17,1) / 0.00000 00000 00014 12D0/
      DATA BC(18,1) /-0.00000 00000 00000 23D0/
      DATA BC(19,1) /-0.00000 00000 00001 98D0/
      DATA BC(20,1) /-0.00000 00000 00000 13D0/
      DATA BC(21,1) / 0.00000 00000 00000 29D0/
      DATA BC(22,1) / 0.00000 00000 00000 03D0/
      DATA BC(23,1) /-0.00000 00000 00000 05D0/

      DATA BC( 0,2) / 0.99363 49867 16925 14D0/
      DATA BC( 1,2) /-0.00646 71526 00616 03D0/
      DATA BC( 2,2) /-0.00010 60188 22351 55D0/
      DATA BC( 3,2) /-0.00000 41406 57716 24D0/
      DATA BC( 4,2) /-0.00000 02916 95418 21D0/
      DATA BC( 5,2) /-0.00000 00365 71574 33D0/
      DATA BC( 6,2) /-0.00000 00075 81590 37D0/
      DATA BC( 7,2) /-0.00000 00019 23008 52D0/
      DATA BC( 8,2) /-0.00000 00004 20438 80D0/
      DATA BC( 9,2) /-0.00000 00000 39372 04D0/
      DATA BC(10,2) / 0.00000 00000 19007 44D0/
      DATA BC(11,2) / 0.00000 00000 10137 64D0/
      DATA BC(12,2) / 0.00000 00000 01331 30D0/
      DATA BC(13,2) /-0.00000 00000 00676 92D0/
      DATA BC(14,2) /-0.00000 00000 00311 72D0/
      DATA BC(15,2) / 0.00000 00000 00011 87D0/
      DATA BC(16,2) / 0.00000 00000 00040 21D0/
      DATA BC(17,2) / 0.00000 00000 00004 78D0/
      DATA BC(18,2) /-0.00000 00000 00004 74D0/
      DATA BC(19,2) /-0.00000 00000 00001 16D0/
      DATA BC(20,2) / 0.00000 00000 00000 59D0/
      DATA BC(21,2) / 0.00000 00000 00000 21D0/
      DATA BC(22,2) /-0.00000 00000 00000 08D0/
      DATA BC(23,2) /-0.00000 00000 00000 03D0/

      DATA CC( 0,1) / 0.99353 64122 76093 39D0/
      DATA CC( 1,1) /-0.00631 44392 60798 63D0/
      DATA CC( 2,1) / 0.00014 30095 80961 13D0/
      DATA CC( 3,1) /-0.00000 57870 60592 03D0/
      DATA CC( 4,1) / 0.00000 03265 50333 20D0/
      DATA CC( 5,1) /-0.00000 00231 23231 95D0/
      DATA CC( 6,1) / 0.00000 00019 39555 14D0/
      DATA CC( 7,1) /-0.00000 00001 85897 89D0/
      DATA CC( 8,1) / 0.00000 00000 19868 42D0/
      DATA CC( 9,1) /-0.00000 00000 02326 79D0/
      DATA CC(10,1) / 0.00000 00000 00294 68D0/
      DATA CC(11,1) /-0.00000 00000 00039 95D0/
      DATA CC(12,1) / 0.00000 00000 00005 75D0/
      DATA CC(13,1) /-0.00000 00000 00000 87D0/
      DATA CC(14,1) / 0.00000 00000 00000 14D0/
      DATA CC(15,1) /-0.00000 00000 00000 02D0/

      DATA CC( 0,2) / 1.00914 95380 72789 40D0/
      DATA CC( 1,2) / 0.00897 12068 42483 60D0/
      DATA CC( 2,2) /-0.00017 13895 98261 54D0/
      DATA CC( 3,2) / 0.00000 65547 92549 82D0/
      DATA CC( 4,2) /-0.00000 03595 19190 49D0/
      DATA CC( 5,2) / 0.00000 00250 24412 19D0/
      DATA CC( 6,2) /-0.00000 00020 74924 13D0/
      DATA CC( 7,2) / 0.00000 00001 97223 67D0/
      DATA CC( 8,2) /-0.00000 00000 20946 47D0/
      DATA CC( 9,2) / 0.00000 00000 02440 93D0/
      DATA CC(10,2) /-0.00000 00000 00307 91D0/
      DATA CC(11,2) / 0.00000 00000 00041 61D0/
      DATA CC(12,2) /-0.00000 00000 00005 97D0/
      DATA CC(13,2) / 0.00000 00000 00000 91D0/
      DATA CC(14,2) /-0.00000 00000 00000 14D0/
      DATA CC(15,2) / 0.00000 00000 00000 02D0/

      LEX=.FALSE.
      GO TO 8

#if !defined(CERNLIB_DOUBLE)
      ENTRY  EBSIR3(X,NU)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBIR3(X,NU)
#endif
      LEX=.TRUE.

    8 MU=ABS(NU)
      IF(MU .NE. 1 .AND. MU .NE. 2 .OR. NU .LT. 0 .AND. X .LE. 0
     1   .OR. NU .GT. 0 .AND. X .LT. 0) THEN
       S=0
       WRITE(ERRTXT,101) X,NU
       IF(.NOT.LEX) CALL MTLPRT(NAMEI ,'C340.1',ERRTXT)
       IF(     LEX) CALL MTLPRT(NAMEIE,'C340.1',ERRTXT)
      ELSEIF(X .EQ. 0) THEN
       S=0
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
       R=1/X
       W=SQRT(RPI2*R)
       H=16*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 2 I = 23,0,-1
       B0=BC(I,MU)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       S=W*(B0-H*B2)
       IF(.NOT.LEX) S=EXP(X)*S
       T=0
       IF(NU .LT. 0) THEN
        H=10*R-1
        ALFA=H+H
        B1=0
        B2=0
        DO 3 I = 15,0,-1
        B0=CC(I,MU)+ALFA*B1-B2
        B2=B1
    3   B1=B0
        R=EXP(-X)
        T=W3*W*R*(B0-H*B2)
        IF(LEX) T=R*T
       END IF
       S=S+T
      END IF
      GO TO 99

#if !defined(CERNLIB_DOUBLE)
      ENTRY  BSKR3(X,NU)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DBSKR3(X,NU)
#endif
      LEX=.FALSE.
      GO TO 9

#if !defined(CERNLIB_DOUBLE)
      ENTRY  EBSKR3(X,NU)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBKR3(X,NU)
#endif
      LEX=.TRUE.

    9 MU=ABS(NU)
      IF(MU .NE. 1 .AND. MU .NE. 2 .OR. X .LE. 0) THEN
       S=0
       WRITE(ERRTXT,101) X,NU
       IF(.NOT.LEX) CALL MTLPRT(NAMEK ,'C340.1',ERRTXT)
       IF(     LEX) CALL MTLPRT(NAMEKE,'C340.1',ERRTXT)
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
       IF(MU .EQ. 2) S=BK1/B
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
       R=1/X
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 13 I = 15,0,-1
       B0=CC(I,MU)+ALFA*B1-B2
       B2=B1
   13  B1=B0
       S=SQRT(PIH*R)*(B0-H*B2)
       IF(.NOT.LEX) S=EXP(-X)*S
      END IF
#if !defined(CERNLIB_DOUBLE)
   99  BSIR3=S
#endif
#if defined(CERNLIB_DOUBLE)
   99 DBSIR3=S
#endif
      RETURN
  101 FORMAT('INCORRECT ARGUMENT OR INDEX, X = ',1P,E15.6,' NU = ',I5)
      END
