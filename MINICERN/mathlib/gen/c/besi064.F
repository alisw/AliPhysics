*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:00  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION  BESI0(X)
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION DBESI0(X)
#include "gen/imp64.inc"
#endif
      LOGICAL LEX
      CHARACTER NAME0*(*),NAME1*(*),NAME0E*(*),NAME1E*(*)
      CHARACTER*80 ERRTXT
      DIMENSION CI(0:24,0:1),CK(0:16,0:1)

#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME0 = 'BESK0', NAME0E = 'EBESK0')
      PARAMETER (NAME1 = 'BESK1', NAME1E = 'EBESK1')
#endif
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME0 = 'BESK0/DBESK0', NAME0E = 'EBESK0/DBESK0')
      PARAMETER (NAME1 = 'BESK1/DBESK1', NAME1E = 'EBESK1/DBESK1')
#endif
      PARAMETER (EPS=1D-15)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (CE = 0.57721 56649 01532 86D0)
      PARAMETER (PIH = PI/2, RPIH = 2/PI, RPI2 = 1/(2*PI))

      DATA CI( 0,0) /+1.00827 92054 58740 032D0/
      DATA CI( 1,0) /+0.00844 51226 24920 943D0/
      DATA CI( 2,0) /+0.00017 27006 30777 567D0/
      DATA CI( 3,0) /+0.00000 72475 91099 959D0/
      DATA CI( 4,0) /+0.00000 05135 87726 878D0/
      DATA CI( 5,0) /+0.00000 00568 16965 808D0/
      DATA CI( 6,0) /+0.00000 00085 13091 223D0/
      DATA CI( 7,0) /+0.00000 00012 38425 364D0/
      DATA CI( 8,0) /+0.00000 00000 29801 672D0/
      DATA CI( 9,0) /-0.00000 00000 78956 698D0/
      DATA CI(10,0) /-0.00000 00000 33127 128D0/
      DATA CI(11,0) /-0.00000 00000 04497 339D0/
      DATA CI(12,0) /+0.00000 00000 01799 790D0/
      DATA CI(13,0) /+0.00000 00000 00965 748D0/
      DATA CI(14,0) /+0.00000 00000 00038 604D0/
      DATA CI(15,0) /-0.00000 00000 00104 039D0/
      DATA CI(16,0) /-0.00000 00000 00023 950D0/
      DATA CI(17,0) /+0.00000 00000 00009 554D0/
      DATA CI(18,0) /+0.00000 00000 00004 443D0/
      DATA CI(19,0) /-0.00000 00000 00000 859D0/
      DATA CI(20,0) /-0.00000 00000 00000 709D0/
      DATA CI(21,0) /+0.00000 00000 00000 087D0/
      DATA CI(22,0) /+0.00000 00000 00000 112D0/
      DATA CI(23,0) /-0.00000 00000 00000 012D0/
      DATA CI(24,0) /-0.00000 00000 00000 018D0/

      DATA CI( 0,1) /+0.97580 06023 26285 926D0/
      DATA CI( 1,1) /-0.02446 74429 63276 385D0/
      DATA CI( 2,1) /-0.00027 72053 60763 829D0/
      DATA CI( 3,1) /-0.00000 97321 46728 020D0/
      DATA CI( 4,1) /-0.00000 06297 24238 640D0/
      DATA CI( 5,1) /-0.00000 00659 61142 154D0/
      DATA CI( 6,1) /-0.00000 00096 13872 919D0/
      DATA CI( 7,1) /-0.00000 00014 01140 901D0/
      DATA CI( 8,1) /-0.00000 00000 47563 167D0/
      DATA CI( 9,1) /+0.00000 00000 81530 681D0/
      DATA CI(10,1) /+0.00000 00000 35408 148D0/
      DATA CI(11,1) /+0.00000 00000 05102 564D0/
      DATA CI(12,1) /-0.00000 00000 01804 409D0/
      DATA CI(13,1) /-0.00000 00000 01023 594D0/
      DATA CI(14,1) /-0.00000 00000 00052 678D0/
      DATA CI(15,1) /+0.00000 00000 00107 094D0/
      DATA CI(16,1) /+0.00000 00000 00026 120D0/
      DATA CI(17,1) /-0.00000 00000 00009 561D0/
      DATA CI(18,1) /-0.00000 00000 00004 713D0/
      DATA CI(19,1) /+0.00000 00000 00000 829D0/
      DATA CI(20,1) /+0.00000 00000 00000 743D0/
      DATA CI(21,1) /-0.00000 00000 00000 080D0/
      DATA CI(22,1) /-0.00000 00000 00000 117D0/
      DATA CI(23,1) /+0.00000 00000 00000 011D0/
      DATA CI(24,1) /+0.00000 00000 00000 019D0/

      DATA CK( 0,0) /+0.98840 81742 30825 800D0/
      DATA CK( 1,0) /-0.01131 05046 46928 281D0/
      DATA CK( 2,0) /+0.00026 95326 12762 724D0/
      DATA CK( 3,0) /-0.00001 11066 85196 665D0/
      DATA CK( 4,0) /+0.00000 06325 75108 500D0/
      DATA CK( 5,0) /-0.00000 00450 47337 641D0/
      DATA CK( 6,0) /+0.00000 00037 92996 456D0/
      DATA CK( 7,0) /-0.00000 00003 64547 179D0/
      DATA CK( 8,0) /+0.00000 00000 39043 756D0/
      DATA CK( 9,0) /-0.00000 00000 04579 936D0/
      DATA CK(10,0) /+0.00000 00000 00580 811D0/
      DATA CK(11,0) /-0.00000 00000 00078 832D0/
      DATA CK(12,0) /+0.00000 00000 00011 360D0/
      DATA CK(13,0) /-0.00000 00000 00001 727D0/
      DATA CK(14,0) /+0.00000 00000 00000 275D0/
      DATA CK(15,0) /-0.00000 00000 00000 046D0/
      DATA CK(16,0) /+0.00000 00000 00000 008D0/

      DATA CK( 0,1) /+1.03595 08587 72358 331D0/
      DATA CK( 1,1) /+0.03546 52912 43331 114D0/
      DATA CK( 2,1) /-0.00046 84750 28166 889D0/
      DATA CK( 3,1) /+0.00001 61850 63810 053D0/
      DATA CK( 4,1) /-0.00000 08451 72048 124D0/
      DATA CK( 5,1) /+0.00000 00571 32218 103D0/
      DATA CK( 6,1) /-0.00000 00046 45554 607D0/
      DATA CK( 7,1) /+0.00000 00004 35417 339D0/
      DATA CK( 8,1) /-0.00000 00000 45757 297D0/
      DATA CK( 9,1) /+0.00000 00000 05288 133D0/
      DATA CK(10,1) /-0.00000 00000 00662 613D0/
      DATA CK(11,1) /+0.00000 00000 00089 048D0/
      DATA CK(12,1) /-0.00000 00000 00012 726D0/
      DATA CK(13,1) /+0.00000 00000 00001 921D0/
      DATA CK(14,1) /-0.00000 00000 00000 305D0/
      DATA CK(15,1) /+0.00000 00000 00000 050D0/
      DATA CK(16,1) /-0.00000 00000 00000 009D0/

      NU=0
      LEX=.FALSE.
      GO TO 6

#if !defined(CERNLIB_DOUBLE)
      ENTRY  EBESI0(X)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBSI0(X)
#endif
      NU=0
      LEX=.TRUE.
      GO TO 6

#if !defined(CERNLIB_DOUBLE)
      ENTRY  BESI1(X)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DBESI1(X)
#endif
      NU=1
      LEX=.FALSE.
      GO TO 6

#if !defined(CERNLIB_DOUBLE)
      ENTRY  EBESI1(X)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBSI1(X)
#endif
      NU=1
      LEX=.TRUE.

    6 V=ABS(X)
      IF(V .LT. 8) THEN
       Y=(HF*V)**2
       XL=NU+2
       A0=1
       A1=1+2*Y/((XL+1)*(XL-1))
       A2=1+Y*(4+3*Y/((XL+2)*XL))/((XL+3)*(XL-1))
       B0=1
       B1=1-Y/(XL+1)
       B2=1-Y*(1-Y/(2*(XL+2)))/(XL+3)
       W1=3+XL
       V1=3-XL
       V3=XL-1
       V2=V3+V3
       C=0
       DO 3 N = 3,30
       C0=C
       FN=N
       W1=W1+2
       W2=W1-1
       W3=W2-1
       W4=W3-1
       W5=W4-1
       W6=W5-1
       V1=V1+1
       V2=V2+1
       V3=V3+1
       U1=FN*W4
       E=V3/(U1*W3)
       U2=E*Y
       F1=1+Y*V1/(U1*W1)
       F2=(1+Y*V2/(V3*W2*W5))*U2
       F3=-Y*Y*U2/(W4*W5*W5*W6)
       A=F1*A2+F2*A1+F3*A0
       B=F1*B2+F2*B1+F3*B0
       C=A/B
       IF(ABS(C0-C) .LT. EPS*ABS(C)) GO TO 4
       A0=A1
       A1=A2
       A2=A
       B0=B1
       B1=B2
       B2=B
    3  CONTINUE
    4  H=C
       IF(NU .EQ. 1) H=HF*X*H
       IF(LEX) H=EXP(-V)*H
      ELSE
       R=1/V
       H=16*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 24,0,-1
       B0=CI(I,NU)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=SQRT(RPI2*R)*(B0-H*B2)
       IF(NU*X .LT. 0) H=-H
       IF(.NOT.LEX) H=EXP(V)*H
      ENDIF
      GO TO 9

#if !defined(CERNLIB_DOUBLE)
      ENTRY  BESK0(X)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DBESK0(X)
#endif
      NU=0
      LEX=.FALSE.
      GO TO 8

#if !defined(CERNLIB_DOUBLE)
      ENTRY  EBESK0(X)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBSK0(X)
#endif
      NU=0
      LEX=.TRUE.
      GO TO 8

#if !defined(CERNLIB_DOUBLE)
      ENTRY  BESK1(X)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DBESK1(X)
#endif
      NU=1
      LEX=.FALSE.
      GO TO 8

#if !defined(CERNLIB_DOUBLE)
      ENTRY  EBESK1(X)
#endif
#if defined(CERNLIB_DOUBLE)
      ENTRY DEBSK1(X)
#endif
      NU=1
      LEX=.TRUE.

    8 IF(X .LE. 0) THEN
       H=0
       WRITE(ERRTXT,101) X
       IF(NU .EQ. 0 .AND. .NOT.LEX) CALL MTLPRT(NAME0 ,'C313.1',ERRTXT)
       IF(NU .EQ. 0 .AND.      LEX) CALL MTLPRT(NAME0E,'C313.1',ERRTXT)
       IF(NU .EQ. 1 .AND. .NOT.LEX) CALL MTLPRT(NAME1 ,'C313.1',ERRTXT)
       IF(NU .EQ. 1 .AND.      LEX) CALL MTLPRT(NAME1E,'C313.1',ERRTXT)
      ELSEIF(X .LT. 1) THEN
       B=HF*X
       BK=-(LOG(B)+CE)
       F=BK
       P=HF
       Q=HF
       C=1
       D=B**2
       BK1=P
       DO 11 N = 1,15
       FN=N
       RFN=1/FN
       P=P*RFN
       Q=Q*RFN
       F=(F+P+Q)*RFN
       C=C*D*RFN
       G=C*(P-FN*F)
       H=C*F
       BK=BK+H
       BK1=BK1+G
       IF(BK1*H+ABS(G)*BK .LE. EPS*BK*BK1) GO TO 12
   11  CONTINUE
   12  H=BK
       IF(NU .EQ. 1) H=BK1/B
       IF(LEX) H=EXP(X)*H
      ELSEIF(X .LE. 5) THEN
       XN=4*NU**2
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
   25  H=C/SQRT(RPIH*X)
       IF(.NOT.LEX) H=EXP(-X)*H
      ELSE
       R=1/X
       H=10*R-1
       ALFA=H+H
       B1=0
       B2=0
       DO 23 I = 16,0,-1
       B0=CK(I,NU)+ALFA*B1-B2
       B2=B1
   23  B1=B0
       H=SQRT(PIH*R)*(B0-H*B2)
       IF(.NOT.LEX) H=EXP(-X)*H
      ENDIF
    9 CONTINUE
#if !defined(CERNLIB_DOUBLE)
       BESI0=H
#endif
#if defined(CERNLIB_DOUBLE)
      DBESI0=H
#endif
      RETURN
  101 FORMAT(' NON-POSITIVE ARGUMENT X = ',1P,E15.6)
      END
