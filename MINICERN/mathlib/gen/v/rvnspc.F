*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION RVNSPC(R,RHO,D)

      DIMENSION AZ(4),BZ(4),PZ(4),QZ(4)

      PARAMETER (Z1 = 1)
      PARAMETER (C1 = 4*Z1/3, C2 = 2*Z1/3, C3 = 4*Z1/9, C4 = Z1/3)
      PARAMETER (PI = 3.14159 265D0, CL = 1.38629 436D0)
      PARAMETER (SF = 4*PI/3, SFH = 2*PI/3, C0 = 2*PI/3-8*Z1/9)
      PARAMETER (PIH = PI/2, Z10 = 10, HF = Z1/2, CA = Z10**(-4))

      DATA AZ
     1/1.45133 8556D-2, 3.74253 9571D-2, 3.58998 0090D-2,
     2 9.66633 8350D-2/
      DATA BZ
     1/4.41839 8230D-3, 3.32852 1016D-2, 6.88029 5505D-2,
     2 1.24985 9468D-1/
      DATA PZ
     1/1.73631 4854D-2, 4.75740 4429D-2, 6.26076 1942D-2,
     2 4.43251 5145D-1/
      DATA QZ
     1/5.26378 9328D-3, 4.06946 8414D-2, 9.20010 9374D-2,
     2 2.49983 6641D-1/

      RS=ABS(R)
      RC=ABS(RHO)
      DA=ABS(D)
      DR=RS-RC
      RS2=RS**2
      RS3=RS*RS2
      IF(RC .EQ. 0 .OR. RS .EQ. 0 .OR. DA .GE. RS+RC) THEN
       V=0
      ELSEIF(DR .LE. DA .AND. DA .LE. -DR) THEN
       V=SF*RS3
      ELSEIF(DA .EQ. RC .AND. RS .EQ. 2*DA) THEN
       V=C0*RS3
      ELSEIF(DA .EQ. 0) THEN
       V=RS3
       IF(RS .GT. RC) V=V-SQRT(RS2-RC**2)**3
       V=SF*V
      ELSE
       BP=DA+RC
       BM=DA-RC
       BP2=BP**2
       A=MAX(RS2,BP2)
       B=MIN(RS2,BP2)
       C=BM**2
       AB=A-B
       AC=A-C
       BC=B-C
       S=BP*BM
       IF(DA .NE. DR) THEN
        Y=AB/AC
        YL=LOG(Y)
        PA=AZ(1)
        PB=BZ(1)
        DO 1 I = 2,4
        PA=PA*Y+AZ(I)
    1   PB=PB*Y+BZ(I)
        HK=CL+PA*Y-YL*(HF+PB*Y)
        PA=PZ(1)
        PB=QZ(1)
        DO 2 I = 2,4
        PA=PA*Y+PZ(I)
    2   PB=PB*Y+QZ(I)
        HE=1+(PA-YL*PB)*Y
        IF(C .NE. 0) THEN
         YKP=SQRT(AB/AC)
         EE=YKP
         AM0=1
         PP=B/C
         IF(PP .GT. 0) THEN
          CC=1
          PP=SQRT(PP)
          DD=1/PP
         ELSE
          GG=1-PP
          FF=YKP**2-PP
          PP=SQRT(FF/GG)
          DD=-BC/(AC*GG*PP)
          CC=0
         ENDIF
    3    FF=CC
         CC=DD/PP+CC
         GG=EE/PP
         DD=2*(FF*GG+DD)
         PP=GG+PP
         GG=AM0
         AM0=YKP+AM0
         IF(ABS(GG-YKP) .GT. CA*GG) THEN
          YKP=2*SQRT(EE)
          EE=YKP*AM0
          GO TO 3
         ENDIF
         H3=PIH*(CC*AM0+DD)/(AM0*(AM0+PP))
        ENDIF
       ENDIF
       IF(DA .LT. DR) THEN
        IF(C .EQ. 0) THEN
         V=SFH*RS3+C3*SQRT(A)*(AB*HK-2*(A+AB)*HE)
        ELSE
         V=C1*(H3*A**2*S/C-HK*(A*S-C4*AB*AC)-HE*AC*(S+C2*(AB+AC)))/
     1     SQRT(AC)
         IF(RC .GT. DA) V=V+SF*RS3
        ENDIF
       ELSEIF(DA .EQ. DR) THEN
        V=C1*(RS3*ATAN2(2*SQRT(DA*RC),BM)-SQRT(AC)*(S+C2*AC))
       ELSE
        IF(C .EQ. 0) THEN
         V=SFH*RS3+C3*(AB*(B-2*AB)*HK+2*A*(AB-B)*HE)/SQRT(A)
        ELSE
         V=C1*(H3*B**2*S/C+HK*(S*(AB-B)+C4*AB*(BC-2*AB))
     1    -HE*AC*(S-C2*(AB-BC)))/SQRT(AC)
         IF(RC .GT. DA) V=V+SF*RS3
        ENDIF
       ENDIF
      ENDIF
      RVNSPC=V
      RETURN
      END
#endif
