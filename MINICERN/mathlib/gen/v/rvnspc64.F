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
      FUNCTION DVNSPC(R,RHO,D)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='DVNSPC')
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION RVNSPC(R,RHO,D)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RVNSPC')
#endif
C
C     Based on
C     F. LAMARCHE and C. LEROY, Evaluation of the volume of
C     a sphere with a cylinder by elliptic integrals,
C     Computer Phys. Comm. 59 (1990) 359-369
C
      PARAMETER (Z1 = 1)
      PARAMETER (C1 = 4*Z1/3, C2 = 2*Z1/3, C3 = 4*Z1/9, C4 = Z1/3)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (SF = 4*PI/3, SFH = 2*PI/3, C0 = 2*PI/3-8*Z1/9)

      RS=ABS(R)
      RC=ABS(RHO)
      DA=ABS(D)
      DR=RS-RC
      RS2=RS**2
      RS3=RS2*RS
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
       XK2=BC/AC
       XK=SQRT(XK2)
       IF(DA .LT. DR) THEN
        IF(C .EQ. 0) THEN
#if defined(CERNLIB_DOUBLE)
         V=SFH*RS3+C3*SQRT(A)*(AB*DELIKC(XK)-2*(A+AB)*DELIEC(XK))
        ELSE
         V=C1*(DELI3C(SQRT(1-XK2),XK2,B/C)*A**2*S/C
     1    -DELIKC(XK)*(A*S-C4*AB*AC)
     2    -DELIEC(XK)*AC*(S+C2*(AB+AC)))/SQRT(AC)
#endif
#if !defined(CERNLIB_DOUBLE)
         V=SFH*RS3+C3*SQRT(A)*(AB*RELIKC(XK)-2*(A+AB)*RELIEC(XK))
        ELSE
         V=C1*(RELI3C(SQRT(1-XK2),XK2,B/C)*A**2*S/C
     1    -RELIKC(XK)*(A*S-C4*AB*AC)
     2    -RELIEC(XK)*AC*(S+C2*(AB+AC)))/SQRT(AC)
#endif
         IF(RC .GT. DA) V=V+SF*RS3
        ENDIF
       ELSEIF(DA .EQ. DR) THEN
        V=C1*(RS3*ATAN2(2*SQRT(DA*RC),BM)-SQRT(AC)*(S+C2*AC))
       ELSE
        IF(C .EQ. 0) THEN
#if defined(CERNLIB_DOUBLE)
         V=SFH*RS3+C3*(AB*(B-2*AB)*DELIKC(XK)+2*A*(AB-B)*DELIEC(XK))/
     1     SQRT(A)
        ELSE
         V=C1*(DELI3C(SQRT(1-XK2),XK2,B/C)*B**2*S/C
     1    +DELIKC(XK)*(S*(AB-B)+C4*AB*(BC-2*AB))
     2    -DELIEC(XK)*AC*(S-C2*(AB-BC)))/SQRT(AC)
#endif
#if !defined(CERNLIB_DOUBLE)
         V=SFH*RS3+C3*(AB*(B-2*AB)*RELIKC(XK)+2*A*(AB-B)*RELIEC(XK))/
     1     SQRT(A)
        ELSE
         V=C1*(RELI3C(SQRT(1-XK2),XK2,B/C)*B**2*S/C
     1    +RELIKC(XK)*(S*(AB-B)+C4*AB*(BC-2*AB))
     2    -RELIEC(XK)*AC*(S-C2*(AB-BC)))/SQRT(AC)
#endif
         IF(RC .GT. DA) V=V+SF*RS3
        ENDIF
       ENDIF
      ENDIF
#if defined(CERNLIB_DOUBLE)
      DVNSPC=V
#endif
#if !defined(CERNLIB_DOUBLE)
      RVNSPC=V
#endif
      RETURN
      END
