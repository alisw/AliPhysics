*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE DEFS1(I,J,K)
C
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (15-JAN-1984)
C
#include "geant321/mxgkgh.inc"
#include "geant321/s_blank.inc"
C
      DATA PI/3.141592653589793/
      PX=PV(1,I)
      PY=PV(2,I)
      PZ=PV(3,I)
      PT2 = PV(1,J)**2+PV(2,J)**2
      IF(PT2.GT.0.) THEN
      CALL LENGTX(J,P)
      COST=PV(3,J)/P
      SINT1=SQRT(ABS((1.-COST)*(1.+COST)))
      SINT2=SQRT(PT2)/P
      SINT=0.5*(SINT1+SINT2)
      PH=PI*0.5
      IF(PV(2,J).LT.0.) PH=PI*1.5
      IF(ABS(PV(1,J)).GT.1.E-6) PH=ATAN2(PV(2,J),PV(1,J))
      COSP=COS(PH)
      SINP=SIN(PH)
      PV(1,K)= COST*COSP*PX-     SINP*PY+SINT*COSP*PZ
      PV(2,K)= COST*SINP*PX+     COSP*PY+SINT*SINP*PZ
      PV(3,K)=-SINT     *PX             +COST     *PZ
      ELSE
      PV(1,K)=PX
      PV(2,K)=PY
      PV(3,K)=PZ
C --- Take the case of theta=pi into account (MR/NVE 27-sep-1990) ---
      IF (PV(3,J) .LT. 0.) PV(3,K)=-PZ
      ENDIF
      END
