#include "isajet/pilot.h"
      SUBROUTINE SSHCC
C-----------------------------------------------------------------------
C     Calculates the decay widths of all neutral Higgses into all
C     possible pairs of charginos.
C
C     Bisset's CHGINO
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sspar.inc"
#include "isajet/sssm.inc"
#include "isajet/sstype.inc"
C
      DOUBLE PRECISION PI,SR2,XM,THETX,YM,THETY,SGL,CGL,SGR,CGR
     $,MW1,MW2,THETM,THETP,G2,GP2,BETA,ALPHA,T1,MH,M1,M2,LAMB
     $,DWID,TEMP,TEMPXY
      DOUBLE PRECISION MHI(3),IDHI(3),SHP(3),SHM(3),SH(3),PH(3)
      DOUBLE PRECISION SSDLAM
      REAL WID
      INTEGER NUMH,IDHHA
C
      PI=4.*ATAN(1.D0)
      SR2=SQRT(2.D0)
      XM=1./TAN(GAMMAL)
      THETX=SIGN(1.D0,XM)
      YM=1./TAN(GAMMAR)
      THETY=SIGN(1.D0,YM)
      SGL=1/(DSQRT(1+XM**2))
      CGL=SGL*XM
      SGR=1/(DSQRT(1+YM**2))
      CGR=SGR*YM
      MW1=DBLE(ABS(AMW1SS))
      MW2=DBLE(ABS(AMW2SS))
      THETM=SIGN(1.,AMW1SS)
      THETP=SIGN(1.,AMW2SS)
      G2=4*PI*ALFAEM/SN2THW
      GP2=4*PI*ALFAEM/(1-SN2THW)
      BETA=ATAN(1.0/RV2V1)
      ALPHA=ALFAH
C
C          Arrays for loops
C
      MHI(1)=AMHL
      MHI(2)=AMHH
      MHI(3)=AMHA
      IDHI(1)=ISHL
      IDHI(2)=ISHH
      IDHI(3)=ISHA
C          The following came from Bisset's MASZIN, but with L,H,P
C          replaced by a generic H and a subscript.
      TEMPXY=0.5*THETX*THETY*(-THETP)
      SHP(1)=SIN(ALPHA)*CGR*SGL+COS(ALPHA)*CGL*SGR
      SHP(1)=SHP(1)*TEMPXY
      SHM(1)=SIN(ALPHA)*SGR*CGL+COS(ALPHA)*SGL*CGR
      SHM(1)=SHM(1)*0.5*THETM
      SH(1)=-THETX*SGR*SGL*SIN(ALPHA)*THETM
      PH(1)=-SH(1)
      T1=THETX*CGL*CGR*COS(ALPHA)*THETM
      SH(1)=SH(1)+T1
      PH(1)=PH(1)-T1
      T1=THETY*SGL*SGR*COS(ALPHA)*THETP
      SH(1)=SH(1)-T1
      PH(1)=PH(1)-T1
      T1=THETY*CGL*CGR*SIN(ALPHA)*THETP
      SH(1)=SH(1)+T1
      PH(1)=PH(1)+T1
      SH(1)=0.5*SH(1)
      PH(1)=0.5*PH(1)
      SHP(2)=COS(ALPHA)*CGR*SGL-SIN(ALPHA)*CGL*SGR
      SHP(2)=SHP(2)*TEMPXY
      SHM(2)=COS(ALPHA)*SGR*CGL-SIN(ALPHA)*SGL*CGR
      SHM(2)=SHM(2)*0.5*THETM
      SH(2)=-THETX*SGR*SGL*COS(ALPHA)*THETM
      PH(2)=-SH(2)
      T1=THETX*CGL*CGR*SIN(ALPHA)*THETM
      SH(2)=SH(2)-T1
      PH(2)=PH(2)+T1
      T1=THETY*SGL*SGR*SIN(ALPHA)*THETP
      SH(2)=SH(2)+T1
      PH(2)=PH(2)+T1
      T1=THETY*CGL*CGR*COS(ALPHA)*THETP
      SH(2)=SH(2)+T1
      PH(2)=PH(2)+T1
      SH(2)=0.5*SH(2)
      PH(2)=0.5*PH(2)
      SHP(3)=SIN(BETA)*CGR*SGL+COS(BETA)*CGL*SGR
      SHP(3)=SHP(3)*0.5*THETX*THETY*(-THETP)
      SHM(3)=SIN(BETA)*SGR*CGL+COS(BETA)*SGL*CGR
      SHM(3)=SHM(3)*0.5*THETM
      SH(3)=-THETX*SGR*SGL*SIN(BETA)*THETM
      PH(3)=-SH(3)
      T1=THETX*CGL*CGR*COS(BETA)*THETM
      SH(3)=SH(3)+T1
      PH(3)=PH(3)-T1
      T1=THETY*SGL*SGR*COS(BETA)*THETP
      SH(3)=SH(3)+T1
      PH(3)=PH(3)+T1
      T1=THETY*CGL*CGR*SIN(BETA)*THETP
      SH(3)=SH(3)-T1
      PH(3)=PH(3)-T1
      SH(3)=0.5*SH(3)
      PH(3)=0.5*PH(3)
C
C          Loop over neutral Higgs
C
      DO 100 NUMH=1,3
        MH=MHI(NUMH)
        IDHHA=IDHI(NUMH)
C          w1 + w1
        M1=ABS(AMW1SS)
        M2=M1
        IF(MH.GT.M1+M2) THEN
          LAMB=SSDLAM(MH**2,M1**2,M2**2)
          TEMP=1-4*M1**2/MH**2
          DWID=G2*MH*SHM(NUMH)**2/(4.0*PI)
          DWID=DWID*SQRT(TEMP**3)
          WID=DWID
          CALL SSSAVE(IDHHA,WID,ISW1,-ISW1,0,0,0)
        ENDIF
C          w2 + w2
        M1=ABS(AMW2SS)
        M2=M1
        IF(MH.GT.M1+M2) THEN
          TEMP=1-4*M1**2/MH**2
          DWID=G2*MH*SHP(NUMH)**2/(4*PI)
          DWID=DWID*SQRT(TEMP**3)
          WID=DWID
          CALL SSSAVE(IDHHA,WID,ISW2,-ISW2,0,0,0)
        ENDIF
C          w1 + w2
        M1=ABS(AMW1SS)
        M2=ABS(AMW2SS)
        IF(MH.GT.M1+M2) THEN
          LAMB=SSDLAM(MH**2,M1**2,M2**2)
          DWID=PH(NUMH)**2*(MH**2-(M1-M2)**2)
          DWID=DWID+SH(NUMH)**2*(MH**2-(M1+M2)**2)
          DWID=DWID*G2*SQRT(LAMB)/(16.0*PI*(MH**3))
          WID=DWID
          CALL SSSAVE(IDHHA,WID,ISW1,-ISW2,0,0,0)
          CALL SSSAVE(IDHHA,WID,-ISW1,ISW2,0,0,0)
        ENDIF
100   CONTINUE
C
      RETURN
      END
