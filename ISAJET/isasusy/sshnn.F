#include "isajet/pilot.h"
      SUBROUTINE SSHNN
C-----------------------------------------------------------------------
C     Calculates the decay widths of all neutral Higgses into all
C     possible pairs of neutralinos, and the decay widths of the
C     charged Higgs into any neutralino and any chargino
C
C     Bisset's NEUINO
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sspar.inc"
#include "isajet/sssm.inc"
#include "isajet/sstype.inc"
C
      DOUBLE PRECISION XIJ,XJI,DIJ,TEMP,DWZN,TEMP2,T2,RWZ,SWZ
     $,PI,SR2,XM,THETX,YM,THETY,SGL,CGL,SGR,CGR,MW1,MW2,THETM,THETP
     $,G2,GP2,BETA,ALPHA,MH,M1,M2
      DOUBLE PRECISION SN1,SN2,DWID,LAMB
      DOUBLE PRECISION A(4,4),MHI(3)
      DOUBLE PRECISION SSDLAM
      REAL WID
      INTEGER II,NUMH,I1,I2,IZ,IW,ID1,ID2,IDHHA
      INTEGER IDHI(3),IDZI(4)
C
C          Mass matrix parameters
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
C          The following was calculated in Bisset's MASZIN
      DO 10 II=1,4
         TEMP=SQRT(G2)*ZMIXSS(3,II)+SQRT(GP2)*ZMIXSS(4,II)
         TEMP=TEMP/SR2
         A(1,II)=-TEMP*SGR-SQRT(G2)*ZMIXSS(1,II)*CGR
         A(2,II)=TEMP*CGR-SQRT(G2)*ZMIXSS(1,II)*SGR
         A(3,II)=-TEMP*SGL+SQRT(G2)*ZMIXSS(2,II)*CGL
         A(4,II)=TEMP*CGL+SQRT(G2)*ZMIXSS(2,II)*SGL
10    CONTINUE
C
C          Arrays for loops
C
      MHI(1)=AMHL
      MHI(2)=AMHH
      MHI(3)=AMHA
      IDHI(1)=ISHL
      IDHI(2)=ISHH
      IDHI(3)=ISHA
      IDZI(1)=ISZ1
      IDZI(2)=ISZ2
      IDZI(3)=ISZ3
      IDZI(4)=ISZ4
C
C          Loop over neutral Higgs decays h(numh) into neutralino
C          pairs zi(i1) and zi(i2)
C
      DO 100 NUMH=1,3
        MH=MHI(NUMH)
        IDHHA=IDHI(NUMH)
        DO 110 I1=1,4
          M1=ABS(AMZISS(I1))
          ID1=IDZI(I1)
          DO 120 I2=I1,4
            M2=ABS(AMZISS(I2))
            ID2=IDZI(I2)
            IF(M1+M2.GE.MH) GO TO 120
            LAMB=SSDLAM(MH**2,M1**2,M2**2)
            IF(I1.EQ.I2) THEN
              DIJ = 0.5
            ELSE
              DIJ = 1.0
            ENDIF
            TEMP=-0.5*SIGN(1.,AMZISS(I1))*SIGN(1.,AMZISS(I2))
            XIJ=TEMP*(SQRT(G2)*ZMIXSS(3,I2)-SQRT(GP2)*ZMIXSS(4,I2))
            XJI=TEMP*(SQRT(G2)*ZMIXSS(3,I1)-SQRT(GP2)*ZMIXSS(4,I1))
            IF(NUMH.EQ.1) THEN
              XIJ=XIJ*(ZMIXSS(2,I1)*SIN(ALPHA)-ZMIXSS(1,I1)*COS(ALPHA))
              XJI=XJI*(ZMIXSS(2,I2)*SIN(ALPHA)-ZMIXSS(1,I2)*COS(ALPHA))
            ELSEIF (NUMH .EQ. 2) THEN
              XIJ=XIJ*(ZMIXSS(2,I1)*COS(ALPHA)+ZMIXSS(1,I1)*SIN(ALPHA))
              XJI=XJI*(ZMIXSS(2,I2)*COS(ALPHA)+ZMIXSS(1,I2)*SIN(ALPHA))
            ELSEIF(NUMH.EQ.3) THEN
              XIJ=-XIJ*(ZMIXSS(2,I1)*SIN(BETA)-ZMIXSS(1,I1)*COS(BETA))
              XJI=-XJI*(ZMIXSS(2,I2)*SIN(BETA)-ZMIXSS(1,I2)*COS(BETA))
            ENDIF
            DWID=DIJ*(XIJ+XJI)**2
            DWID=DWID*SQRT(LAMB)/(8.0*PI*(MH**3))
            IF(NUMH.EQ.1.OR.NUMH.EQ.2) THEN
              TEMP2 = ((MH**2)-(M1-2.0*TEMP*M2)**2)
            ELSEIF(NUMH.EQ.3) THEN
              TEMP2=((MH**2)-(M1+2.0*TEMP*M2)**2)
            ENDIF
            DWID=DWID*TEMP2
            WID=DWID
            CALL SSSAVE(IDHHA,WID,ID1,ID2,0,0,0)
120       CONTINUE
110     CONTINUE
100   CONTINUE
C
C          Loop over h+ decays into wi(iw) + zi(iz)
C
      MH=AMHC
      DO 210 IW=1,2
        IF(IW.EQ.1) THEN
          M1=ABS(AMW1SS)
          ID1=ISW1
          SN1=SIGN(1.,AMW1SS)
        ELSE
          M1=ABS(AMW2SS)
          ID1=ISW2
          SN1=SIGN(1.,AMW2SS)
        ENDIF
        DO 220 IZ=1,4
          M2=ABS(AMZISS(IZ))
          ID2=IDZI(IZ)
          SN2=SIGN(1.,AMZISS(IZ))
          IF(M1+M2.GE.MH) GO TO 220
          LAMB=SSDLAM(MH**2,M1**2,M2**2)
          T2=MH**2-M1**2-M2**2
          IF(IW.EQ.1) THEN
            RWZ=COS(BETA)*A(2,IZ)*SN1
            TEMP=SIN(BETA)*A(4,IZ)*SN2
            SWZ=0.5*(RWZ+TEMP)
            RWZ=0.5*(RWZ-TEMP)
          ELSE
            RWZ=COS(BETA)*A(1,IZ)*THETY*SN1
            TEMP=SIN(BETA)*A(3,IZ)*THETX*SN2
            SWZ=0.5*(RWZ+TEMP)
            RWZ=0.5*(RWZ-TEMP)
          ENDIF
          DWID=RWZ**2+SWZ**2
          DWID=DWID*T2
          TEMP=2*M1*M2*(RWZ**2-SWZ**2)
          DWID=(DWID-TEMP)/(8.0*PI*(MH**3))
          DWID=DWID*SQRT(LAMB)
          WID=DWID
          CALL SSSAVE(ISHC,WID,ID1,ID2,0,0,0)
220     CONTINUE
210   CONTINUE
      RETURN
      END
