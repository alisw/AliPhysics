#include "isajet/pilot.h"
      SUBROUTINE SSHFF
C-----------------------------------------------------------------------
C
C     Calculate all decays higgs -> f fbar, including QCD radiative
C     corrections for quarks.
C
C     Bisset's SETFAC, WDHFFN, QCDRAD
C
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstype.inc"
C
      DOUBLE PRECISION PI,SR2,G2,DWID,MHIH,BETA,BEFAC,ALFAC,MH,MF
     $,MFRUN,FACTOR,ALAM,MF1,MF2,SUM,MF1RUN,MF2RUN,COLOR,TEMP1
     $,QCDFAC
      DOUBLE PRECISION MFIFF(9),MFIF1(6),MFIF2(6)
      DOUBLE PRECISION SSDLAM,SSMQCD,SSHFF1
      REAL WID
      INTEGER IH,IDIH,IFF,IDF,ID1,ID2
      INTEGER IDIFF(9),IDIF1(6),IDIF2(6)
C
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      BETA=ATAN(1./RV2V1)
      G2=4.0*PI*ALFAEM/SN2THW
C
C          Loop over HL, HH, HA and fermions
C
      MFIFF(1)=AME
      IDIFF(1)=IDE
      MFIFF(2)=AMMU
      IDIFF(2)=IDMU
      MFIFF(3)=AMTAU
      IDIFF(3)=IDTAU
      MFIFF(4)=AMDN
      IDIFF(4)=IDDN
      MFIFF(5)=AMST
      IDIFF(5)=IDST
      MFIFF(6)=AMBT
      IDIFF(6)=IDBT
      MFIFF(7)=AMUP
      IDIFF(7)=IDUP
      MFIFF(8)=AMCH
      IDIFF(8)=IDCH
      MFIFF(9)=AMTP
      IDIFF(9)=IDTP
C
      DO 100 IH=1,3
        IF(IH.EQ.1) THEN
          MH=AMHL
          IDIH=ISHL
          BEFAC=COS(BETA)
          ALFAC=SIN(ALFAH)
        ELSEIF(IH.EQ.2) THEN
          MH=AMHH
          IDIH=ISHH
          BEFAC=COS(BETA)
          ALFAC=COS(ALFAH)
        ELSE
          MH=AMHA
          IDIH=ISHA
          BEFAC=1/TAN(BETA)
          ALFAC=1.
        ENDIF
C
C          Down type fermions
C
        DO 110 IFF=1,6
          MF=MFIFF(IFF)
          IDF=IDIFF(IFF)
          FACTOR=1.-4.*MF**2/MH**2
          IF(FACTOR.LE.0) GO TO 110
          FACTOR=SQRT(FACTOR)
          IF(IFF.GE.4) THEN
             COLOR=3.
             MFRUN=SSMQCD(MF,MH)
             QCDFAC=SSHFF1(MH,MF,IH)
          ELSE
             COLOR=1.
             MFRUN=MF
             QCDFAC=1.
          ENDIF
          DWID=G2*MFRUN**2*MH*ALFAC**2/(32.*PI*AMW**2*BEFAC**2)
          IF(IH.EQ.1.OR.IH.EQ.2) THEN
            DWID=DWID*FACTOR**3
          ELSEIF(IH.EQ.3) THEN
            DWID=DWID*FACTOR
          ENDIF
          DWID=DWID*COLOR*QCDFAC
          WID=DWID
          CALL SSSAVE(IDIH,WID,IDF,-IDF,0,0,0)
110     CONTINUE
C
C          Up type fermions
C
        IF(IH.EQ.1) THEN
          BEFAC=SIN(BETA)
          ALFAC=COS(ALFAH)
        ELSEIF(IH.EQ.2) THEN
          BEFAC=SIN(BETA)
          ALFAC=SIN(ALFAH)
        ELSE
          BEFAC=TAN(BETA)
          ALFAC=1.
        ENDIF
        DO 120 IFF=7,9
          MF=MFIFF(IFF)
          IDF=IDIFF(IFF)
          FACTOR=1.-4.*MF**2/MH**2
          IF(FACTOR.LE.0) GO TO 120
          FACTOR=SQRT(FACTOR)
          MFRUN=SSMQCD(MF,MH)
          QCDFAC=SSHFF1(MH,MF,IH)
          DWID=G2*MFRUN**2*MH*ALFAC**2/(32.*PI*AMW**2*BEFAC**2)
          IF(IH.EQ.1.OR.IH.EQ.2) THEN
            DWID=DWID*FACTOR**3
          ELSEIF(IH.EQ.3) THEN
            DWID=DWID*FACTOR
          ENDIF
          DWID=3.*DWID*QCDFAC
          WID=DWID
          CALL SSSAVE(IDIH,WID,IDF,-IDF,0,0,0)
120     CONTINUE
100   CONTINUE
C
C           HC decays. F1 has Iz=+1/2, F2 has Iz=-1/2
C
      MFIF1(1)=0
      IDIF1(1)=IDNE
      MFIF2(1)=AME
      IDIF2(1)=IDE
      MFIF1(2)=0
      IDIF1(2)=IDNM
      MFIF2(2)=AMMU
      IDIF2(2)=IDMU
      MFIF1(3)=0
      IDIF1(3)=IDNT
      MFIF2(3)=AMTAU
      IDIF2(3)=IDTAU
      MFIF1(4)=AMUP
      IDIF1(4)=IDUP
      MFIF2(4)=AMDN
      IDIF2(4)=IDDN
      MFIF1(5)=AMCH
      IDIF1(5)=IDCH
      MFIF2(5)=AMST
      IDIF2(5)=IDST
      MFIF1(6)=AMTP
      IDIF1(6)=IDTP
      MFIF2(6)=AMBT
      IDIF2(6)=IDBT
      MH=AMHC
C
      DO 200 IFF=1,6
        MF1=MFIF1(IFF)
        MF2=MFIF2(IFF)
        ID1=IDIF1(IFF)
        ID2=IDIF2(IFF)
        SUM=MF1+MF2
        ALAM=SSDLAM(MH**2,MF1**2,MF2**2)
        IF(ALAM.LE.0.OR.SUM.GE.MH) GO TO 200
        IF(IFF.LE.3) THEN
          MF1RUN=MF1
          MF2RUN=MF2
          COLOR=1
        ELSE
          MF1RUN=SSMQCD(MF1,MH)
          MF2RUN=SSMQCD(MF2,MH)
          COLOR=3
        ENDIF
        TEMP1=MF1RUN**2*1./TAN(BETA)**2+MF2RUN**2*TAN(BETA)**2
        TEMP1=TEMP1*(MH**2-MF1**2-MF2**2)-4.*MF1**2*MF2**2
        IF (TEMP1.LT.0.0) GO TO 200
        DWID=G2*COLOR*SQRT(ALAM)*TEMP1/MH**3/(32.0*PI*AMW**2)
        WID=DWID
        CALL SSSAVE(ISHC,WID,ID1,-ID2,0,0,0)
200   CONTINUE
C
      RETURN
      END
