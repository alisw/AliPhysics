#include "isajet/pilot.h"
      SUBROUTINE SSHGL
C-----------------------------------------------------------------------
C
C     Calculate H -> gl gl decays including both SM particles and
C     SUSY particles in loop.
C
C     This subroutine uses the tau variable of the Higgs Hunters'
C     Guide.  Many other authors, including the paper cited in 
C     Higgs Hunters' Guide (PR. D. 38(11): 3481) and Collider Physics
C     by Barger and Phillips use the variable lambda
C          LAMBDA = ( MASS OF PARTICLE IN LOOP / MASS OF HIGGS )**2
C          TAU = 4.0 * LAMBDA 
C
C     Bisset's HGLGL
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstype.inc"
C
      DOUBLE PRECISION ETAH,IITOT,RITOT,TAU,IFFF,RFFF
     $,IFHALF,RFHALF,IF1,RF1,IF0,RF0,TW2,RHF,RHSF,RHSFL,RHSFR
     $,IIHF,RIHF,IIHSFL,RIHSFL,IIHSFR,RIHSFR,AS,SUMISQ,DW
     $,RHSF1,RHSF2,IIHSF1,IIHSF2,RIHSF1,RIHSF2
      DOUBLE PRECISION PI,SR2,XM,THETX,YM,THETY,SGL,CGL,SGR,CGR
     $,MW1,MW2,THETM,THETP,G2,BETA,ALPHA,SW2,CW2,MH,AMSQ
      DOUBLE PRECISION MFL(3),MFD(3),MFU(3)
      DOUBLE PRECISION SSALFS
      REAL WID
      REAL ASMB,MBMB,MBQ,ASMT,MTMT,MTQ,SUALFS
      DOUBLE PRECISION SSMQCD
      INTEGER IJ,II,NUMOUT,NUMH,IDHHA
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
      G2=4.0*PI*ALFAEM/SN2THW
      BETA=ATAN(1.0/RV2V1)
      ALPHA=ALFAH
      SW2=SN2THW
      CW2=1.-SN2THW
C
C          Loop over neutral Higgs bosons
C
      DO 100 NUMH=1,3
        IF(NUMH.EQ.1) THEN
          MH=AMHL
          IDHHA=ISHL
        ELSEIF(NUMH.EQ.2) THEN
          MH=AMHH
          IDHHA=ISHH
        ELSE
          MH=AMHA
          IDHHA=ISHA
        ENDIF
        ETAH=1.0
        IITOT=0.0
        RITOT=0.0
C
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(MH))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(MH))

C
      MFL(1)=DBLE(AME)
      MFL(2)=DBLE(AMMU)
      MFL(3)=DBLE(AMTAU)
      MFD(1)=DBLE(AMDN)
      MFD(2)=DBLE(AMST)
      MFD(3)=DBLE(MBQ)
      MFU(1)=DBLE(AMUP)
      MFU(2)=DBLE(AMCH)
      MFU(3)=DBLE(MTQ)
C
C
C          Down-type quark loops
C
        DO 20 II=1,3
          TAU=4.0*MFD(II)**2/MH**2                  
          CALL SSHGM1(TAU,IFFF,RFFF)         
          IFHALF=0.5*TAU*(1.0-TAU*ETAH)*IFFF
          RFHALF=0.5*TAU*(ETAH+(1.0-TAU*ETAH)*RFFF)
          IF(NUMH.EQ.1) THEN
            RHF=SIN(ALPHA)/COS(BETA)
          ELSEIF(NUMH.EQ.2) THEN
            RHF=COS(ALPHA)/COS(BETA)
          ELSE
            RHF=TAN(BETA)
          ENDIF
          IIHF=RHF*IFHALF
          RIHF=RHF*RFHALF
          IITOT=IITOT+IIHF
          RITOT=RITOT+RIHF
20      CONTINUE 
C
C          Up-type quark loops
C
        DO 30 II=1,2
          TAU=4.0*MFU(II)**2/MH**2                  
          CALL SSHGM1(TAU,IFFF,RFFF)         
          IFHALF=0.5*TAU*(1.0-TAU*ETAH)*IFFF
          RFHALF=0.5*TAU*(ETAH+(1.0-TAU*ETAH)*RFFF)
          IF(NUMH.EQ.1) THEN
            RHF=COS(ALPHA)/SIN(BETA)
          ELSEIF(NUMH.EQ.2) THEN
            RHF=-SIN(ALPHA)/SIN(BETA)
          ELSE
            RHF=TAN(BETA)
          ENDIF
          IIHF=RHF*IFHALF
          RIHF=RHF*RFHALF
          IITOT=IITOT+IIHF
          RITOT=RITOT+RIHF
30      CONTINUE 
C
        TAU=4.0*MTQ**2/MH**2                  
        CALL SSHGM1(TAU,IFFF,RFFF)         
        IFHALF=0.5*TAU*(1.0-TAU*ETAH)*IFFF
        RFHALF=0.5*TAU*(ETAH+(1.0-TAU*ETAH)*RFFF)
        IF(NUMH.EQ.1) THEN
          RHF=COS(ALPHA)/SIN(BETA)
        ELSEIF(NUMH.EQ.2) THEN
          RHF=-SIN(ALPHA)/SIN(BETA)
        ELSE
          RHF=1.0/TAN(BETA)
        ENDIF
        IIHF=RHF*IFHALF
        RIHF=RHF*RFHALF
        IITOT=IITOT+IIHF
        RITOT=RITOT+RIHF
C
C          Down-type squark loops
C          Mixing between the sbottom squarks is included, so  
C          masses used here are the mixed masses (AMB1SS & AMB2SS)
C          First do d_L and s_L squarks
        TW2=SW2/CW2
        DO 50 II=1,2
          IF(NUMH.EQ.1) THEN
            RHSF=2.0*(MFD(II)/AMW)**2*SIN(ALPHA)/COS(BETA)
            RHSFL=(-1.0-TW2/3.0)*SIN(BETA-ALPHA)-RHSF
          ELSEIF(NUMH.EQ.2) THEN
            RHSF=2.0*(MFD(II)/AMW)**2*COS(ALPHA)/COS(BETA)
            RHSFL=(-1.0-TW2/3.0)*(-1.0)*COS(BETA-ALPHA)-RHSF
          ELSE
            RHSF=0
            RHSFL=0
          ENDIF
          IF (II.EQ.1) AMSQ=AMDLSS
          IF (II.EQ.2) AMSQ=AMSLSS
          TAU=4.0*AMSQ**2/MH**2                  
          CALL SSHGM1(TAU,IFFF,RFFF)         
          IF0=-TAU*TAU*IFFF          
          RF0=TAU*(1.0-TAU*RFFF)          
          IIHSFL=RHSFL*IF0*(AMW/AMSQ)**2/8.0
          RIHSFL=RHSFL*RF0*(AMW/AMSQ)**2/8.0
          IITOT=IITOT+IIHSFL
          RITOT=RITOT+RIHSFL
50      CONTINUE
c       Next, do R squarks
        DO 51 II=1,2
          IF(NUMH.EQ.1) THEN
            RHSF=2.0*(MFD(II)/AMW)**2*SIN(ALPHA)/COS(BETA)
            RHSFR=(-2.0*TW2/3.0)*SIN(BETA-ALPHA)-RHSF
          ELSEIF(NUMH.EQ.2) THEN
            RHSF=2.0*(MFD(II)/AMW)**2*COS(ALPHA)/COS(BETA)
            RHSFR=(-2.0*TW2/3.0)*(-1.0)*COS(BETA-ALPHA)-RHSF
          ELSE
            RHSF=0
            RHSFR=0
          ENDIF
          IF (II.EQ.1) AMSQ=AMDRSS
          IF (II.EQ.2) AMSQ=AMSRSS
          TAU=4.0*AMSQ**2/MH**2                  
          CALL SSHGM1(TAU,IFFF,RFFF)         
          IF0=-TAU*TAU*IFFF          
          RF0=TAU*(1.0-TAU*RFFF)          
          IIHSFR=RHSFR*IF0*(AMW/AMSQ)**2/8.0
          RIHSFR=RHSFR*RF0*(AMW/AMSQ)**2/8.0
          IITOT=IITOT+IIHSFR
          RITOT=RITOT+RIHSFR
51      CONTINUE
        IF(NUMH.EQ.1) THEN
          RHSF=2.0*(MBQ/AMW)**2*SIN(ALPHA)/COS(BETA)
          RHSFL=(-1.0-TW2/3.0)*SIN(BETA-ALPHA)-RHSF
          RHSFR=(-2.0*TW2/3.0)*SIN(BETA-ALPHA)-RHSF
        ELSEIF(NUMH.EQ.2) THEN
          RHSF=2.0*(MBQ/AMW)**2*COS(ALPHA)/COS(BETA)
          RHSFL=(-1.0-TW2/3.0)*(-1.0)*COS(BETA-ALPHA)-RHSF
          RHSFR=(-2.0*TW2/3.0)*(-1.0)*COS(BETA-ALPHA)-RHSF
        ELSE
          RHSF=0
          RHSFL=0
          RHSFR=0
        ENDIF
        RHSF1=RHSFL*COS(THETAB)-RHSFR*SIN(THETAB)
        RHSF2=RHSFL*SIN(THETAB)+RHSFR*COS(THETAB)
        TAU=4.0*AMB1SS**2/MH**2 
        CALL SSHGM1(TAU,IFFF,RFFF)         
        IF0=-TAU*TAU*IFFF          
        RF0=TAU*(1.0-TAU*RFFF)    
        IIHSF1=RHSF1*IF0*(AMW/AMB1SS)**2/8.0
        RIHSF1=RHSF1*RF0*(AMW/AMB1SS)**2/8.0
        IITOT=IITOT+IIHSF1
        RITOT=RITOT+RIHSF1
        TAU=4.0*AMB2SS**2/MH**2 
        CALL SSHGM1(TAU,IFFF,RFFF)         
        IF0=-TAU*TAU*IFFF          
        RF0=TAU*(1.0-TAU*RFFF)    
        IIHSF2=RHSF2*IF0*(AMW/AMB2SS)**2/8.0
        RIHSF2=RHSF2*RF0*(AMW/AMB2SS)**2/8.0
        IITOT=IITOT+IIHSF2
        RITOT=RITOT+RIHSF2
C
C          Up-type squark loops
C          Mixing between the stop squarks is included, so  
C          masses used here are the mixed masses (AMT1SS & AMT2SS)
C          First do u_L and c_L
        DO 60 II=1,2
          IF(NUMH.EQ.1) THEN
            RHSF=2.0*(MFU(II)/AMW)**2*COS(ALPHA)/SIN(BETA)
            RHSFL=(1.0-TW2/3.0)*SIN(BETA-ALPHA)-RHSF
          ELSEIF(NUMH.EQ.2) THEN
            RHSF=2.0*(MFU(II)/AMW)**2
            RHSF=RHSF*(-1.0)*SIN(ALPHA)/SIN(BETA)
            RHSFL=(1.0-TW2/3.0)*(-1.0)*COS(BETA-ALPHA)-RHSF
          ELSE
            RHSF=0
            RHSFL=0
          ENDIF
          IF (II.EQ.1) AMSQ=AMULSS
          IF (II.EQ.2) AMSQ=AMCLSS
          TAU=4.0*(AMSQ)**2/MH**2                  
          CALL SSHGM1(TAU,IFFF,RFFF)         
          IF0=-TAU*TAU*IFFF          
          RF0=TAU*(1.0-TAU*RFFF)          
          IIHSFL=RHSFL*IF0*(AMW/AMSQ)**2/8.0
          RIHSFL=RHSFL*RF0*(AMW/AMSQ)**2/8.0
          IITOT=IITOT+IIHSFL
          RITOT=RITOT+RIHSFL
60      CONTINUE
C          Next, do u_R and c_R
        DO 61 II=1,2
          IF(NUMH.EQ.1) THEN
            RHSF=2.0*(MFU(II)/AMW)**2*COS(ALPHA)/SIN(BETA)
            RHSFR=(4.0*TW2/3.0)*SIN(BETA-ALPHA)-RHSF
          ELSEIF(NUMH.EQ.2) THEN
            RHSF=2.0*(MFU(II)/AMW)**2
            RHSF=RHSF*(-1.0)*SIN(ALPHA)/SIN(BETA)
            RHSFR=(4.0*TW2/3.0)*(-1.0)*COS(BETA-ALPHA)-RHSF
          ELSE
            RHSF=0
            RHSFR=0
          ENDIF
          IF (II.EQ.1) AMSQ=AMURSS
          IF (II.EQ.2) AMSQ=AMCRSS
          TAU=4.0*(AMSQ)**2/MH**2                  
          CALL SSHGM1(TAU,IFFF,RFFF)         
          IF0=-TAU*TAU*IFFF          
          RF0=TAU*(1.0-TAU*RFFF)          
          IIHSFR=RHSFR*IF0*(AMW/AMSQ)**2/8.0
          RIHSFR=RHSFR*RF0*(AMW/AMSQ)**2/8.0
          IITOT=IITOT+IIHSFR
          RITOT=RITOT+RIHSFR
61      CONTINUE
C
        IF(NUMH.EQ.1) THEN
          RHSF=2.0*(MTQ/AMW)**2*COS(ALPHA)/SIN(BETA)
          RHSFL=(1.0-TW2/3.0)*SIN(BETA-ALPHA)-RHSF
          RHSFR=(4.0*TW2/3.0)*SIN(BETA-ALPHA)-RHSF
        ELSEIF(NUMH.EQ.2) THEN
          RHSF=2.0*(MTQ/AMW)**2
          RHSF=RHSF*(-1.0)*SIN(ALPHA)/SIN(BETA)
          RHSFL=(1.0-TW2/3.0)*(-1.0)*COS(BETA-ALPHA)-RHSF
          RHSFR=(4.0*TW2/3.0)*(-1.0)*COS(BETA-ALPHA)-RHSF
        ELSE
          RHSF=0
          RHSFL=0
          RHSFR=0
        ENDIF
        RHSF1=RHSFL*COS(THETAT)-RHSFR*SIN(THETAT)
        RHSF2=RHSFL*SIN(THETAT)+RHSFR*COS(THETAT)
        TAU=4.0*AMT1SS**2/MH**2 
        CALL SSHGM1(TAU,IFFF,RFFF)         
        IF0=-TAU*TAU*IFFF          
        RF0=TAU*(1.0-TAU*RFFF)    
        IIHSF1=RHSF1*IF0*(AMW/AMT1SS)**2/8.0
        RIHSF1=RHSF1*RF0*(AMW/AMT1SS)**2/8.0
        IITOT=IITOT+IIHSF1
        RITOT=RITOT+RIHSF1
        TAU=4.0*AMT2SS**2/MH**2 
        CALL SSHGM1(TAU,IFFF,RFFF)         
        IF0=-TAU*TAU*IFFF          
        RF0=TAU*(1.0-TAU*RFFF)    
        IIHSF2=RHSF2*IF0*(AMW/AMT2SS)**2/8.0
        RIHSF2=RHSF2*RF0*(AMW/AMT2SS)**2/8.0
        IITOT=IITOT+IIHSF2
        RITOT=RITOT+RIHSF2
C
C           IITOT and RITOT now contain the total imaginary and
C           real parts of the I function
C
        SUMISQ=IITOT**2+RITOT**2
        AS=SSALFS(MH**2)
        DW=AS**2*G2*MH**3/(32.0*(PI**3)*AMW**2) 
        WID=DW*SUMISQ
        CALL SSSAVE(IDHHA,WID,IDGL,IDGL,0,0,0)
100   CONTINUE
C
       RETURN
       END
