#include "isajet/pilot.h"
        SUBROUTINE SSSTBF
C-----------------------------------------------------------------------
C
C        This program gives stop squark branching fractions to gauginos
C        according to Baer and Tata.
C        If no other modes are allowed, stop -> c z_i through loops is
C        used as the default.
C
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/ssmode.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstype.inc"
C
      COMPLEX ZI,ZONE,ZA,ZB,ZPP,ZPM,ZAUIZ,ZBUIZ
      DOUBLE PRECISION SSALFS,SSMQCD
      REAL SSXLAM
      REAL WID,AWD(2),BW(2),FB,FT,XM,YM,THX,THY,AU1,MZ1,WIDC1
      REAL PI,SR2,G,GP,TANB,COTB,MPL,MMI,AH
      REAL AUIZ,MZIZ,SINT,COST,AS,BS,SNZI,THIZ
      INTEGER IZ,ISTOP,IDSTOP
      REAL AMSTOP,BWP(2),A
      REAL MW1,MW2,SNW1,SNW2,CS2THW,BETA,TN2THW,SINB,COSB
      REAL ASMB,MBMB,MBQ,ASMT,MTMT,MTQ,SUALFS
      INTEGER ISZIZ(4)
      DATA ZONE/(1.,0.)/,ZI/(0.,1.)/
C
C          Partly duplicated from SSMASS
C
      CS2THW=1.-SN2THW
      TN2THW=SN2THW/CS2THW
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      TANB=1./RV2V1
      COTB=RV2V1
      BETA=ATAN(TANB)
C          Reconstruct masses from SSMASS
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(AMT1SS))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(AMT1SS))
      FB=G*MBQ/SR2/AMW/COS(BETA)
      FT=G*MTQ/SR2/AMW/SIN(BETA)
      MW1=ABS(AMW1SS)
      MW2=ABS(AMW2SS)
      SNW1=SIGN(1.,AMW1SS)
      SNW2=SIGN(1.,AMW2SS)
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
C
      AWD(1)=-G*SNW1*SIN(GAMMAR)
      AWD(2)=-G*SNW2*THY*COS(GAMMAR)
      BW(1)=-FT*SNW1*COS(GAMMAR)
      BW(2)=FT*SNW2*THY*SIN(GAMMAR)
      BWP(1)=-FB*COS(GAMMAL)
      BWP(2)=FB*THX*SIN(GAMMAL)
      MMI=AMW1SS
      MPL=AMW2SS
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
C
C          Compute stop_i branching fractions to charm + zi if no other
C          modes are allowed. WIDC1 is an unknown amplitude arbitrarily
C          set equal to 1.0.
C
      ISZIZ(1)=ISZ1
      ISZIZ(2)=ISZ2
      ISZIZ(3)=ISZ3
      ISZIZ(4)=ISZ4
      AU1=-G/SR2*ZMIXSS(3,1)-GP/3./SR2*ZMIXSS(4,1)
      MZ1=ABS(AMZISS(1))
      DO 100 ISTOP=1,2
        IF(ISTOP.EQ.1) THEN
          AMSTOP=AMT1SS
          IDSTOP=ISTP1
          WIDC1=1.0E-6
        ELSE
          AMSTOP=AMT2SS
          IDSTOP=ISTP2
          WIDC1=1.0E-6
        ENDIF
        IF(AMSTOP.LT.(MW1+AMBT).AND.AMSTOP.GT.(AMCH+MZ1)) THEN
          DO 110 IZ=1,4
            MZIZ=ABS(AMZISS(IZ))
            AUIZ=-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ)
            IF (AMT1SS.GT.(AMCH+MZIZ)) THEN
              WID=AUIZ**2*(AMSTOP**2-MZIZ**2)/AU1**2
     $        /(AMSTOP**2-MZ1**2)*WIDC1
              CALL SSSAVE(IDSTOP,WID,ISZIZ(IZ),IDCH,0,0,0)
            END IF
110       CONTINUE
        ELSEIF(AMSTOP.LT.(MW1+AMBT).AND.AMSTOP.LE.(AMCH+MZ1)) THEN
          WRITE(LOUT,1000) ISTOP
1000      FORMAT(' ERROR IN SSSTBF: NO ALLOWED MODE FOR STOP',I2)
        END IF
100   CONTINUE
C
C          stop_i -> gluino + top
C
      IF (AMT1SS.GT.(AMGLSS+AMTP)) THEN
        WID=2*SSALFS(DBLE(AMT1SS**2))*AMT1SS*((1.-AMGLSS**2/AMT1SS**2-
     $  AMTP**2/AMT1SS**2)-2*SIN(2*THETAT)*AMTP*AMGLSS/AMT1SS**2)
     $  *SQRT(SSXLAM(1.,AMGLSS**2/AMT1SS**2,AMTP**2/AMT1SS**2))/3.
        CALL SSSAVE(ISTP1,WID,ISGL,IDTP,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMGLSS+AMTP)) THEN
        WID=2*SSALFS(DBLE(AMT2SS**2))*AMT2SS*((1.-AMGLSS**2/AMT2SS**2-
     $  AMTP**2/AMT2SS**2)+2*SIN(2*THETAT)*AMTP*AMGLSS/AMT2SS**2)
     $  *SQRT(SSXLAM(1.,AMGLSS**2/AMT2SS**2,AMTP**2/AMT2SS**2))/3.
        CALL SSSAVE(ISTP2,WID,ISGL,IDTP,0,0,0)
      END IF
C
C          stop_1 -> top + zino_i
C
      DO 200 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        SNZI=SIGN(1.,AMZISS(IZ))
        IF (SNZI.EQ.1.) THEN
           THIZ=0.
        ELSE
           THIZ=1.
        END IF
        ZAUIZ=ZI**(THIZ-1.)*SNZI
     $  *(-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBUIZ=ZI**(THIZ-1.)*4*GP*ZMIXSS(4,IZ)/3./SR2
        ZPP=ZI**THIZ
        ZPM=(-ZI)**THIZ
        ZA=((ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*COST
     $  -(ZI*ZBUIZ-ZPM*FT*ZMIXSS(1,IZ))*SINT)/2.
        ZB=((-ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*COST
     $  -(ZI*ZBUIZ+ZPM*FT*ZMIXSS(1,IZ))*SINT)/2.
        AS=ZA*CONJG(ZA)
        BS=ZB*CONJG(ZB)
        IF (AMT1SS.GT.(AMTP+MZIZ)) THEN
          WID=(AS*(AMT1SS**2-(AMTP+MZIZ)**2)+BS
     $    *(AMT1SS**2-(AMTP-MZIZ)**2))/8./PI/AMT1SS
     $    *SQRT(SSXLAM(1.,AMTP**2/AMT1SS**2,MZIZ**2/AMT1SS**2))
          CALL SSSAVE(ISTP1,WID,ISZIZ(IZ),IDTP,0,0,0)
        END IF
200   CONTINUE
C
C          Wino decays
C
      IF (AMT1SS.GT.(AMBT+MW1)) THEN
        A=AWD(1)*COST-BW(1)*SINT
        AS=A*A
        WID=AMT1SS*((AS+BWP(1)**2*COST**2)*(1.-MW1**2/AMT1SS**2-
     $   AMBT**2/AMT1SS**2)-4*MW1*AMBT*BWP(1)*COST*A/AMT1SS**2)
     $   *SQRT(SSXLAM(1.,MW1**2/AMT1SS**2,AMBT**2/AMT1SS**2))/16./PI
        CALL SSSAVE(ISTP1,WID,ISW1,IDBT,0,0,0)
      END IF
      IF (AMT1SS.GT.(AMBT+MW2)) THEN
        A=AWD(2)*COST-BW(2)*SINT
        AS=A*A
        WID=AMT1SS*((AS+BWP(2)**2*COST**2)*(1.-MW2**2/AMT1SS**2-
     $   AMBT**2/AMT1SS**2)-4*MW2*AMBT*BWP(2)*COST*A/AMT1SS**2)
     $   *SQRT(SSXLAM(1.,MW2**2/AMT1SS**2,AMBT**2/AMT1SS**2))/16./PI
        CALL SSSAVE(ISTP1,WID,ISW2,IDBT,0,0,0)
      END IF
C 
      IF (AMT2SS.GT.(AMBT+MW1)) THEN
        A=AWD(1)*SINT+BW(1)*COST
        AS=A*A
        WID=AMT2SS*((AS+BWP(1)**2*SINT**2)*(1.-MW1**2/AMT2SS**2-
     $   AMBT**2/AMT2SS**2)-4*MW1*AMBT*BWP(1)*SINT*A/AMT2SS**2)
     $   *SQRT(SSXLAM(1.,MW1**2/AMT2SS**2,AMBT**2/AMT2SS**2))/16./PI
        CALL SSSAVE(ISTP2,WID,ISW1,IDBT,0,0,0)
      END IF
      IF (AMT2SS.GT.(AMBT+MW2)) THEN
        A=AWD(2)*SINT+BW(2)*COST
        AS=A*A
        WID=AMT2SS*((AS+BWP(2)**2*SINT**2)*(1.-MW2**2/AMT2SS**2-
     $   AMBT**2/AMT2SS**2)-4*MW2*AMBT*BWP(2)*SINT*A/AMT2SS**2)
     $  *SQRT(SSXLAM(1.,MW2**2/AMT2SS**2,AMBT**2/AMT2SS**2))/16./PI
        CALL SSSAVE(ISTP2,WID,ISW2,IDBT,0,0,0)
      END IF
C
C          stop_2 -> stop_1 + X modes
C
      IF (AMT2SS.GT.(AMT1SS+AMZ)) THEN
        WID=G**2*COST**2*SINT**2
     $  *(SQRT(SSXLAM(AMT2SS**2,AMZ**2,AMT1SS**2)))**3
     $  /64./PI/CS2THW/AMT2SS**3/AMZ**2
        CALL SSSAVE(ISTP2,WID,IDZ,ISTP1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMT1SS+AMHL)) THEN
        AH=G*AMW*SIN(BETA-ALFAH)*(1.-5.*TN2THW/3.)*SINT*COST/2.
     $  +G*AMTP*COS(2.*THETAT)*(TWOM1*SIN(ALFAH)+AAT*COS(ALFAH))/2.
     $  /AMW/SIN(BETA)
        WID=AH**2/16./PI/AMT2SS**3
     $  *SQRT(SSXLAM(AMT2SS**2,AMHL**2,AMT1SS**2))
        CALL SSSAVE(ISTP2,WID,ISHL,ISTP1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMT1SS+AMHH)) THEN
        AH=-G*AMW*COS(BETA-ALFAH)*(1.-5.*TN2THW/3.)*SINT*COST/2.
     $  +G*AMTP*COS(2.*THETAT)*(TWOM1*COS(ALFAH)-AAT*SIN(ALFAH))/2.
     $  /AMW/SIN(BETA)
        WID=AH**2/16./PI/AMT2SS**3
     $  *SQRT(SSXLAM(AMT2SS**2,AMHH**2,AMT1SS**2))
        CALL SSSAVE(ISTP2,WID,ISHH,ISTP1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMT1SS+AMHA)) THEN
        AH=G*AMTP*(TWOM1-AAT/TANB)/2./AMW
        WID=AH**2/16./PI/AMT2SS**3
     $  *SQRT(SSXLAM(AMT2SS**2,AMHA**2,AMT1SS**2))
        CALL SSSAVE(ISTP2,WID,ISHA,ISTP1,0,0,0)
      END IF
C
C          t_i --> b_i + W decays
C
      IF (AMT1SS.GT.(AMB1SS+AMW)) THEN
        WID=G**2*COST**2*COSB**2*(SSXLAM(AMT1SS**2,AMB1SS**2,
     $AMW**2))**1.5/32./PI/AMT1SS**3/AMW**2
        CALL SSSAVE(ISTP1,WID,IDW,ISBT1,0,0,0)
      END IF
C
      IF (AMT1SS.GT.(AMB2SS+AMW)) THEN
        WID=G**2*COST**2*SINB**2*(SSXLAM(AMT1SS**2,AMB2SS**2,
     $AMW**2))**1.5/32./PI/AMT1SS**3/AMW**2
        CALL SSSAVE(ISTP1,WID,IDW,ISBT2,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMB1SS+AMW)) THEN
        WID=G**2*SINT**2*COSB**2*(SSXLAM(AMT2SS**2,AMB1SS**2,
     $AMW**2))**1.5/32./PI/AMT2SS**3/AMW**2
        CALL SSSAVE(ISTP2,WID,IDW,ISBT1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMB2SS+AMW)) THEN
        WID=G**2*SINT**2*SINB**2*(SSXLAM(AMT2SS**2,AMB2SS**2,
     $AMW**2))**1.5/32./PI/AMT2SS**3/AMW**2
        CALL SSSAVE(ISTP2,WID,IDW,ISBT2,0,0,0)
      END IF
C
C          t_i --> b_i + H+ decays
C
      IF (AMT1SS.GT.(AMB1SS+AMHC)) THEN
        A=G/SR2/AMW*(AMTP*AMBT*(COTB+TANB)*SINT*SINB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $COST*COSB-AMTP*(TWOM1-AAT*COTB)*SINT*COSB-AMBT*
     $(TWOM1-AAB*TANB)*SINB*COST)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMT1SS**2,AMB1SS**2,AMHC**2))/
     $      16./PI/AMT1SS**3
        CALL SSSAVE(ISTP1,WID,ISHC,ISBT1,0,0,0)
      END IF
C
      IF (AMT1SS.GT.(AMB2SS+AMHC)) THEN
        A=G/SR2/AMW*(-AMTP*AMBT*(COTB+TANB)*SINT*COSB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $COST*SINB-AMTP*(TWOM1-AAT*COTB)*SINT*SINB+AMBT*
     $(TWOM1-AAB*TANB)*COST*COSB)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMT1SS**2,AMB2SS**2,AMHC**2))/
     $      16./PI/AMT1SS**3
        CALL SSSAVE(ISTP1,WID,ISHC,ISBT2,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMB1SS+AMHC)) THEN
        A=G/SR2/AMW*(-AMTP*AMBT*(COTB+TANB)*COST*SINT+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $SINT*COSB+AMTP*(TWOM1-AAT*COTB)*COST*COSB-AMBT*
     $(TWOM1-AAB*TANB)*SINT*SINB)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMT2SS**2,AMB1SS**2,AMHC**2))/
     $      16./PI/AMT2SS**3
        CALL SSSAVE(ISTP2,WID,ISHC,ISBT1,0,0,0)
      END IF
C
      IF (AMT2SS.GT.(AMB2SS+AMHC)) THEN
        A=G/SR2/AMW*(AMTP*AMBT*(COTB+TANB)*COST*COSB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $SINT*SINB+AMTP*(TWOM1-AAT*COTB)*SINB*COST+AMBT*
     $(TWOM1-AAB*TANB)*COSB*SINT)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMT2SS**2,AMB2SS**2,AMHC**2))/
     $      16./PI/AMT2SS**3
        CALL SSSAVE(ISTP2,WID,ISHC,ISBT2,0,0,0)
      END IF
C
C
C          stop_2 -> top + zino_i
C
      DO 500 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        SNZI=SIGN(1.,AMZISS(IZ))
        IF (SNZI.EQ.1.) THEN
           THIZ=0.
        ELSE
           THIZ=1.
        END IF
        ZAUIZ=ZI**(THIZ-1.)*SNZI
     $  *(-G/SR2*ZMIXSS(3,IZ)-GP/3./SR2*ZMIXSS(4,IZ))
        ZBUIZ=ZI**(THIZ-1.)*4*GP*ZMIXSS(4,IZ)/3./SR2
        ZPP=ZI**THIZ
        ZPM=(-ZI)**THIZ
        ZA=((ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*SINT
     $  +(ZI*ZBUIZ-ZPM*FT*ZMIXSS(1,IZ))*COST)/2.
        ZB=((-ZI*ZAUIZ-ZPP*FT*ZMIXSS(1,IZ))*SINT
     $  +(ZI*ZBUIZ+ZPM*FT*ZMIXSS(1,IZ))*COST)/2.
        AS=ZA*CONJG(ZA)
        BS=ZB*CONJG(ZB)
        IF (AMT2SS.GT.(AMTP+MZIZ)) THEN
          WID=(AS*(AMT2SS**2-(AMTP+MZIZ)**2)+BS
     $    *(AMT2SS**2-(AMTP-MZIZ)**2))/8./PI/AMT2SS
     $    *SQRT(SSXLAM(1.,AMTP**2/AMT2SS**2,MZIZ**2/AMT2SS**2))
          CALL SSSAVE(ISTP2,WID,ISZIZ(IZ),IDTP,0,0,0)
        END IF
500   CONTINUE
C 
C          Normalize branching ratios 
C
       CALL SSNORM(ISTP1)
       CALL SSNORM(ISTP2)
C
       RETURN
       END
