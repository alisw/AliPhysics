*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_OLDNAME)
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*=== decay ============================================================*
*
      SUBROUTINE DECAY(IHAD,ISTAB)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*  Decay89: slight revision by A. Ferrari                              *
*----------------------------------------------------------------------*
*
#include "geant321/finpar2.inc"
#include "geant321/metlsp.inc"
#include "geant321/part.inc"
#include "geant321/decayc.inc"
      COMMON /FKDREI/ TEST(12)
      COMMON /FKGAMR/REDU,AMO,AMM(15)
      COMMON /FKPRUN/ISYS
      REAL RNDM(1)
C
C
      REDU=2.D0
      DO 801 I=1,IHAD
         ITS(I) = NREF(I)
         PLS(I) = SQRT(PXF(I)**2 + PYF(I)**2 + PZF(I)**2)
         IF (PLS(I) .NE. 0.D0) THEN
            CXS(I) = PXF(I)/PLS(I)
            CYS(I) = PYF(I)/PLS(I)
            CZS(I) = PZF(I)/PLS(I)
         END IF
         ELS(I) = HEF(I)
  801 CONTINUE
      IST = IHAD
      IR  = 0
  200 CONTINUE
C*****TEST STABLE OR UNSTABLE
C$$$$$ISTAB=1/2/3 MEANS  STRONG + WEAK DECAYS / ONLY STRONG DECAYS /
C*****STRONG DECAYS + WEAK DECAYS FOR CHARMED PARTICLES AND TAU LEPTONS
      IF(ISTAB.EQ.1) GOTO 793
      IF(ISTAB.EQ.2) GOTO 737
      IF(ISTAB.EQ.3) GOTO 738
  793 IF(ITS(IST).EQ.135.OR.ITS(IST).EQ.136) GOTO 202
      IF(ITS(IST).GE.1.AND.ITS(IST).LE.7) GOTO 202
      GOTO 300
  738 IF(ITS(IST).GE.1.AND.ITS(IST).LE.23) GOTO 202
      IF(ITS(IST).GE. 97.AND.ITS(IST).LE.103) GOTO 202
      IF(ITS(IST).EQ.109.AND.ITS(IST).EQ.115) GOTO 202
      IF(ITS(IST).GE.133.AND.ITS(IST).LE.136) GOTO 202
      GOTO 300
  737 IF(ITS(IST).GE.  1.AND.ITS(IST).LE. 30) GOTO 202
      IF(ITS(IST).GE. 97.AND.ITS(IST).LE.103) GOTO 202
      IF(ITS(IST).GE.115.AND.ITS(IST).LE.122) GOTO 202
      IF(ITS(IST).GE.131.AND.ITS(IST).LE.136) GOTO 202
      IF(ITS(IST).EQ.109) GO TO 202
      IF(ITS(IST).GE.137.AND.ITS(IST).LE.160) GOTO 202
      GO TO 300
  202 IR = IR + 1
      NREF(IR) = ITS(IST)
      ITT = ITS(IST)
      AMF(IR) = AM(ITT)
      ANF(IR) = ANAME(ITT)
      ICHF(IR)  = ICH(ITT)
      IBARF(IR) = IBAR(ITT)
      HEF(IR) = ELS(IST)
      PXF(IR) = CXS(IST)*PLS(IST)
      PYF(IR) = CYS(IST)*PLS(IST)
      PZF(IR) = CZS(IST)*PLS(IST)
      IST = IST - 1
      IF(IST .GE. 1) GO TO 200
      GO TO 500
  300 IT = ITS(IST)
      GAM = ELS(IST)/AM(IT)
      BGAM = PLS(IST)/AM(IT)
      ECO = AM(IT)
      KZ1 = K1(IT)
  310 CONTINUE
      CALL GRNDM(RNDM,1)
      VV = RNDM(1) - 1.D-17
      IIK = KZ1 - 1
  301 IIK = IIK + 1
      IF (VV.GT.WT(IIK)) GO TO 301
C  IIK IS THE DECAY CHANNEL
      IT1 = NZK(IIK,1)
      IT2 = NZK(IIK,2)
      IF (IT2-1 .LT. 0) GO TO 110
      IT3 = NZK(IIK,3)
C  IT1,IT2, IT3 ARE THE PRODUCED PARTICLES FROM  IT
      IF(IT3 .EQ. 0) GO TO 400
      CALL THREPD(ECO,ECM1,ECM2,ECM3,PCM1,PCM2,PCM3,COD1,COF1,SIF1,COD2,
     &            COF2,SIF2,COD3,COF3,SIF3,AM(IT1),AM(IT2),AM(IT3))
      GO TO 411
  400 CALL TWOPAD(ECO,ECM1,ECM2,PCM1,PCM2,COD1,COF1,SIF1,COD2,COF2,SIF2,
     &            AM(IT1),AM(IT2))
  411 CONTINUE
  110 CONTINUE
      ITS(IST) = IT1
      IF (IT2-1 .LT. 0) GO TO 305
      ITS(IST+1) = IT2
      ITS(IST+2) = IT3
      RX = CXS(IST)
      RY = CYS(IST)
      RZ = CZS(IST)
      CALL TRAFO(GAM,BGAM,RX,RY,RZ,COD1,COF1,SIF1,PCM1,ECM1,
     &           PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IST = IST + 1
      CALL TRAFO(GAM,BGAM,RX,RY,RZ,COD2,COF2,SIF2,PCM2,ECM2,
     &           PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IF (IT3 .LE. 0) GO TO 305
      IST = IST + 1
      CALL TRAFO(GAM,BGAM,RX,RY,RZ,COD3,COF3,SIF3,PCM3,ECM3,
     &           PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
  305 CONTINUE
      GO TO 200
  500 CONTINUE
      IDAPU = IDMAX3
      IF(IR .GT. IDMAX3) WRITE(ISYS,928)IDAPU
  928 FORMAT(' NUMBER OF STAB. FINAL PART. IS GREATER THAN',I5)
      IHAD = IR
      RETURN
      END
#endif
