*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE HADRIN.FOR
*COPY HADRIN
*
*=== hadrin ===========================================================*
*
*
*----------------------------------------------------------------------*
*    hadrin89: slight modifications by A. Ferrari
*----------------------------------------------------------------------*
*
      SUBROUTINE HADRIN (N,PLAB,ELAB,CX,CY,CZ,ITTA)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
#include "geant321/metlsp.inc"
#include "geant321/finlsp.inc"
      PARAMETER ( AMPROT = 0.93827231D+00 )
*
#include "geant321/reac.inc"
#include "geant321/redver.inc"
#include "geant321/split.inc"
*
****  INTEGER*2ICH,IBAR,K1,K2,NZK,NRK,IEII,IKII,NURE
      COMMON / FKGAMR / REDU, AMO, AMM (15)
C
      COMMON / FKABLT / AM   (110), GA   (110), TAU  (110), ICH   (110),
     &                  IBAR (110), K1   (110), K2   (110)
      COMMON / FKCODS / COD1, COD2, COD3, COF1, COF2, COF3, SIF1, SIF2,
     &                  SIF3, ECM1, ECM2, ECM3, PCM1, PCM2, PCM3
      COMMON / FKRUN    / RUNTES, EFTES
*
      SAVE UMODA, ITPRF, NNN
      DIMENSION ITPRF(110)
      REAL RNDM(1)
      DATA NNN/0/
      DATA UMODA/0.D0/
      DATA ITPRF/-1,-1,5*1,-1,-1,1,1,1,-1,-1,-1,-1,6*1,-1,-1,-1,85*1/
C
C-----------------------------
C*** INPUT VARIABLES LIST:
C*** SAMPLING OF HADRON NUCLEON INTERACTION FOR (ABOUT) 0.1 LE PLAB LE 6
C*** GEV/C LABORATORY MOMENTUM REGION
C*** N    - PROJECTILE HADRON INDEX
C*** PLAB - LABORATORY MOMENTUM OF N (GEV/C)
C*** ELAB - LABORATORY ENERGY OF N (GEV)
C*** CX,CY,CZ - DIRECTION COSINES OF N IN THE LABORATORY SYSTEM
C*** ITTA - TARGET NUCLEON INDEX
C*** OUTPUT VARIABLES LIST OF PARTICLE CHARACTERISTICS IN /FINLSP/
C  IR COUNTS THE NUMBER OF PRODUCED PARTICLES
C*** ITR - PARTICLE INDEX, CXR,CYR,CZR - DIRECTION COSINES (LAB. SYST.)
C*** ELR,PLR LAB. ENERGY AND LAB. MOMENTUM OF THE SAMPLED PARTICLE
C*** RESPECT., UNITS (GEV/C AND GEV)
C----------------------------
      LOWP=0
      IF (ITPRF(  N  ).LT.0) GO TO 99999
      WRITE(LUNOUT,99998) N
*     STOP    Commented out A. Fasso' 1989
      IR=0
      RETURN
99998 FORMAT (3(5H ****/),
     *45H FALSE USE OF THE PARTICLE TYPE INDEX, VALUE ,I4,3(5H ****/))
99999 CONTINUE
      IATMPT=0
         IF (ABS(PLAB-5.D0).GE.4.99999D0) THEN
            WRITE(LUNOUT,99996) PLAB
*           STOP     Commented out A. Fasso' 1989
            IR=0
            RETURN
99996       FORMAT (3(5H ****/),64
     *H PROJECTILE HADRON MOMENTUM OUTSIDE OF THE ALLOWED REGION, PLAB=,
     *1E15.5/,3(5H ****/))
      END IF
      UMODAT=N*1.11111D0+ITTA*2.19291D0
      IF (UMODAT.NE.UMODA.OR.(AMPROT-AM(1)).GT.1.D-6)CALL CALUMO(N,ITTA)
      UMODA=UMODAT
 1009 CONTINUE
      IATMPT=0
      LOWP=LOWP+1
 1000 CONTINUE
      IMACH=0
      REDU=2.D0
      IW1=0
      IF (LOWP.GT.20) GO TO 8
      NNN=N
      IF (NNN.EQ.N) GO TO 4322
      RUNTES=0.D0
      EFTES=0.D0
 4322 CONTINUE
      IS=1
      IR=0
      IST=1
      NSTAB=25
      IF ( ITTA .EQ. 1 ) THEN
         IRE=NURE(N,1)
      ELSE
         IRE=NURE(N,2)
      END IF
C
C-----------------------------
C*** IE,AMT,ECM,SI DETERMINATION
C----------------------------
      CALL FKSIGI(IRE,PLAB,N,IE,AMT,AMN,ECM,SI,ITTA)
      IF ( AMPROT - AM(1) .GT. 1.D-6 ) THEN
         IANTH = 1
         SI = 1.D0
      ELSE
         IANTH = -1
      END IF
      ECMMH=ECM
C
C-----------------------------
C    ENERGY INDEX
C  IRE CHARACTERIZES THE REACTION
C  IE IS THE ENERGY INDEX
C----------------------------
      IF (SI.LT.1.D-6) GO TO 8
      IF (N .LE.NSTAB) GO TO 1
      RUNTES=RUNTES+1.D0
      IF (RUNTES.LT.20.D0) WRITE(LUNOUT,602)N
 602   FORMAT(3H N=,I10,30H THE PROEKTILE IS A RESONANCE  )
      IF(IBAR(N).EQ.1) N=8
      IF(IBAR(N).EQ.-1)  N=9
   1  CONTINUE
      IMACH=IMACH+1
      IF (IMACH.GT.10) GO TO 8
      ECM =ECMMH
      AMN2=AMN**2
      AMT2=AMT**2
      ECMN=(ECM**2+AMN2-AMT2)/(2.D0*ECM)
      IF(ECMN.LE.AMN) ECMN=AMN
      PCMN=SQRT(ECMN**2-AMN2)
      GAM=(ELAB+AMT)/ECM
      BGAM=PLAB/ECM
      IF (IANTH.GE.0) ECM=2.1D0
C
C-----------------------------
C*** RANDOM CHOICE OF REACTION CHANNEL
C----------------------------
      IST=0
      CALL GRNDM(RNDM,1)
      VV = RNDM (1)
C
C-----------------------------
C***  PLACE REDUCED VERSION
C----------------------------
      IIEI=IEII(IRE)
      IDWK=IEII(IRE+1)-IIEI
      IIWK=IRII(IRE)
      IIKI=IKII(IRE)
C
C-----------------------------
C***  SHRINKAGE TO THE CONSIDERED ENERGY REGION FOR THE USE OF WEIGHTS
C----------------------------
      HECM=ECM
* The following cards assure that Ecm =< Umax + DUmax for this reaction
* where:
*       Umax  = max cms energy at which data are tabulated
*       DUmax = width of the last interval for the tabulated data
      HUMO=2.D0*UMO(IIEI+IDWK)-UMO(IIEI+IDWK-1)
      IF (HUMO.LT.ECM) ECM=HUMO
C
C-----------------------------
C*** INTERPOLATION PREPARATION
C----------------------------
* Cms energy of the upper limit of the considered energy interval
      ECMO=UMO(IE)
* Cms energy of the lower limit of the considered energy interval
      ECM1=UMO(IE-1)
* Width of the considered interval
      DECM=ECMO-ECM1
* Width from actual value to the upper limit (note that if Ecm > Ecmo
* it can be negative but its | | is always less than decm for the
* above condition on Humo
      DEC=ECMO-ECM
C
C-----------------------------
C*** RANDOM LOOP
C----------------------------
* Ik : index of the exit channel
      IK=0
      VFW=VV
      WKK=0.D+00
      WICOR=0.D+00
  111 CONTINUE
      IK=IK+1
* Save the weight accumulated up to now
      VFWO=WKK
* Get the index for the weight of ikth channel at ieth energy (upper
* limit of the interval in energy)
      IWK=IIWK+(IK-1)*IDWK+IE-IIEI
* Cumulative Weight of channels 1...ik at energy Ie
      WOK=WK(IWK)
* Difference for the cumulative weight of channels 1...ik at Ie and Ie-1
      WDK=WOK-WK(IWK-1)
* This card is not clear at all, it should zeroes the weight difference
* if we are in the first interval (that means take only the weights
* at the upper boundary of the interval)
      IF (PLAB.LT.PLABF(IIEI+2)) WDK=0.D+00
C
C-----------------------------
C*** TESTVARIABLE WICO/WICOR: IF CHANNEL IK HAS THE SAME WEIGHTS LIKE IK
C    GO TO NEXT CHANNEL, BECAUSE WKK((IK))-WKK((IK-1))=0, IK CAN NOT
C    CONTRIBUTE
C----------------------------
      WICO=WOK*1.23459876D0+WDK*1.735218469D0
      IF (WICO.EQ.WICOR) GO TO 111
* This card zeroes the weight difference if we are beyond the last
* interval
      IF (UMO(IIEI+IDWK).LT.HECM) WDK=0.D+00
* Save wico
      WICOR=WICO
C
C-----------------------------
C*** INTERPOLATION IN CHANNEL WEIGHTS
C----------------------------
* Set Eklim to a negative value to flag for Iefun it is a
* cms energy and not a lab momentum: this is the cms threshold
* for the exit channel ik
      EKLIM=-THRESH(IIKI+IK)
* Iefun returns the energy index of upper limit of the interval
* containing the threshold
      IELIM=IEFUN(EKLIM,IRE)
* Compute the difference between the upper limit and the
* threshold
      DELIM=UMO(IELIM)+EKLIM+ANGLSQ
* Dete is twice the difference between the actual energy value and
* the average value between Ecmo and the threshold
      DETE=(ECM-(ECMO-EKLIM)*.5D0)*2.D0
      IF ( DELIM*DELIM-DETE*DETE .GT. 0.D+00 ) THEN
         DECC=DELIM
      ELSE
         DECC=DECM
      END IF
      WKK=WOK-WDK*DEC/(DECC+ANGLSQ)
C
C-----------------------------
C*** RANDOM CHOICE
C----------------------------
C
      IF (VV.GT.WKK) GO TO 111
C
C***IK IS THE REACTION CHANNEL
C----------------------------
      INRK=IKII(IRE)+IK
      ECM=HECM
      I1001 =0
C
 1001 CONTINUE
      IT1=NRK(1,INRK)
      AM1=AMGA(IT1)
      IT2=NRK(2,INRK)
      AM2=AMGA(IT2)
      AMS=AM1+AM2
      I1001=I1001+1
      IF (I1001.GT.50) GO TO 1
C
      IF (IT2*AMS.GT.IT2*ECM) GO TO 1001
      IT11=IT1
      IT22=IT2
      IF (IANTH.GE.0) ECM=ELAB+AMT+0.00000001D0
      AM11=AM1
      AM22=AM2
      IF (IT2.GT.0) GO TO 401
C
C-----------------------------
C  INCLUSION OF DIRECT RESONANCES
C  RANDOM CHOICE OF DECAY CHANNELS OF THE DIRECT RESONANCE  IT1
C------------------------
      KZ1=K1(IT1)
      IST=IST+1
      IECO=0
*   Here was the mistake in the pseudo-masses treatment!!!!!
*     ECO=ECM
      ECO=ECMMH
      GAM=(ELAB+AMT)/ECO
      BGAM=PLAB/ECO
      CXS(1)=CX
      CYS(1)=CY
      CZS(1)=CZ
      GO TO 310
  401 CONTINUE
      CALL GRNDM(RNDM,1)
      IF ( RNDM(1) .LT. 0.5D0) GO TO 902
      IT1=IT22
      IT2=IT11
      AM1=AM22
      AM2=AM11
  902   CONTINUE
C
C-----------------------------
C   THE FIRST PARTICLE IS DEFINED TO BE THE FORWARD GOING ONE AT SMALL T
      IBN=IBAR(N)
      IB1=IBAR(IT1)
      IT11=IT1
      IT22=IT2
      AM11=AM1
      AM22=AM2
      IF(IB1.EQ.IBN) GO TO 901
      IT1=IT22
      IT2=IT11
      AM1=AM22
      AM2=AM11
 901  CONTINUE
C-----------------------------
C***IT1,IT2 ARE THE CREATED PARTICLES
C***MOMENTA AND DIRECTION COSINA IN THE CM - SYSTEM
C------------------------
      CALL TWOPAR(ECM1,ECM2,PCM1,PCM2,COD1,COD2,COF1,COF2,SIF1,SIF2,
     *IT1,IT2,ECM,ECMN,PCMN,N,AM1,AM2)
      IST=IST+1
      ITS(IST)=IT1
      AMM(IST)=AM1
C
C-----------------------------
C***TRANSFORMATION INTO LAB SYSTEM AND ROTATION
C----------------------------
      CALL TRAFO(GAM,BGAM,CX,CY,CZ,COD1,COF1,SIF1,PCM1,ECM1,PLS(IST),
     *CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IST=IST+1
      ITS(IST)=IT2
      AMM(IST)=AM2
      CALL TRAFO(GAM,BGAM,CX,CY,CZ,COD2,COF2,SIF2,PCM2,ECM2,PLS(IST),CXS
     *(IST),CYS(IST),CZS(IST),ELS(IST))
  200 CONTINUE
C
C-----------------------------
C***TEST   STABLE OR UNSTABLE
C----------------------------
      IF(ITS(IST).GT.NSTAB) GO TO 300
      IR=IR+1
C
C-----------------------------
C***IR IS THE NUMBER OF THE FINAL STABLE PARTICLE
C----------------------------
      IF (REDU.LT.0.D0) GO TO 1009
      ITR(IR)=ITS(IST)
       PLR(IR)=PLS(IST)
      CXR(IR)=CXS(IST)
      CYR(IR)=CYS(IST)
      CZR(IR)=CZS(IST)
      ELR(IR)=ELS(IST)
      IST=IST-1
      IF(IST.GE.1) GO TO 200
         GO TO 500
  300 CONTINUE
C
C  RANDOM CHOICE OF DECAY CHANNELS
C----------------------------
C
      IT=ITS(IST)
      ECO=AMM(IST)
      GAM=ELS(IST)/ECO
      BGAM=PLS(IST)/ECO
      IECO=0
      KZ1=K1(IT)
  310 CONTINUE
      IECO=IECO+1
      CALL GRNDM(RNDM,1)
      VV=RNDM(1)
      IIK=KZ1-1
  301 CONTINUE
      IIK=IIK+1
      IF (VV.GT.WT(IIK)) GO TO 301
C
C  IIK IS THE DECAY CHANNEL
C----------------------------
      IT1=NZK(IIK,1)
      I310=0
 1310 CONTINUE
      I310=I310+1
      AM1=AMGA(IT1)
      IT2=NZK(IIK,2)
      AM2=AMGA(IT2)
      IF (IT2-1.LT.0) GO TO 110
      IT3=NZK(IIK,3)
      AM3=AMGA(IT3)
      AMS=AM1+AM2+AM3
C
C  IF  IIK-KIN.LIM.GT.ACTUAL TOTAL CM-ENERGY, DO AGAIN RANDOM IIK-CHOICE
C----------------------------
      IF (IECO.LE.10) GO TO 1002
      IATMPT=IATMPT+1
* Note: we can go to 8 also for too many iterations
      IF (IATMPT.GT.3) GO TO 8
      GO TO 1000
 1002 CONTINUE
      IF (I310.GT.50) GO TO 310
      IF (AMS.GT.ECO) GO TO 1310
C
C  FOR THE DECAY CHANNEL
C  IT1,IT2, IT3 ARE THE PRODUCED PARTICLES FROM  IT
C----------------------------
      IF (REDU.LT.0.D0) GO TO 1009
      ITWTHC=0
      REDU=2.D0
      IF(IT3.EQ.0) GO TO 400
 4001 CONTINUE
      ITWTH=1
      CALL THREPD(ECO,ECM1,ECM2,ECM3,PCM1,PCM2,PCM3,COD1,COF1,SIF1,
     *COD2,COF2,SIF2,COD3,COF3,SIF3,AM1,AM2,AM3)
      GO TO 411
  400 CONTINUE
      CALL TWOPAD(ECO,ECM1,ECM2,PCM1,PCM2,COD1,COF1,SIF1,COD2,COF2,SIF2,
     *AM1,AM2)
      ITWTH=-1
      IT3=0
  411 CONTINUE
      ITWTHC=ITWTHC+1
      IF (REDU.GT.0.D0) GO TO 110
      REDU=2.D0
      IF (ITWTHC.GT.100) GO TO 1009
      IF (ITWTH) 400,400,4001
  110 CONTINUE
      ITS(IST  )=IT1
      IF (IT2-1.LT.0) GO TO 305
      ITS(IST+1)  =IT2
      ITS(IST+2)=IT3
      RX=CXS(IST)
      RY=CYS(IST)
      RZ=CZS(IST)
      AMM(IST)=AM1
      CALL TRAFO(GAM,BGAM,RX,RY,RZ,COD1,COF1,SIF1,PCM1,ECM1,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IST=IST+1
      AMM(IST)=AM2
      CALL TRAFO(GAM,BGAM,RX,RY,RZ,COD2,COF2,SIF2,PCM2,ECM2,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
      IF (IT3.LE.0) GO TO 305
      IST=IST+1
      AMM(IST)=AM3
      CALL TRAFO(GAM,BGAM,RX,RY,RZ,COD3,COF3,SIF3,PCM3,ECM3,
     *PLS(IST),CXS(IST),CYS(IST),CZS(IST),ELS(IST))
  305 CONTINUE
      GO TO 200
  500 CONTINUE
  631 CONTINUE
      RETURN
    8 CONTINUE
C
C----------------------------
C
C   ZERO CROSS SECTION CASE
C
C----------------------------
C
         IR=1
         ITR(1)=N
         CXR(1)=CX
         CYR(1)=CY
         CZR(1)=CZ
         ELR(1)=ELAB
         PLR(1)=PLAB
      RETURN
      END
