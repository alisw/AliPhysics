*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.43  by  S.Giani
*-- Author :
*$ CREATE HADRIV.FOR
*COPY HADRIV
*
*=== hadriv ===========================================================*
*
*
*----------------------------------------------------------------------*
*                                                                      *
*    Modified version of Hadrin created by Alfredo Ferrari, INFN-Milan *
*                                                                      *
*    Last change  on  20-jun-92  by  Alfredo Ferrari, INFN - MIlan     *
*                                                                      *
*    hadriv: this is a modified version of Hadrin, used by the Eventv  *
*            package for hadron-hadron interactions below 5 GeV        *
*                                                                      *
*----------------------------------------------------------------------*
*
      SUBROUTINE HADRIV ( N, PLAB, ELAB, CX, CY, CZ, ITTA )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/finlsp.inc"
#include "geant321/hadflg.inc"
#include "geant321/metlsp.inc"
#include "geant321/reac.inc"
#include "geant321/redver.inc"
#include "geant321/split.inc"
*
      COMMON / FKGAMR / REDU, AMO, AMM (15)
C
      COMMON / FKABLT / AM   (110), GA   (110), TAU  (110), ICH   (110),
     &                  IBAR (110), K1   (110), K2   (110)
      COMMON / FKCODS / COD1, COD2, COD3, COF1, COF2, COF3, SIF1, SIF2,
     &                  SIF3, ECM1, ECM2, ECM3, PCM1, PCM2, PCM3
      COMMON / FKRUN    / RUNTES, EFTES
*
      PARAMETER ( AMPROT = 0.93827231D+00 )
      DIMENSION WCHANN (40), WCUMCH (0:40), IKIK (40)
      DIMENSION ITPRF(110)
      REAL RNDM(1)
      LOGICAL LSWAP
*
      SAVE IKIK,ITPRF
      DATA NNN / 0 /
      DATA IKIK /  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
     &            16, 17, 18 ,19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
     &            29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 /
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
      INEWHD = N + 1000 * ITTA
      IF ( INEWHD .NE. IOLDHD ) THEN
         CALL CALUMV (N,ITTA)
      ELSE IF ( (AMPROT-AM(1)) .GT. 1.D-6 ) THEN
         CALL CALUMV (N,ITTA)
         INEWHD = - INEWHD
      END IF
      IOLDHD = INEWHD
 1009 CONTINUE
      IATMPT=0
      LOWP=LOWP+1
 1000 CONTINUE
      IMACH=0
      REDU=2.D0
      IW1=0
      IF (LOWP.GT.20) GO TO 8
      NNN=N
*  +-------------------------------------------------------------------*
*  |  The following condition is never verified
      IF (NNN.NE.N) THEN
         RUNTES=0.D0
         EFTES=0.D0
      END IF
*  |
*  +-------------------------------------------------------------------*
      IS=1
      IR=0
      IST=1
      NSTAB=25
*  +-------------------------------------------------------------------*
*  |  Select the reaction channel Ire: proton target
      IF ( ITTA .EQ. 1 ) THEN
         IRE = NURE(N,1)
*  |
*  +-------------------------------------------------------------------*
*  |  neutron target
      ELSE
         IRE = NURE(N,2)
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Elastic scattering index:
      IELSCT = MIN (N,ITTA) + 1000 * MAX (N,ITTA)
C
C-----------------------------
C*** IE,AMT,ECM,SI DETERMINATION
C----------------------------
      CALL FKSIGI(IRE,PLAB,N,IE,AMT,AMN,ECM,SI,ITTA)
*  +------------------------------------------------------------------*
*  |  Check if masses have been changed for annihilation treated with
*  |  pseudo masses
      IF ( AMPROT - AM(1) .GT. 1.D-6 ) THEN
         IANTH = 1
         SI = 1.D0
*  |
*  +------------------------------------------------------------------*
*  |
      ELSE
         IANTH = -1
      END IF
*  |
*  +------------------------------------------------------------------*
      ECMMH=ECM
C
C-----------------------------
C    ENERGY INDEX
C  IRE CHARACTERIZES THE REACTION
C  IE IS THE ENERGY INDEX
C----------------------------
      IF (SI.LT.1.D-10) GO TO 8
* The following condition should be always verified
      IF (N .LE.NSTAB) GO TO 1
      RUNTES=RUNTES+1.D0
      IF (RUNTES.LT.20.D0) WRITE(LUNOUT,602)N
 602   FORMAT(3H N=,I10,30H THE PROEKTILE IS A RESONANCE  )
      IF(IBAR(N).EQ.1) N=8
      IF(IBAR(N).EQ.-1)  N=9
*  **** Come here ( 1 continue ) every time we unsuccessfully selected
*       for more than 50 times the mass of the produced resonaces ****
   1  CONTINUE
      IMACH=IMACH+1
      IF (IMACH.GT.10) GO TO 8
      ECM =ECMMH
      AMN2=AMN**2
      AMT2=AMT**2
*  CMS energy of the projectile
      ECMN=(ECM**2+AMN2-AMT2)/(2.D0*ECM)
*  It should never happen
      IF(ECMN.LE.AMN) ECMN=AMN
*  CMS momentum of the projectile (and of the target)
      PCMN=SQRT(ECMN**2-AMN2)
      GAM=(ELAB+AMT)/ECM
      BGAM=PLAB/ECM
      IF (IANTH.GE.0) ECM=2.1D0
*  From this point starts the random choiche of the reaction channel:
*  it was extensively modified by A. Ferrari
      IST=0
* Initial energy index for the reaction IRE (index 0)
      IIEI=IEII(IRE)
* Number of energy intervals for the reaction IRE
      IDWK=IEII(IRE+1)-IIEI
* Initial index for the exit channel weights
      IIWK=IRII(IRE)
* Initial index (for 0) of the exit channels
      IIKI=IKII(IRE)
* Number of exit channels of reaction ire
      IKE =IKII(IRE+1)-IIKI
* *** Shrinkage to the considered energy region for the use of weights
      HECM=ECM
* The following cards assure that Ecm =< Umax + DUmax for this reaction
* where:
*       Umax  = max cms energy at which data are tabulated
*       DUmax = width of the last interval for the tabulated data
      HUMO=2.D0*UMO(IIEI+IDWK)-UMO(IIEI+IDWK-1)
      IF (HUMO.LT.ECM) ECM=HUMO
* *** Interpolation preparation
* Cms energy of the upper limit of the considered energy interval
      ECMO=UMO(IE)
* Cms energy of the lower limit of the considered energy interval
      ECM1=UMO(IE-1)
* Width of the considered interval
      DECM=ECMO-ECM1
* Width from actual value to the lower limit (note that if Ecm > Ecmo
* it can be larger than Decm but it is always less than 2xdecm for the
* above condition on Humo
      DEC0=ECM -ECM1
* Set to 1 the default total weight
      WACCUM = 1.D+00
      WCUMCH (0) = 0.D+00
      WCUMIE = 0.D+00
      WCUMI0 = 0.D+00
      WCSUM0 = 0.D+00
      CALL GRNDM(RNDM,1)
      RNDMIK = RNDM (1)
      IOUT1  = NRK (1,IIKI+1)
      IOUT2  = NRK (2,IIKI+1)
      IELCHK = MIN ( IOUT1, IOUT2 ) + 1000 * MAX ( IOUT1, IOUT2 )
*  +-------------------------------------------------------------------*
*  |  Look for "inverse" reactions for which the elastic channel is not
*  |  the first: for these reactions we must exchange the weight of
*  |  the elastic channel with the charge exchange one at minimum, to
*  |  fulfill the detailed balance theorem
      IF ( IELCHK .NE. IELSCT ) THEN
         LSWAP  = .TRUE.
         IELCHA = IKCHXG (IRE)
         ICXCHA = 1
         IKIK (1) = IELCHA
         IKIK (IELCHA) = 1
*  |
*  +-------------------------------------------------------------------*
*  |  Loop to find the elastic channel
      ELSE
         LSWAP  = .FALSE.
         IELCHA = 1
         ICXCHA = IKCHXG (IRE)
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Loop on the exit channels:
*  |                              Jkk = index for weights
*  |                              Ik  = index for reaction
*  |  usually they are the same, but for "inverse" reactions where the
*  |  elastic and the charge exchange channels are exchanged, for these
*  |  two channels they are crossed
      DO 2000 JKK = 1, IKE
*  |  Ik : index of the exit channel under consideration
         IK = IKIK (JKK)
*  |  Get the index for the weight of ikth channel at ieth energy (upper
*  |  limit of the interval in energy)
         IWK = IIWK + (JKK-1) * IDWK + IE - IIEI
*  |  Cumulative Weight of channels 1...ik at energy Ie
         WIEK   = WK (IWK)
*  |  +----------------------------------------------------------------*
*  |  |  Check if we are in the first interval: is so all the weights
*  |  |  are set to zero for Ie-1, for all channels, so set them to the
*  |  |  same value as for Ie to get the proper normalization to 1,
*  |  |  then possible thresholds are accounted for after
         IF ( IE .LE. IIEI+2 ) THEN
            WIEM1K = WIEK
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Cumulative Weight of channels 1...ik at energy Ie-1
         ELSE
            WIEM1K = WK (IWK-1)
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  Compute the weight at Ie for this channel
         WIEK   = WIEK   - WCUMIE
*  |  Compute the weight at Ie-1 for this channel
         WIEM1K = WIEM1K - WCUMI0
*  |  +----------------------------------------------------------------*
*  |  |  This channel is not open at energies Ie and Ie-1
         IF ( WIEM1K + WIEK .LE. ANGLGB ) THEN
            WCHANN (JKK) = 0.D+00
            WCUMCH (JKK) = WCUMCH (JKK-1)
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
*  |  |  Set Eklim to a negative value to flag for Iefun it is a
*  |  |  cms energy and not a lab momentum: this is the cms threshold
*  |  |  for the exit channel ik
            EKLIM = - THRESH (IIKI+IK)
*  |  |  Iefun returns the energy index of upper limit of the interval
*  |  |  containing the threshold
            IELIM = IEFUN (EKLIM,IRE)
            EKLIM = - EKLIM
            WCHAN0 = WIEM1K + ( WIEK - WIEM1K ) * DEC0 / DECM
            WCSUM0 = WCSUM0 + WCHAN0
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check if we are below threshold
            IF ( ECM .LE. EKLIM ) THEN
               WCHANN (JKK) = 0.D+00
               WCUMCH (JKK) = WCUMCH (JKK-1)
               WCUMIE = WCUMIE + WIEK
               WCUMI0 = WCUMI0 + WIEM1K
               RNDRED = WCHAN0
               WACCUM = WACCUM - RNDRED
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  We are above threshold
            ELSE
               WCUMIE = WCUMIE + WIEK
               WCUMI0 = WCUMI0 + WIEM1K
               ECMD   = MAX ( EKLIM, ECM1 )
               DEC    = ECM  - ECMD
               DECC   = ECMO - ECMD
               WCHANN (JKK) = WIEM1K + ( WIEK - WIEM1K ) * DEC / DECC
*  |  |  |  If we are beyond the last tabulated point and the xsec for
*  |  |  |  this channel is going down it can happen that Wchann < 0
*  |  |  |  set it to 0 and correct according to the usual formalism
               WCHANN (JKK) = MAX ( WCHANN (JKK), ZERZER )
               RNDRED = WCHAN0 - WCHANN (JKK)
               WACCUM = WACCUM - RNDRED
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Ielflg check: first of all check if this channel
*  |  |  |  |  is the elastic one and if so if reduction must be applied
*  |  |  |  |  to
               IF ( IK .EQ. IELCHA .AND. IELFLG .NE. 0 ) THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  Elastic scattering reduced due to collision in the
*  |  |  |  |  |  nucleus and not with a free proton/neutron
                  IF ( IELFLG .LT. 0 ) THEN
                     IF ( IBAR (N) .NE. 0 ) THEN
                        REDUC = PLAB / PPAMXB
                        REDUC = PAUMXB * REDUC / ( 1.D+00 + REDUC**2 )
                     ELSE
                        REDUC = PLAB / PPAMXM
                        REDUC = PAUMXM * REDUC / ( 1.D+00 + REDUC**2 )
                     END IF
                     RNDRED = WCHANN (JKK) * ( 1.D+00 - REDUC )
                     WCHANN (JKK) = WCHANN (JKK) - RNDRED
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  Elastic scattering forbidden
                  ELSE
                     RNDRED = WCHANN (JKK)
                     WCHANN (JKK) = 0.D+00
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  WACCUM = WACCUM - RNDRED
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Icxflg check: first of all check if this channel
*  |  |  |  |  is the charge exchange one and if so if reduction must
*  |  |  |  |  be applied to
               ELSE IF ( IK .EQ. ICXCHA .AND. ICXFLG .NE. 0 )THEN
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |  Charge exchange reduced due to collision in the
*  |  |  |  |  |  nucleus and not with a free proton/neutron
                  IF ( ICXFLG .LT. 0 ) THEN
                     IF ( IBAR (N) .NE. 0 ) THEN
                        REDUC = PLAB / PPAMXB
                        REDUC = PAUMXB * REDUC / ( 1.D+00 + REDUC**2 )
                     ELSE
                        REDUC = PLAB / PPAMXM
                        REDUC = PAUMXM * REDUC / ( 1.D+00 + REDUC**2 )
                     END IF
                     RNDRED = WCHANN (JKK) * ( 1.D+00 - REDUC )
                     WCHANN (JKK) = WCHANN (JKK) - RNDRED
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Charge exchange forbidden
                  ELSE
                     RNDRED = WCHANN (JKK)
                     WCHANN (JKK) = 0.D+00
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  WACCUM = WACCUM - RNDRED
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               WCUMCH (JKK) = WCUMCH (JKK-1) + WCHANN (JKK)
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            RNDMCH = RNDMIK * WACCUM
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check if a possible decrease/increase of this channel
*  |  |  |  opened one of the already examinated channels
            IF ( RNDMCH .LT. WCUMCH (JKK-1) ) THEN
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Loop on the previous channels
               DO 1900 JPK = 1, JKK - 1
                  IF ( RNDMCH .LT. WCUMCH (JPK) ) THEN
                     IK = IKIK (JPK)
                     GO TO 2100
                  END IF
 1900          CONTINUE
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check if this one is the right channel
            ELSE IF ( RNDMCH .LT. WCUMCH (JKK) ) THEN
               GO TO 2100
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 2000 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      IK = IELCHA
      WRITE (LUNERR,*)' **** Hadriv: elastic channel selected when',
     &                ' prohibited ! ****',N,ELAB,PLAB
*  Finally we selected channel ik
 2100 CONTINUE
* *** Ik is the reaction channel ***
*  +-------------------------------------------------------------------*
*  |  Set to the default values the array Ikik
      IF ( LSWAP ) THEN
         IKIK (1) = 1
         IKIK (IELCHA) = IELCHA
      END IF
*  |
*  +-------------------------------------------------------------------*
      INRK=IKII(IRE)+IK
*  First resonance to be created
      IT1=NRK(1,INRK)
*  Second resonance to be created
      IT2=NRK(2,INRK)
      ECM=HECM
      I1001 =0
*  +-------------------------------------------------------------------*
*  |  Rejection loop for the choiche of the resonance masses
 1001 CONTINUE
         IF (I1001.GT.50) GO TO 1
*  |  Selection of the resonance mass according to its width
         AM1=AMGA(IT1)
*  |  Selection of the resonance mass according to its width
         AM2=AMGA(IT2)
         AMS=AM1+AM2
         I1001=I1001+1
      IF ( IT2*AMS .GT. IT2*ECM ) GO TO 1001
*  |--<--<--<--<--<--<  Loop back if m1+m2 > Ecm
*  +-------------------------------------------------------------------*
      IT11=IT1
      IT22=IT2
      IF (IANTH.GE.0) ECM=ELAB+AMT+0.00000001D0
      AM11=AM1
      AM22=AM2
*  +-------------------------------------------------------------------*
*  |  Direct (single) resonances
      IF (IT2.LE.0) THEN
*  | Random choice of decay channels of the direct resonance  it1
         KZ1=K1(IT1)
         IST=IST+1
         IECO=0
*  |   Here was the mistake in the pseudo-masses treatment!!!!!
*        ECO=ECM
         ECO=ECMMH
         GAM=(ELAB+AMT)/ECO
         BGAM=PLAB/ECO
         CXS(1)=CX
         CYS(1)=CY
         CZS(1)=CZ
         GO TO 310
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      CALL GRNDM(RNDM,1)
      IF ( RNDM(1) .LT. 0.5D+00 ) THEN
         IT1=IT22
         IT2=IT11
         AM1=AM22
         AM2=AM11
      END IF
*  |
*  +-------------------------------------------------------------------*
C
C-----------------------------
C   THE FIRST PARTICLE IS DEFINED TO BE THE FORWARD GOING ONE AT SMALL T
      IBN=IBAR(N)
      IB1=IBAR(IT1)
      IT11=IT1
      IT22=IT2
      AM11=AM1
      AM22=AM2
*  +-------------------------------------------------------------------*
*  |
      IF(IB1.NE.IBN) THEN
         IT1=IT22
         IT2=IT11
         AM1=AM22
         AM2=AM11
      END IF
*  |
*  +-------------------------------------------------------------------*
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
         IF (IECO.GT.10) THEN
            IATMPT=IATMPT+1
* Note: we can go to 8 also for too many iterations
            IF (IATMPT.GT.3) GO TO 8
            GO TO 1000
         END IF
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
      ITS(IST)=IT1
      IF (IT2.LE.0) GO TO 305
      ITS(IST+1)=IT2
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
