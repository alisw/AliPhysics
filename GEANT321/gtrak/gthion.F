#include "geant321/pilot.h"
*CMZ :  3.21/02 03/07/94  17.58.49  by  S.Giani
*-- Author :
      SUBROUTINE GTHION
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Heavy  ion  type track. Computes step size and propagates    *
C.    *    particle through step.                                      *
C.    *                                                                *
C.    *   The ionisation energy loss is calculated here (mean +        *
C.    *       fluctuations)                                            *
C.    *   The fluctuations are the same for ILOSS=1,2,3 and            *
C.    *       there is no fluctuation for ILOSS=4.                     *
C.    *                                                                *
C.    *   ==>Called by : GTRACK                                        *
C.    *       Authors    R.Brun, F.Bruyant, M.Maire, L.Urban ***       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
 
#if !defined(CERNLIB_OLD)
#include "geant321/gcvolu.inc"
#include "geant321/gcvdma.inc"
#endif
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
      DOUBLE PRECISION GKR,DEMEAN,STOPP,STOPMX,STOPRG,STOPC,EKIPR
      DOUBLE PRECISION ONE,XCOEF1,XCOEF2,XCOEF3,YCOEF1,YCOEF2,YCOEF3
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
      PARAMETER (THRESH=0.7,ONE=1)
      PARAMETER (TWOTHR=2*ONE/3,AMU=0.9314943)
      PARAMETER (DME=7.84572E-8,CNORM=2.5)
      REAL VNEXT(6)
      DIMENSION RNDM(2)
      SAVE RMASS,CUTPRO,IKCUT,STOPC,FACFLU,CHAR23
C.
C.    ------------------------------------------------------------------
*
* *** Particle below energy threshold ? short circuit
*
      IF (GEKIN.LE.CUTHAD) GO TO 100
*
* *** Update local pointers if medium has changed
*
      IF (IUPD.EQ.0) THEN
         IUPD  = 1
         JLOSS = LQ(JMA-3)
         JRANG = LQ(JMA-16) + NEK1
         JCOEF = LQ(JMA-18) + 3*NEK1
         RMASS  = PMASS/AMASS
         CUTPRO = MAX(CUTHAD*RMASS,ELOW(1))
         IKCUT = GEKA*LOG10(CUTPRO) + GEKB
         GKR   = (CUTPRO - ELOW(IKCUT))/(ELOW(IKCUT+1) - ELOW(IKCUT))
         STOPC = (1.-GKR)*Q(JRANG+IKCUT) + GKR*Q(JRANG+IKCUT+1)
         FACFLU = DME*(Z*DENS/A)
         CHAR23 = ONE/CHARGE**TWOTHR
         IF(IMCKOV.EQ.1) THEN
            JTCKOV = LQ(JTM-3)
            JABSCO = LQ(JTCKOV-1)
            JEFFIC = LQ(JTCKOV-2)
            JINDEX = LQ(JTCKOV-3)
            JCURIN = LQ(JTCKOV-4)
            NPCKOV = Q(JTCKOV+1)
         ENDIF
      ENDIF
*
* *** Compute energy dependent parameters
*
      GAMASS=GETOT+AMASS
      BET2=GEKIN*GAMASS/(GETOT*GETOT)
      BET=SQRT(BET2)
      W1=1.034-0.1777*EXP(-0.08114*CHARGE)
      W2=BET*CHAR23
      W3=121.4139*W2+0.0378*SIN(190.7165*W2)
      CHARG1=CHARGE*(1.-W1*EXP(-W3))
*
*              the effective charge  CHARG1
*            can be negative only for very low energy and
*     for CHARGE > 20 ( very low energy : T/A < 20 keV/nucleon)
*              in this case short circuit
*
      IF(CHARG1.LT.0.) GOTO 100
      CHARG2=CHARG1**2
*
      OMCMOL=Q(JPROB+21)*CHARG2
      CHCMOL=Q(JPROB+25)*ABS(CHARG1)
      IF(FIELDM.NE.0.) THEN
         CFLD=3333.*DEGRAD*TMAXFD/ABS(FIELDM*CHARG1)
      ELSE
         CFLD=BIG
      ENDIF
*
* *** Compute current step size
*
      STEP   = STEMAX
      IPROC  = 103
      GEKRT1 = 1. -GEKRAT
*
*  **   Step limitation due to hadron interaction ?
*
      IF (IHADR.GT.0) THEN
#if !defined(CERNLIB_USRJMP)
         CALL GUPHAD
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT0(JUPHAD)
#endif
         IF (SHADR.LT.STEP) THEN
            IF (SHADR.LE.0.) SHADR = PREC
            STEP  = SHADR
            IPROC = 12
         ENDIF
      ENDIF
*
*  ** Step limitation due to delta-ray production ?
*       (Cannot be tabulated easily because dependent on AMASS)
*
      IF (IDRAY.GT.0) THEN
         STEPDR = BIG
         IF (GEKIN.GT.DCUTM) THEN
            TMAX   = EMASS*GEKIN*GAMASS/(0.5*AMASS*AMASS+EMASS*GETOT)
            IF (TMAX.GT.DCUTM) THEN
               Y    = DCUTM/TMAX
               SIG  = (1.-Y+BET2*Y*LOG(Y))/DCUTM
*              extra term for spin 1/2
               IF (AMASS.GT.0.9) SIG=SIG+0.5*(TMAX-DCUTM)/(GETOT*GETOT)
               SIG = SIG*Q(JPROB+17)*CHARG2*EMASS/BET2
*
               IF (SIG.GT.0.) THEN
                  STEPDR = 1./SIG
                  SDRAY  = STEPDR*ZINTDR
                  IF (SDRAY.LE.STEP) THEN
                     STEP  = SDRAY
                     IPROC = 10
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
*
      IF (STEP.LE.0.) THEN
         STEP  = 0.
         GO TO 110
      ENDIF
*
*  **   Step limitation due to energy loss (stopping range) ?
*
      IF (ILOSL.GT.0) THEN
         IF(GEKRAT.LT.THRESH) THEN
            I1 = MAX(IEKBIN-1,1)
         ELSE
            I1 = MIN(IEKBIN,NEKBIN-1)
         ENDIF
         I1 = 3*(I1-1)+1
         XCOEF1 = Q(JCOEF+I1)
         XCOEF2 = Q(JCOEF+I1+1)
         XCOEF3 = Q(JCOEF+I1+2)
         IF(XCOEF1.NE.0) THEN
            STOPP = -XCOEF2+SIGN(ONE,XCOEF1)* SQRT(XCOEF2
     +      **2 -(XCOEF3-GEKIN*RMASS/XCOEF1))
         ELSE
            STOPP = - (XCOEF3-GEKIN*RMASS)/XCOEF2
         ENDIF
         STOPMX = (STOPP - STOPC)/(RMASS*CHARG2)
         IF (STOPMX.LT.MIN(STEP,STMIN)) THEN
            STEP = STOPMX
            IPROC = 0
            IF(STEP.LE.0.)THEN
               GO TO 100
            ENDIF
            GO TO 10
         ENDIF
         EKF = (1. - DEEMAX)*GEKIN*RMASS
         IF (EKF.LT.ELOW(1)) THEN
            EKF = ELOW(1)
         ELSEIF (EKF.GE.ELOW(NEK1)) THEN
            EKF = ELOW(NEK1)*0.99
         ENDIF
         IKF=GEKA*LOG10(EKF)+GEKB
         GKR=(EKF-ELOW(IKF))/(ELOW(IKF+1)-ELOW(IKF))
         IF(GKR.LT.THRESH) THEN
            IK1 = MAX(IKF-1,1)
         ELSE
            IK1 = MIN(IKF,NEKBIN-1)
         ENDIF
         IK1 = 3*(IK1-1)+1
         YCOEF1=Q(JCOEF+IK1)
         YCOEF2=Q(JCOEF+IK1+1)
         YCOEF3=Q(JCOEF+IK1+2)
         IF(YCOEF1.NE.0.) THEN
            SLOSP = -YCOEF2+SIGN(ONE,YCOEF1)*SQRT(YCOEF2**2- (YCOEF3-
     +      EKF/YCOEF1))
         ELSE
            SLOSP = - (YCOEF3-EKF)/YCOEF2
         ENDIF
         SLOSP = STOPP - SLOSP
         SLOSS = MAX(STMIN, SLOSP/(RMASS*CHARG2) )
         IF (SLOSS.LT.STEP) THEN
            STEP = SLOSS
            IPROC = 0
         ENDIF
      ENDIF
*
*  **   Step limitation due to energy loss in magnetic field ?
*
      IF (IFIELD.NE.0) THEN
         SFIELD = CFLD*VECT(7)
         SFIELD=MAX(SFIELD, STMIN)
         IF (SFIELD.LT.STEP) THEN
            STEP  = SFIELD
            IPROC = 0
         ENDIF
      ENDIF
*
*  **   Step limitation due to multiple scattering ?
*
      IF (IMULL.GT.0) THEN
         SMULS=MIN(2232.*RADL*((VECT(7)**2)/(GETOT*CHARG1))**2,10.*RADL)
         SMULS  = MAX(STMIN, SMULS )
         IF (SMULS.LT.STEP) THEN
            STEP  = SMULS
            IPROC = 0
         ENDIF
      ENDIF
*
   10 CONTINUE
*
*  **   Step limitation due to Cerenkov production ?
*
      IF (IMCKOV.GT.0) THEN
         CALL GNCKOV
         STCKOV = MXPHOT/MAX(3.*DNDL,1E-10)
         SMULS  = MAX(STMIN, STCKOV)
         IF (SMULS.LT.STEP) THEN
            STEP  = STCKOV
            IPROC = 0
         ENDIF
      ENDIF
*
*  **   Step limitation due to geometry ?
*
      IF (STEP.GE.0.95*SAFETY) THEN
         CALL GTNEXT
         IF (IGNEXT.NE.0) THEN
            STEP   = SNEXT + PREC
            IPROC = 0
         ENDIF
*
*        Update SAFETY in stack companions, if any
         IF (IQ(JSTAK+3).NE.0) THEN
            DO 20 IST = IQ(JSTAK+3),IQ(JSTAK+1)
               JST    = JSTAK + 3 + (IST-1)*NWSTAK
               Q(JST+11) = SAFETY
   20       CONTINUE
            IQ(JSTAK+3) = 0
         ENDIF
      ELSE
         IQ(JSTAK+3) = 0
      ENDIF
*
* *** Linear transport when no field or very short step
*
      IF (IFIELD.EQ.0.OR.STEP.LE.PREC) THEN
*
         IF (IGNEXT.NE.0) THEN
            DO 30 I = 1,3
               VECTMP  = VECT(I) +STEP*VECT(I+3)
               IF(VECTMP.EQ.VECT(I)) THEN
*
* *** Correct for machine precision
*
                  IF(VECT(I+3).NE.0.) THEN
                     VECTMP =
     +               VECT(I)+ABS(VECT(I))*SIGN(1.,VECT(I+3))*EPSMAC
                     IF(NMEC.GT.0) THEN
                        IF(LMEC(NMEC).EQ.104) NMEC=NMEC-1
                     ENDIF
                     NMEC=NMEC+1
                     LMEC(NMEC)=104
#if defined(CERNLIB_DEBUG)
                     WRITE(CHMAIL, 10000)
                     CALL GMAIL(0,0)
                     WRITE(CHMAIL, 10100) GEKIN, NUMED, STEP, SNEXT
                     CALL GMAIL(0,0)
10000 FORMAT(' Boundary correction in GTHION: ',
     +       '    GEKIN      NUMED       STEP      SNEXT')
10100 FORMAT(31X,E10.3,1X,I10,1X,E10.3,1X,E10.3,1X)
#endif
                  ENDIF
               ENDIF
               VECT(I) = VECTMP
   30       CONTINUE
            INWVOL = 2
            NMEC = NMEC +1
            LMEC(NMEC) = 1
         ELSE
            DO 40 I = 1,3
               VECT(I)  = VECT(I) +STEP*VECT(I+3)
   40       CONTINUE
         ENDIF
      ELSE
*
* ***   otherwise, swim particle in magnetic field
*
#if !defined(CERNLIB_OLD)
         if(mycoun.gt.1.and.nfmany.gt.0.and.step.ge.safety)then
           nlevel=manyle(nfmany)
           do 99 i=1,nlevel
             names(i)=manyna(nfmany,i)
             number(i)=manynu(nfmany,i)
 99        continue
           call glvolu(nlevel,names,number,ier)
           if(ier.ne.0)print *,'Fatal error in GLVOLU'
           ingoto=0
         endif
#endif
         NMEC = NMEC +1
         LMEC(NMEC) = 4
*
#if !defined(CERNLIB_USRJMP)
   50    CALL GUSWIM (CHARG1, STEP, VECT, VOUT)
#endif
#if defined(CERNLIB_USRJMP)
   50    CALL JUMPT4(JUSWIM, CHARG1, STEP, VECT, VOUT)
#endif
*
*  ** When near to boundary, take proper action (cut-step,crossing...)
*
         IF(STEP.GE.SAFETY)THEN
            INEAR = 0
            IF (IGNEXT.NE.0) THEN
               DO 60 I = 1,3
                  VNEXT(I+3) = VECT(I+3)
                  VNEXT(I) = VECT(I) +SNEXT*VECT(I+3)
   60          CONTINUE
               DO 70 I = 1,3
                  IF (ABS(VOUT(I)-VNEXT(I)).GT.EPSIL) GO TO 80
   70          CONTINUE
               INEAR = 1
            ENDIF
*
   80       CALL GINVOL (VOUT, ISAME)
            IF (ISAME.EQ.0)THEN
               IF ((INEAR.NE.0).OR.(STEP.LT.EPSIL)) THEN
                  INWVOL = 2
                  NMEC = NMEC +1
                  LMEC(NMEC) = 1
               ELSE
*              Cut step
                  STEP = 0.5*STEP
                  IF (LMEC(NMEC).NE.24) THEN
                     NMEC = NMEC +1
                     LMEC(NMEC) = 24
                  ENDIF
                  GO TO 50
               ENDIF
            ENDIF
         ENDIF
*
         DO 90 I = 1,6
            VECT(I) = VOUT(I)
   90    CONTINUE
*
      ENDIF
*
* *** Correct the step due to multiple scattering
      IF (IMULL.NE.0) THEN
         STMULS = STEP
         CORR=0.0001*CHARG2*(STEP/RADL)*(GETOT/(VECT(7)*VECT(7)))**2
         IF (CORR.GT.0.25) CORR = 0.25
         STEP  = (1.+CORR)*STEP
      ENDIF
*
      SLENG = SLENG + STEP
*
* *** Generate Cherenkov photons if required
*
      IF(IMCKOV.EQ.1) THEN
         CALL GGCKOV
         NMEC=NMEC+1
         LMEC(NMEC)=105
      ENDIF
*
* *** apply energy loss : find the kinetic energy corresponding
*      to the new stopping range = stopmx - step
*
      IF (ILOSL.NE.0) THEN
         NMEC = NMEC +1
         LMEC(NMEC) = 3
         STOPRG = STOPP - STEP*RMASS*CHARG2
         IF (STOPRG.LE.STOPC) THEN
            STEP = STOPMX
            GO TO 100
         ENDIF
         IF(XCOEF1.NE.0.) THEN
            EKIPR = XCOEF1*(XCOEF3+STOPRG*(2.*XCOEF2+STOPRG))
         ELSE
            EKIPR = XCOEF2*STOPRG+XCOEF3
         ENDIF
         DEMEAN=GEKIN - EKIPR/RMASS
         IF(DEMEAN.LE.5.*GEKIN*EPSMAC) THEN
            DEMEAN=(GEKRT1*Q(JLOSS+IEKBIN)+GEKRAT*Q(JLOSS+IEKBIN+1))
     +             *STEP*CHARG2
         ENDIF
*
*        fluctuations : differ from that of 'ordinary' hadrons
*
         IF (ILOSS.EQ.4.OR.IEKBIN.LE.IKCUT+1) THEN
            DESTEP = DEMEAN
         ELSE
*
*     Charge exchange fluctuations + Gaussian 'Landau' fluctuations
*           (it is the same for ILOSS=1,2,3 !)
*
            SIGMA2=CNORM*CHARG1*(1.-CHARG1/CHARGE)
            SIGMA2=MAX(SIGMA2,0.)
            TA = RMASS*GEKIN
            TAM=TA/AMU
            SIGMA2=SIGMA2+2.+TAM*(2.+TAM)
*
            SIGMA2=FACFLU*CHARG2*STEP*SIGMA2
            IF(SIGMA2.GT.0.0) THEN
                SIGMA=SQRT(SIGMA2)
            ELSE
                SIGMA= 0.0
            END IF
*
*     Check if we are in 'Gaussian' regime ...
*
            CAPPA=(1.+TAM)/(TAM*(2.+TAM)*EMASS)
            CAPPA=0.5*CAPPA**2*FACFLU*CHARG2*STEP
*
*     ... if not , correct SIGMA !
 
            IF( (CAPPA.LT.10.) .AND. (CAPPA.GT.0.0) ) THEN
               SIGMA=SIGMA/(0.97+0.03*SQRT(10./CAPPA))
            ENDIF
*
            CALL GRNDM(RNDM,2)
            DEFLUC=SIGMA*SIN(TWOPI*RNDM(1))*SQRT(-2.*LOG(RNDM(2)))
            DESTEP=DEMEAN+DEFLUC
         ENDIF
*
*     protection against negative destep
*
         IF(DESTEP.LT.0.) DESTEP=DEMEAN
*                          IF (DESTEP.LT.0.) DESTEP = 0.
         GEKINT = GEKIN -DESTEP
         IF (GEKINT.LE.(1.01*CUTHAD)) GO TO 100
         DESTEL = DESTEP
         GEKIN  = GEKINT
         GETOT  = GEKIN +AMASS
         VECT(7)= SQRT((GETOT+AMASS)*GEKIN)
         CALL GEKBIN
      ENDIF
*
* *** Apply multiple scattering.
*
      IF (IMULL.NE.0) THEN
         NMEC = NMEC +1
         LMEC(NMEC) = 2
*      check charge dependence ...........!!!!!!!  (later..)
         CALL GMULTS
      ENDIF
*
* *** Update time of flight
*
      SUMLIF = SUMLIF -STEP*AMASS/VECT(7)
      TOFG   = TOFG +STEP*GETOT/(VECT(7)*CLIGHT)
      IF (TOFG.GE.TOFMAX) THEN
         ISTOP = 4
         NMEC  = NMEC +1
         LMEC(NMEC) = 22
         GO TO 999
      ENDIF
*
* *** Update interaction probabilities
*
      IF (IHADR.GT.0) ZINTHA = ZINTHA*(1.-STEP/SHADR)
      IF (IDRAY.GT.0) ZINTDR = ZINTDR -STEP/STEPDR
*
      GO TO 110
*
*  **   Special treatment for overstopped tracks
*
  100 DESTEP = GEKIN
      DESTEL = DESTEP
      GEKIN  = 0.
      GETOT  = AMASS
      VECT(7)= 0.
      INWVOL = 0
      ISTOP  = 2
      NMEC = NMEC + 1
      LMEC(NMEC) = 30
      IF (IHADR.EQ.0) GO TO 999
      IPROC = 12
*
* *** apply slected process if any
*
  110 IF (IPROC.EQ.0) GO TO 999
      NMEC = NMEC +1
      LMEC(NMEC) = IPROC
*
*  **   Hadron interaction ?
*
      IF (IPROC.EQ.12) THEN
#if !defined(CERNLIB_USRJMP)
         CALL GUHADR
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT0(JUHADR)
#endif
*   *   Check time cut-off for decays at rest
         IF (LMEC(NMEC).EQ.5) THEN
            TOFG   = TOFG +SUMLIF/CLIGHT
            SUMLIF = 0.
            IF (TOFG.GE.TOFMAX) THEN
               NGKINE = 0
               ISTOP  = 4
               LMEC(NMEC) = 22
            ENDIF
         ENDIF
*
*  **   Delta-ray ?
*
      ELSE IF (IPROC.EQ.10) THEN
         CALL GDRAY
      ENDIF
  999 END
