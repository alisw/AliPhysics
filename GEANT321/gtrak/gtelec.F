#include "geant321/pilot.h"
*CMZ :  3.21/04 22/02/95  16.08.31  by  S.Giani
*-- Author :
      SUBROUTINE GTELEC
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Electron type track. Computes step size and propagates       *
C.    *    particle through step.                                      *
C.    *                                                                *
C.    *   ==> Called by : GTRACK                                       *
C.    *       Authors    R.Brun, F.Bruyant, M.Maire L.Urban ********   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcking.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
 
#if !defined(CERNLIB_OLD)
#include "geant321/gcvolu.inc"
#include "geant321/gcvdma.inc"
#endif
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
      DOUBLE PRECISION DEMEAN,STOPRG,STOPMX,STOPC
      DOUBLE PRECISION ONE,XCOEF1,XCOEF2,XCOEF3,ZERO
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
      PARAMETER (ONE=1,ZERO=0.)
      REAL VNEXT(6)
      SAVE IKCUT,STOPC
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DIMENSION RNDM(3)
      PARAMETER ( TLIM = 0.0002)
      IABAN = NINT(DPHYS1)
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C.
C.    ------------------------------------------------------------------
*
* *** Particle below energy threshold ? short circuit
*
      IF (GEKIN.LE.CUTELE) GO TO 100
*
* *** Update local pointers if medium or particle code has changed
      IF (IUPD.EQ.0) THEN
         IUPD  = 1
         JMULOF= LQ(JTM-1)
         IF (CHARGE.LT.0.) THEN
            JBREM = LQ(JMA-9)
            JLOSS = LQ(JMA-1)
            JDRAY = LQ(JMA-11)
            JRANG = LQ(JMA-15)
            JCOEF = LQ(JMA-17)
         ELSE
            JBREM = LQ(JMA-9)  +NEK1
            JLOSS = LQ(JMA-1)  +NEK1
            JDRAY = LQ(JMA-11) +NEK1
            JRANG = LQ(JMA-15) +NEK1
            JCOEF = LQ(JMA-17) +3*NEK1
            JANNI = LQ(JMA-7)
         ENDIF
         IF(IMCKOV.EQ.1) THEN
            JTCKOV = LQ(JTM-3)
            JABSCO = LQ(JTCKOV-1)
            JEFFIC = LQ(JTCKOV-2)
            JINDEX = LQ(JTCKOV-3)
            JCURIN = LQ(JTCKOV-4)
            NPCKOV = Q(JTCKOV+1)
         ENDIF
         OMCMOL = Q(JPROB+21)
         CHCMOL = Q(JPROB+25)
         IKCUT  = Q(JMULOF+NEK1+1)
         STOPC  = Q(JMULOF+NEK1+2)
         IF(ISTRA.GT.0) THEN
            JTSTRA = LQ(JMA-19)
            JTSTCO = LQ(JTSTRA-1)
            JTSTEN = LQ(JTSTRA-2)
#if defined(CERNLIB_ASHO)
            IF(ISTRA.EQ.2) THEN
               JTASHO = LQ(JMA-20)
            ENDIF
#endif
         ENDIF
      ENDIF
*
* *** Compute current step size
*
      STEP   = STEMAX
      IPROC  = 103
      GEKRT1 = 1. -GEKRAT
*
*  **   Step limitation due to bremsstrahlung ?
*
      IF (IBREM.GT.0) THEN
         STEPBR = GEKRT1*Q(JBREM+IEKBIN) +GEKRAT*Q(JBREM+IEKBIN+1)
         SBREM  = STEPBR*ZINTBR
         IF (SBREM.LT.STEP) THEN
            STEP  = SBREM
            IPROC = 9
         ENDIF
      ENDIF
*
*  **   Step limitation due to delta-ray production ?
*
      IF (IDRAY.GT.0) THEN
         STEPDR = GEKRT1*Q(JDRAY+IEKBIN) +GEKRAT*Q(JDRAY+IEKBIN+1)
         SDRAY  = STEPDR*ZINTDR
         IF (SDRAY.LT.STEP) THEN
            STEP  = SDRAY
            IPROC = 10
         ENDIF
      ENDIF
*
*  **   Step limitation due to annihilation ?
*
      IF (CHARGE.GT.0.) THEN
         IF (IANNI.GT.0) THEN
            STEPAN = GEKRT1*Q(JANNI+IEKBIN) +GEKRAT*Q(JANNI+IEKBIN+1)
            SANNI  = STEPAN*ZINTAN
            IF (SANNI.LT.STEP) THEN
               STEP  = SANNI
               IPROC = 11
            ENDIF
         ENDIF
      ENDIF
*
      IF (STEP.LE.0.) THEN
         STEP = 0.
         GO TO 90
      ENDIF
*
*  **   Step limitation due to energy-loss,multiple scattering
*             or magnetic field ?
*
      IF (JMULOF.NE.0) THEN
         SMULOF  = GEKRT1*Q(JMULOF+IEKBIN) +GEKRAT*Q(JMULOF+IEKBIN+1)
         IF (SMULOF.LT.STEP) THEN
            STEP  = SMULOF
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
            DO 10 IST = IQ(JSTAK+3),IQ(JSTAK+1)
               JST    = JSTAK + 3 + (IST-1)*NWSTAK
               Q(JST+11) = SAFETY
   10       CONTINUE
            IQ(JSTAK+3) = 0
         ENDIF
      ELSE
         IQ(JSTAK+3) = 0
      ENDIF
*
      IF (IFIELD.EQ.0.OR.STEP.LE.PREC) THEN
*
* *** Linear transport when no field or very short step
*
         IF (IGNEXT.NE.0) THEN
*
* *** Particle is supposed to cross the boundary during step
*
            DO 20 I = 1,3
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
10000 FORMAT(' Boundary correction in GTELEC: ',
     +       '    GEKIN      NUMED       STEP      SNEXT')
10100 FORMAT(31X,E10.3,1X,I10,1X,E10.3,1X,E10.3,1X)
#endif
                  ENDIF
               ENDIF
               VECT(I) = VECTMP
   20       CONTINUE
            INWVOL = 2
            NMEC = NMEC +1
            LMEC(NMEC) = 1
         ELSE
            DO 30 I = 1,3
               VECT(I)  = VECT(I) +STEP*VECT(I+3)
   30       CONTINUE
         ENDIF
      ELSE
*
* *** otherwise, swim particle in magnetic field
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
   40    CALL GUSWIM (CHARGE, STEP, VECT, VOUT)
#endif
#if defined(CERNLIB_USRJMP)
   40    CALL JUMPT4(JUSWIM,CHARGE, STEP, VECT, VOUT)
#endif
*
*  ** When near to boundary, take proper action (cut-step,crossing...)
*
         IF(STEP.GE.SAFETY)THEN
            INEAR = 0
            IF (IGNEXT.NE.0) THEN
               DO 50 I = 1,3
                  VNEXT(I+3) = VECT(I+3)
                  VNEXT(I) = VECT(I) +SNEXT*VECT(I+3)
   50          CONTINUE
               DO 60 I = 1,3
                  IF (ABS(VOUT(I)-VNEXT(I)).GT.EPSIL) GO TO 70
   60          CONTINUE
               INEAR = 1
            ENDIF
*
   70       CALL GINVOL (VOUT, ISAME)
            IF (ISAME.EQ.0) THEN
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
                  GO TO 40
               ENDIF
            ENDIF
         ENDIF
*
*
         DO 80 I = 1,6
            VECT(I) = VOUT(I)
   80    CONTINUE
*
      ENDIF
*
*
* *** Correct the step due to multiple scattering
      IF (IMULL.NE.0) THEN
         STMULS = STEP
         CORR=0.0001*(STEP/RADL)*(GETOT/(VECT(7)*VECT(7)))**2
         IF (CORR.GT.0.25) CORR = 0.25
         STEP  = (1.+CORR)*STEP
      ENDIF
*
      SLENG = SLENG + STEP
*
*  **  Apply synchrotron radiation if required
*
      IF(ISYNC*IFIELD.NE.0) THEN
         CALL GSYNC
         NMEC = NMEC+1
         LMEC(NMEC) = 108
      ENDIF
*
* *** Generate Cherenkov photons if required
*
      IF(IMCKOV.EQ.1) THEN
         CALL GGCKOV
         IF(NGPHOT.NE.0) THEN
            NMEC = NMEC+1
            LMEC(NMEC)=105
         ENDIF
      ENDIF
*
* *** Apply energy loss : find the kinetic energy corresponding
*      to the new stopping range = stopmx - step
*
      IF (ILOSL.NE.0) THEN
         NMEC = NMEC +1
         IF(GEKRAT.LT.0.7) THEN
            I1 = MAX(IEKBIN-1,1)
         ELSE
            I1 = MIN(IEKBIN,NEKBIN-1)
         ENDIF
         I1 = 3*(I1-1)+1
         XCOEF1 = Q(JCOEF+I1)
         XCOEF2 = Q(JCOEF+I1+1)
         XCOEF3 = Q(JCOEF+I1+2)
         IF(XCOEF1.NE.0.) THEN
            STOPMX = -XCOEF2+SIGN(ONE,XCOEF1)*SQRT(XCOEF2**2 - (XCOEF3-
     +      GEKIN/XCOEF1))
         ELSE
            STOPMX = - (XCOEF3-GEKIN)/XCOEF2
         ENDIF
*
         STOPRG = STOPMX - STEP
         IF (STOPRG.LT.STOPC) THEN
            STEP = MAX(STOPMX - STOPC,ZERO)
            GO TO 100
         ENDIF
*
         LMEC(NMEC) = 3
         IF(XCOEF1.NE.0.) THEN
            DEMEAN=GEKIN-XCOEF1*(XCOEF3+STOPRG*(2.*XCOEF2+STOPRG))
         ELSE
            DEMEAN=GEKIN-XCOEF2*STOPRG-XCOEF3
         ENDIF
         IF(DEMEAN.LE.5.*GEKIN*EPSMAC) THEN
            DEMEAN=(GEKRT1*Q(JLOSS+IEKBIN)+GEKRAT*Q(JLOSS+IEKBIN+1))
     +     *STEP
         ENDIF
         IF (ILOSS.EQ.4.OR.IEKBIN.LE.IKCUT+1) THEN
            DESTEP = DEMEAN
         ELSE
            DEMS = DEMEAN
            CALL GFLUCT (DEMS,DESTEP)
         ENDIF
         DESTEP=MAX(DESTEP,0.)
         GEKINT = GEKIN -DESTEP
         IF (GEKINT.LE.(1.01*CUTELE)) GO TO 100
         DESTEL = DESTEP
         GEKIN  = GEKINT
         GETOT  = GEKIN +AMASS
         VECT(7)= SQRT((GETOT+AMASS)*GEKIN)
         CALL GEKBIN
*
*        IABAN = 1 does not distinguish between sensitive and
*        non-sensitive volumes, and can stop particles above the CUTE
*
         IF(IABAN.EQ.1) THEN
*
*        STOP electron/positron if range < safety AND no brems
*
           IF(STOPMX.LE.SAFETY) THEN
             IF(SBREM.GT.STOPMX) THEN
*          + condition in case of Cherenkov generation:
*          STOP if E/p >refractive index (i.e. e+/e- is below threshold)
               IF(IMCKOV.EQ.1) THEN
                 THRIND=GETOT/VECT(7)
*****            IF(THRIND.GT.Q(JINDEX+1)) GOTO 100
                 IF(THRIND.GT.Q(JINDEX+1)) GOTO 98 
               ELSE
                 GOTO 98 
               ENDIF
             ENDIF
           ENDIF
           STOPRG = STOPMX - STEP
           IF (STOPRG.LT.STOPC) THEN
            STEP = MAX(STOPMX - STOPC,ZERO)
            GO TO 98 
           ENDIF
         END IF
*
*        IABAN = 2 distinguishes between sensitive and non-sensitive
*        volumes.
*        In sensitive volumes additional tests are applied before
*        the particle is stopped
*
         IF(IABAN.EQ.2) THEN
           IF(ISVOL.LE.0) THEN
             IF(STOPMX.LE.SAFETY) THEN
*
*        test for brems and annihilation only
*
               IF((IBREM.GT.0).AND.(SBREM.LE.STOPMX)) GOTO 97
               IF((CHARGE.GT.0.).AND.(IANNI.GT.0).and.
     +                                 (SANNI.LE.STOPMX)) GOTO 97
               GOTO 98 
             END IF
           ELSE
*
*        sensitive volume ---> more tests !!
*        is energy below TLIM (=200 keV ) ?
*
             IF(GEKIN.LE.TLIM) THEN
*
*     range of the particle is the overestimated stopping range here
*
                TOSTOP=STOPMX-STEP*DESTEP/DEMEAN-STOPC
*
*        does the track remain in the actual volume?
*
                IF(TOSTOP.LE.SAFETY) THEN
*
*        is there no delta ray, brems, annihilation ?
*
                  IF((IDRAY.GT.0).AND.(SDRAY.LE.TOSTOP)) GOTO 97
                  IF((IBREM.GT.0).AND.(SBREM.LE.TOSTOP)) GOTO 97
                  IF((CHARGE.GT.0.).AND.(IANNI.GT.0).and.
     +                                   (SANNI.LE.TOSTOP)) GOTO 97
*
*        extra condition in case of Cherenkov generation:
*
                  IF(IMCKOV.EQ.1) THEN
                     THRIND=GETOT/VECT(7)
*
*        continue only if e+/e- below threshold
*
                     IF(THRIND.LT.Q(JINDEX+1)) GOTO 97
                   ELSE
*
*        do not make transport if this estimated range negative...
*
                     IF(TOSTOP.LE.0.) GOTO 98 
*
*        estimate final position/direction of the particle
*        from energy loss + multiple scattering
*       ( multiple scattering with path length = range of the particle)
*
                     ALFA=0.18*Z
                     ALFA1=ALFA+1.
                     TGTH2=ALFA*ALFA1/(1.2+ALFA)
                     S=TOSTOP*SQRT(1.+TGTH2)/ALFA1
                     TGTH=SQRT(TGTH2)
                     THET=ATAN(TGTH)
                     IF(THET.Lt.0.) THET=PI-THET
*
*        correct direction
*
                     CT=COS(THET)
                     ST=SQRT(1.-CT*CT)
                     CALL GRNDM(RNDM,1)
                     PHI=TWOPI*RNDM(1)
*
                     D1=ST*COS(PHI)
                     D2=ST*SIN(PHI)
                     D3=CT
                     VMM=SQRT(VECT(4)*VECT(4)+VECT(5)*VECT(5))
                     IF(VMM.NE.0.) THEN
                       PD1=VECT(4)/VMM
                       PD2=VECT(5)/VMM
                       V4=PD1*VECT(6)*D1-PD2*D2+VECT(4)*D3
                       V5=PD2*VECT(6)*D1+PD1*D2+VECT(5)*D3
                       V6=-VMM*D1+VECT(6)*D3
                     ELSE
                       V4=D1
                       V5=D2
                       V6=D3
                     ENDIF
                     VP=1./SQRT(V4*V4+V5*V5+V6*V6)
                     VECT(4)=V4*VP
                     VECT(5)=V5*VP
                     VECT(6)=V6*VP
*
*       transport particle ( assuming FIELD = 0. )
*
                     DO 123 I=1,3
                     VECT(I)=VECT(I)+S*VECT(I+3)
123                  CONTINUE
*
*       put back into GEKIN the original value in order to have
*       a correct DESTEP
*
                  GOTO 98 
                END IF
              END IF
            END IF
*++++++++++++++++++++++++++++++++++++++++++++++++++++++
          END IF
97        CONTINUE
        ENDIF
      ENDIF
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*
* *** Apply multiple scattering.
*
      IF (IMULL.NE.0) THEN
         NMEC = NMEC +1
         LMEC(NMEC) = 2
         CALL GMULTS
      ENDIF
*
* *** Update time of flight
*
      TOFG = TOFG +STEP*GETOT/(VECT(7)*CLIGHT)
      IF (TOFG.GE.TOFMAX) THEN
         ISTOP = 4
         NMEC  = NMEC +1
         LMEC(NMEC) = 22
         GO TO 999
      ENDIF
*
* *** Update interaction probabilities
*
      IF (IBREM.GT.0)    ZINTBR = ZINTBR -STEP/STEPBR
      IF (IDRAY.GT.0)    ZINTDR = ZINTDR -STEP/STEPDR
      IF (CHARGE.GT.0.) THEN
         IF (IANNI.GT.0) ZINTAN = ZINTAN -STEP/STEPAN
      ENDIF
*
* *** Apply the selected process if any
*
   90 IF (IPROC.EQ.0) GO TO 999
      NMEC = NMEC +1
      LMEC(NMEC) = IPROC
*
*  **   Bremsstrahlung ?
*
      IF (IPROC.EQ.9) THEN
         CALL GBREME
*
*  **   Delta ray ?
*
      ELSE IF (IPROC.EQ.10) THEN
*
       IF((IPART.EQ.2).OR.((IPART.EQ.3).AND.(GEKIN.GT.2.*DCUTE))) THEN
         CALL GDRAY
       ELSE
         GOTO 98 
       ENDIF
*
*  **   Positron annihilation ?
*
      ELSE IF (IPROC.EQ.11) THEN
         CALL GANNI
 
      ENDIF
      GO TO 999
*
* *** Special treatment for overstopped tracks
*
  98  GEKIN=GEKIN+DESTEP
  100 DESTEP = GEKIN
      DESTEL = DESTEP
      GEKIN  = 0.
      GETOT  = AMASS
      VECT(7)= 0.
      INWVOL = 0
      NMEC   = NMEC +1
      LMEC(NMEC) = 30
      IF ((CHARGE.LT.0.).OR.(IANNI.EQ.0)) THEN
         ISTOP = 2
      ELSE
         NMEC = NMEC +1
         LMEC(NMEC) = 11
         CALL GANNIR
      ENDIF
  999 IF(NGPHOT.GT.0) THEN
         IF(ITCKOV.EQ.2.AND.ISTOP.EQ.0) THEN
*
*  The electron has produced Cerenkov photons and it is still alive
*  we put it in the stack and we let the photons to be tracked
            NGKINE = NGKINE+1
            GKIN(1,NGKINE) = VECT(4)*VECT(7)
            GKIN(2,NGKINE) = VECT(5)*VECT(7)
            GKIN(3,NGKINE) = VECT(6)*VECT(7)
            GKIN(4,NGKINE) = GETOT
            GKIN(5,NGKINE) = IPART
            TOFD(NGKINE) = 0.
            ISTOP = 1
c----put position as well
            GPOS(1,NGKINE)=VECT(1)
            GPOS(2,NGKINE)=VECT(2)
            GPOS(3,NGKINE)=VECT(3)
         ENDIF
      ENDIF
      END

