*CMZ :          24/11/95  16.28.16  by  S.Ravndal
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERTRCH
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Average charged track is extrapolated by one step           *
C.    *                                                                *
C.    *    ==>Called by : ERTRGO                                       *
C.    *       Original routine : GTHADR                                *
C.    *       Authors   M.Maire, E.Nagy  *********                     *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconst.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "geant321/ertrio.inc"
#include "geant321/erwork.inc"

#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
#if !defined(CERNLIB_SINGLE)&&!defined(CERNLIB_IBM)
      PARAMETER (EPSMAC=5.E-6)
#endif
#if !defined(CERNLIB_SINGLE)&&defined(CERNLIB_IBM)
      PARAMETER (EPSMAC=5.E-5)
#endif
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION GKR,DEMEAN,STOPP1,STOPP2,STOPMX,STOPRG,STOPC
      DOUBLE PRECISION EKIPR
#endif
      REAL VNEXT(6)
      SAVE CFLD,CHARG2,RMASS,CUTPRO,IKCUT,STOPC
C.
C.    ------------------------------------------------------------------
*
*
* *** Update local pointers if medium has changed
*
      IF (IUPD.EQ.0) THEN
         IUPD  = 1
         CHARG2 = CHARGE*CHARGE
         IF (IPART.LE.3) THEN
            CUTEK  = CUTELE
            RMASS  = 1.
            JRANG  = LQ(JMA-15)
         ELSE IF (IPART.LE.6) THEN
            CUTEK  = CUTMUO
            RMASS  = 1.
            JRANG  = LQ(JMA-16)
         ELSE
            CUTEK  = CUTHAD
            RMASS  = PMASS/AMASS
            JRANG  = LQ(JMA-16) + NEK1
         ENDIF
         CUTPRO = MAX(CUTEK*RMASS,ELOW(1))
         IKCUT = GEKA*LOG10(CUTPRO) + GEKB
         GKR   = (CUTPRO - ELOW(IKCUT))/(ELOW(IKCUT+1) - ELOW(IKCUT))
         STOPC = (1.-GKR)*Q(JRANG+IKCUT) + GKR*Q(JRANG+IKCUT+1)
         CFLD  = 0.
         IF (FIELDM.NE.0.) CFLD = 3333.*DEGRAD*TMAXFD/ABS(FIELDM*CHARGE)
      ENDIF
*
* *** Compute current step size
*
      STEP   = BIG
      GEKRT1 = 1. - GEKRAT
*
* *** Step limitation due to energy loss (stopping range) ?
*
      IF (ILOSS*DEEMAX.GT.0.) THEN
         STOPP1 = GEKRT1*Q(JRANG+IEKBIN) + GEKRAT*Q(JRANG+IEKBIN+1)
         STOPMX = (STOPP1 - STOPC)/(RMASS*CHARG2)
         EKF  = (1. - BACKTR*DEEMAX)*GEKIN*RMASS
         IF (EKF.LT.ELOW(1))    EKF = ELOW(1)
         IF (EKF.GE.ELOW(NEK1)) EKF = ELOW(NEK1)*0.99
         IKF=GEKA*LOG10(EKF)+GEKB
         GKR=(EKF-ELOW(IKF))/(ELOW(IKF+1)-ELOW(IKF))
         STOPP2 = (1.-GKR)*Q(JRANG+IKF) + GKR*Q(JRANG+IKF+1)
         SLOSP  = ABS (STOPP1 - STOPP2)
         STEP   = SLOSP/(RMASS*CHARG2)
      ENDIF
*
* *** Step limitation due to energy loss in magnetic field ?
*
      IF (IFIELD*FIELDM.NE.0.) THEN
         SFIELD = CFLD*VECT(7)
         IF (SFIELD.LT.STEP) STEP = SFIELD
      ENDIF
*
* *** Compute point where to store error matrix
*
      LERST  = 0
      STEPER = BIG
      ASCL1  = BIG
      DO 20 IPR = 1,NEPRED
         STEPE  = BIG
         IF (LELENG) STEPE = ERLENG(IPR) - SLENG
         IF (LEPLAN) THEN
            SCAL1 = 0.
            SCAL2 = 0.
            DO 18 I=1,3
               SCAL1 = SCAL1 + ERPLO(I,4,IPR)*(ERPLO(I,3,IPR)-VECT(I))
               SCAL2 = SCAL2 + ERPLO(I,4,IPR)*VECT(I+3)
   18       CONTINUE
            STEPE = SCAL1/SCAL2
         ENDIF
         IF (STEPE.LE.PREC) STEPE = BIG
         IF (STEPE.LT.STEPER) THEN
            STEPER = STEPE
            INLIST = IPR
            IF (LEPLAN) ASCL1  = ABS (SCAL1)
         ENDIF
   20 CONTINUE
      IF (STEPER.LE.STEP)  THEN
         STEP  = STEPER
         LERST = 1
      ENDIF
*
* *** Step limitation due to geometry ?
*
      LNEXT  = 0
      IF (STEP.GE.0.95*SAFETY) THEN
         CALL GTNEXT
         IF (IGNEXT.NE.0) THEN
            STEP   = SNEXT + PREC
            LNEXT = 1
            IF ((STEPER-SNEXT).GT.(2*PREC)) LERST = 0
         ENDIF
      ENDIF
*
* *** Linear transport when no field or very short step
*
      IF (IFIELD.EQ.0.OR.STEP.LE.2*PREC) THEN
        IF (IGNEXT.NE.0) THEN
          DO 25 I = 1,3
            VECTMP = VECT(I) +STEP*VECT(I+3)
            IF(VECTMP.EQ.VECT(I)) THEN
*
* *** Correct for machine precision
*
                  IF(VECT(I+3).NE.0.) THEN
                     VECTMP =
     +               VECT(I)+ABS(VECT(I))*SIGN(1.,VECT(I+3))*EPSMAC
#if defined(DEBUG)
                     WRITE(CHMAIL, 10000)
                     CALL GMAIL(0,0)
                     WRITE(CHMAIL, 10100) GEKIN, NUMED, STEP, SNEXT
                     CALL GMAIL(0,0)
10000 FORMAT(' Boundary correction in ERTRCH: ',
     +       '    GEKIN      NUMED       STEP      SNEXT')
10100 FORMAT(31X,E10.3,1X,I10,1X,E10.3,1X,E10.3,1X)
#endif
                  ENDIF
            ENDIF
            VOUT(I) = VECTMP
   25     CONTINUE
            INWVOL = 2
            NMEC = NMEC +1
            LMEC(NMEC) = 1
        ELSE
            DO 30 I = 1,3
               VOUT(I)  = VECT(I) +STEP*VECT(I+3)
   30       CONTINUE
        ENDIF
        DO 35 I = 4,6
          VOUT(I)  = VECT(I)
   35   CONTINUE
        GOTO 74
      END IF
*
* *** otherwise, swim particle in magnetic field
*
      NMEC = NMEC +1
      LMEC(NMEC) = 4
*
   50 LERST = 0
      LNEXT = 0
      CALL GUSWIM (CHTR , STEP, VECT, VOUT)
*
*     When near to boundary, take proper action (cut-step,crossing...)
      IF (STEP.GE.SAFETY) THEN
         INEAR = 0
         IF (IGNEXT.NE.0) THEN
           DO 51 I = 1,3
               VNEXT(I+3) = VECT(I+3)
               VNEXT(I) = VECT(I) +SNEXT*VECT(I+3)
   51      CONTINUE
           DO 52 I = 1,3
             IF (ABS(VOUT(I)-VNEXT(I)).GT.EPSIL) GO TO 55
   52      CONTINUE
           INEAR = 1
         ENDIF
*
   55    CALL GINVOL (VOUT,ISAME)
         IF (ISAME.EQ.0) THEN
           IF ((INEAR.NE.0).OR.(STEP.LT.EPSIL)) THEN
             INWVOL = 2
             NMEC = NMEC +1
             LMEC(NMEC) = 1
             LNEXT = 1
           ELSE
*              Cut step
             STEP = 0.5*STEP
             IF (LMEC(NMEC).NE.24) THEN
               NMEC = NMEC +1
               LMEC(NMEC) = 24
             ENDIF
             GOTO 50
           ENDIF
         ENDIF
      ENDIF
*
*
*     preset plane reached  ?
   74 CONTINUE
      IF ((LEPLAN).AND.(STEP.GE.ASCL1)) THEN
         SCAL3 = 0.
         DO 28 I=1,3
            SCAL3=SCAL3+ERPLO(I,4,INLIST)*(ERPLO(I,3,INLIST)-VOUT(I))
   28    CONTINUE
         ASCL3 = ABS(SCAL3)
         SSCL1 = ASCL1/SCAL1
         IF (SCAL3*SSCL1.LT. -PREC) THEN
*            Cut step
             STEP  = STEP*(ASCL1/(ASCL1+ASCL3))
             NMEC  = NMEC +1
             LMEC(NMEC) = 24
             GOTO 50
         ELSE
           IF(ASCL3.LE.PREC) LERST = 1
         ENDIF
      ENDIF
*
      DO 75 I=1,6
           VECT(I) = VOUT(I)
   75 CONTINUE
*
      IF (LELENG.AND.(STEP.GE.STEPER)) LERST = 1
*
      SLENG = SLENG + STEP
*
* *** Now apply selected mechanisms
*
      IF (LNEXT.EQ.1) THEN
          INWVOL = 2
          NMEC = NMEC + 1
          LMEC(NMEC) = 1
      ENDIF
*
* *** apply energy loss : find the kinetic energy corresponding
*      to the new stopping range = stopmx -/+ step
*      (take care of the back tracking !)
*
      IF (ILOSS*DEEMAX.GT.0) THEN
         NMEC = NMEC +1
         LMEC(NMEC) = 3
         CALL ERLAND (STEP,Z,A,DENS,VECT(7),GETOT,AMASS,DEDX2)
         DEDX2  = DEDX2*CHARG2*CHARG2
         STOPRG = STOPP1 - BACKTR*STEP*RMASS*CHARG2
         IKF = IEKBIN
         IF (BACKTR.LE.0.) THEN
   95       IF (STOPRG.LT.Q(JRANG+IKF)) THEN
               IKF = IKF - 1
               IF (IKF.GT.1) GO TO 95
            ENDIF
         ELSE
   96       IF (STOPRG.GE.Q(JRANG+IKF+1)) THEN
               IKF = IKF + 1
               IF (IKF.LT.NEK1) GO TO 96
            ENDIF
         ENDIF
         GKR = (STOPRG - Q(JRANG+IKF)) / (Q(JRANG+IKF+1) - Q(JRANG+IKF))
         EKIPR = (1. -GKR)*ELOW(IKF) + GKR*ELOW(IKF+1)
         GEKINT =  EKIPR/RMASS
         IF (GEKINT.GT.CUTEK) THEN
            DESTEP = ABS (GEKIN - GEKINT)
            GEKIN  = GEKINT
            GETOT  = GEKIN + AMASS
            VECT(7)= SQRT((GETOT+AMASS)*GEKIN)
            CALL GEKBIN
         ELSE
            DESTEP = GEKINT
            GEKIN  = 0.
            GETOT  = AMASS
            VECT(7)= 0.
            INWVOL = 0
            ISTOP  = 2
            NMEC = NMEC + 1
            LMEC(NMEC) = 30
         ENDIF
      ENDIF
*
* *** Propagate error matrix
*
      IF (.NOT. LEONLY) CALL ERPROP
*
* *** Store informations
*
      IF(LERST.EQ.1) THEN
         NMEC = NMEC + 1
         LMEC(NMEC) = 27
         CALL ERSTOR
      ENDIF
*
      END
