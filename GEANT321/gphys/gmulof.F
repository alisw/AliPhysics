#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GMULOF
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates table of steps for multiple scattering             *
C.    *     energy loss and magnetic field for electrons,muons         *
C.    *           (cannot be tabuled for hadrons)                      *
C.    *   : smuls  = min (Tbethe , 10*radl)                            *
C.    *   : sloss  = DEEMAX*GEKIN/DEDX                                 *
C.    *   : sfield = CFLD*P                                            *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Authors    R.Brun, Y.Dufour, M.Maire  *********          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gckine.inc"
#include "geant321/gcmate.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcking.inc"
#include "geant321/gctmed.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcphys.inc"
*
      LOGICAL CERKOV
*
*-----------------------------------------------------------------------
*
      SMULS  = BIG
      SLOSS  = BIG
      SFIELD = BIG
      STOPMX = BIG
      STCKOV = BIG
      JPROB  = LQ(JMA-4)
      JMIXT  = LQ(JMA-5)
      OMC    = Q(JPROB+21)
      CHC2   = Q(JPROB+25)**2
      NLMAT=Q(JMA+11)
      NLM=IABS(NLMAT)
      IF (FIELDM.NE.0.) CFLD = 3333.*DEGRAD*TMAXFD/ABS(FIELDM)
*
      IF(ITCKOV.NE.0.AND.IQ(JTM-2).GE.3.AND. LQ(JTM-3)
     +.NE.0.AND.LQ(LQ(JTM-3)-3).NE.0) THEN
*
* ***  In this tracking medium Cerenkov photons are generated and
* ***  tracked. Set to 1 the corresponding flag and calculate the
* ***  relevant pointers.
*
         CERKOV = .TRUE.
         JTCKOV = LQ(JTM-3)
         JABSCO = LQ(JTCKOV-1)
         JEFFIC = LQ(JTCKOV-2)
         JINDEX = LQ(JTCKOV-3)
         JCURIN = LQ(JTCKOV-4)
         NPCKOV = Q(JTCKOV+1)
      ELSE
         CERKOV = .FALSE.
      ENDIF
*
* *** Electrons
*
      JRANG = LQ(JMA-15)
      IKCUT = MAX((GEKA*LOG10(CUTELE) + GEKB),1.)
      GKC   = (CUTELE-ELOW(IKCUT))/(ELOW(IKCUT+1)-ELOW(IKCUT))
      STOPC = (1.-GKC)*Q(JRANG+IKCUT) + GKC*Q(JRANG+IKCUT+1)
      JMULOF = LQ(JTM-1)
      Q(JMULOF+NEK1+1) = IKCUT
      Q(JMULOF+NEK1+2) = STOPC
*
* *** Recompute STMIN ?
*     set STMIN to the range of an electron at energy=CUTELE + 200KeV
*     divided by sqrt(RADL) (important for light materials)
*
      IF(STMIN.LT.0.)THEN
         XES=CUTELE+2.E-4
         IKS = MAX((GEKA*LOG10(XES) + GEKB),1.)
         GKS   = (XES-ELOW(IKS))/(ELOW(IKS+1)-ELOW(IKS))
         STMIN = (1.-GKS)*Q(JRANG+IKS) + GKS*Q(JRANG+IKS+1) - STOPC
         IF(Q(JTM+7).EQ.0.)THEN
            STMIN = 2.*STMIN/SQRT(RADL)
         ELSE
            STMIN = 5.*STMIN/RADL
         ENDIF
      ENDIF
      Q(JTM+14)=STMIN
*
      DO 10 IEKBIN=1,NEK1
         GEKIN = ELOW(IEKBIN)
         GETOT = GEKIN + EMASS
         PMOM2 = GEKIN*(GETOT+EMASS)
         PMOM  = SQRT(PMOM2)
         BETA2 = PMOM2/(GETOT**2)
*
         IF (IMULS.GT.0.) THEN
            IF(JMIXT.LE.0)THEN
               CALL GMOLIO(A,Z,1.,1,DENS,BETA2,1.,OMC)
            ELSE
               CALL GMOLIO(Q(JMIXT+1),Q(JMIXT+NLM+1),Q(JMIXT+2*NLM+1),
     +                    NLM,DENS,BETA2,1.,OMC)
            ENDIF
            PMCH2  = PMOM2/CHC2
            TBETHE = (PMCH2*BETA2)/LOG(OMC*PMCH2)
            TMXCOR = 2232.*RADL*PMOM2*BETA2
            SMULS  = MIN(TBETHE,TMXCOR,10.*RADL)
         ENDIF
*
         IF (IFIELD*FIELDM.NE.0.) THEN
            SFIELD = CFLD*PMOM
         ENDIF
*
         IF (ILOSS*DEEMAX.GT.0.) THEN
            IF (IEKBIN.LE.IKCUT) THEN
               STOPMX = 0.
               SLOSS  = 0.
            ELSE
               STOPMX = Q(JRANG+IEKBIN)
               EKF = (1.-DEEMAX)*GEKIN
               IF (EKF.LE.ELOW(1)) EKF = ELOW(1)
               IKF = MAX((GEKA*LOG10(EKF) + GEKB),1.)
               GKR = (EKF-ELOW(IKF))/(ELOW(IKF+1)-ELOW(IKF))
               SLOSS = STOPMX-(1.-GKR)*Q(JRANG+IKF)-GKR*Q(JRANG+IKF+1)
               IF (SLOSS.LE.0.) SLOSS = 0.
               STOPMX = STOPMX-STOPC
               IF (STOPMX.LE.0.) STOPMX = 0.
            ENDIF
         ENDIF
         IF(CERKOV) THEN
            CHARGE = 1.
            VECT(7) = PMOM
            CALL GNCKOV
            STCKOV = MXPHOT/MAX(3.*DNDL,1E-10)
         ENDIF
*
         STEP = MIN(SMULS,SLOSS,SFIELD,STCKOV)
         IF (STEP.LT.STMIN) THEN
            STEP = MIN(STMIN,STOPMX)
         ENDIF
         Q(JMULOF+IEKBIN) = STEP
   10 CONTINUE
      DO 20 I=1,IKCUT
         Q(JMULOF+I)=0.5*Q(JMULOF+IKCUT+1)
   20 CONTINUE
*
* *** Muons
*
      JRANG = LQ(JMA-16)
      IKCUT = GEKA*LOG10(CUTMUO) + GEKB
      GKC   = (CUTMUO-ELOW(IKCUT))/(ELOW(IKCUT+1)-ELOW(IKCUT))
      STOPC = (1.-GKC)*Q(JRANG+IKCUT) + GKC*Q(JRANG+IKCUT+1)
      JMULOF = LQ(JTM-2)
      Q(JMULOF+NEK1+1)=IKCUT
      Q(JMULOF+NEK1+2)=STOPC
*
      DO 30 IEKBIN=1,NEK1
         GEKIN = ELOW(IEKBIN)
         GETOT = GEKIN + EMMU
         PMOM2 = GEKIN*(GETOT+EMMU)
         PMOM  = SQRT(PMOM2)
         BETA2 = PMOM2/(GETOT**2)
*
         IF (IMULS.GT.0.) THEN
            IF(JMIXT.LE.0)THEN
               CALL GMOLIO(A,Z,1.,1,DENS,BETA2,1.,OMC)
            ELSE
               CALL GMOLIO(Q(JMIXT+1),Q(JMIXT+NLM+1),Q(JMIXT+2*NLM+1),
     +                    NLM,DENS,BETA2,1.,OMC)
            ENDIF
            PMCH2  = PMOM2/CHC2
            TBETHE = (PMCH2*BETA2)/LOG(OMC*PMCH2)
            TMXCOR = 2232.*RADL*PMOM2*BETA2
            SMULS  = MIN(TBETHE,TMXCOR,10.*RADL)
         ENDIF
*
         IF (IFIELD*FIELDM.NE.0.) THEN
            SFIELD = CFLD*PMOM
         ENDIF
*
         IF (ILOSS*DEEMAX.GT.0.) THEN
            IF (IEKBIN.LE.IKCUT) THEN
               STOPMX = 0.
               SLOSS  = 0.
            ELSE
               STOPMX = Q(JRANG+IEKBIN)
               EKF = (1.-DEEMAX)*GEKIN
               IF (EKF.LE.ELOW(1)) EKF = ELOW(1)
               IKF = GEKA*LOG10(EKF) + GEKB
               GKR = (EKF-ELOW(IKF))/(ELOW(IKF+1)-ELOW(IKF))
               SLOSS = STOPMX-(1.-GKR)*Q(JRANG+IKF)-GKR*Q(JRANG+IKF+1)
               IF (SLOSS.LE.0.) SLOSS = 0.
               STOPMX = STOPMX-STOPC
               IF (STOPMX.LE.0.) STOPMX = 0.
            ENDIF
         ENDIF
         IF(CERKOV) THEN
            CHARGE = 1.
            VECT(7) = PMOM
            CALL GNCKOV
            STCKOV = MXPHOT/MAX(3.*DNDL,1E-10)
         ENDIF
*
         STEP = MIN(SMULS,SLOSS,SFIELD,STCKOV)
         IF (STEP.LT.STMIN) THEN
            STEP = MIN(STMIN,STOPMX)
         ENDIF
         Q(JMULOF+IEKBIN) = STEP
   30 CONTINUE
      DO 40 I=1,IKCUT
         Q(JMULOF+I)=0.5*Q(JMULOF+IKCUT+1)
   40 CONTINUE
*
      END
