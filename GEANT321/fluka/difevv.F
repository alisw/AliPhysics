*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.42  by  S.Giani
*-- Author :
*$ CREATE DIFEVV.FOR
*COPY DIFEVV
*
*=== difevv ===========================================================*
*
      SUBROUTINE DIFEVV ( NHAD, KPROJ, KTARG, PPROJ, EPROJ, UMO )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*----------------------------------------------------------------------*
C
C          DIFFRACTIVE HADRON -HADRON COLLISIONS
C     GENERATE HADRON PRODUCTION EVENT IN  KPROJ - KTARG  COLLISION
C     WITH  LAB PROJECTILE MOMENTUM  PPROJ
C
C********************************************************************
C
#include "geant321/auxpar.inc"
#include "geant321/balanc.inc"
#include "geant321/cmsres.inc"
#include "geant321/finpar.inc"
#include "geant321/hadpar.inc"
#include "geant321/inpdat2.inc"
#include "geant321/part.inc"
#include "geant321/qquark.inc"
      COMMON /FKPRIN/ IPRI, INIT
      REAL RNDM(2)
C
C*******************************************************************"
C
C     KINEMATICS
C
C********************************************************************
C
      IPRI = 0
      AMPROJ = AM(KPROJ)
      AMTAR  = AM(KTARG)
* The following are the Lorentz parameters to come from the system
* (projectile + target) rest frame to the starting one, which is the
* one where the target is at rest and the projectile is moving
* along the +z direction with Pproj: from now down to 600 continue
* we are working in the system rest frame !!!
      GAMCM = (EPROJ+AMTAR)/UMO
      BGCM  = PPROJ/UMO
C
*or      IF(IPRI.EQ.1) WRITE(LUNOUT,101)KPROJ,KTARG,PPROJ,AMPROJ,AMTAR,
*or     *EPROJ,UMO,GAMCM,BGCM
*or101   FORMAT(2I5,10F11.5)
C
C
      IBPROJ = IBAR(KPROJ)
      IBTARG = IBAR(KTARG)
C
C
C=====================================================================
C
C     SAMPLE X-VALUES OF QUARK-ANTIQUARK PAIRS
C
C======================================================================
      IF ( KPROJ .GT. 2 ) THEN
         UNOSEA = 5.D+00
      ELSE
         UNOSEA = 3.D+00
      END IF
* Come here if we need to resample xsea and xasea!!
  211 CONTINUE
      TMP005 = 0.05D+00
      XSEA  = BETARN(TMP005,UNOSEA)
      XASEA = BETARN(TMP005,UNOSEA)
      XPIO  = XSEA+XASEA
      IF (XPIO .GE. 1.D+00) GO TO 211
      XHAD = 1.D+00 - XPIO
      CALL GRNDM(RNDM,1)
      ISAM = 2.D+00 * RNDM(1) + 1.D+00
*or      IF (IPRI.EQ.1) WRITE(LUNOUT,371)XSEA,XASEA,XPIO,XHAD,ISAM
*or  371 FORMAT (' XSEA,XASEA,XPIO,XHAD,ISAM',4F10.5,I10)
C=====================================================================
C
C      CALL EQUIVALENT PIO HADRON COLLISIONS
C
C=====================================================================
      GO TO (250,260),ISAM
*  +-------------------------------------------------------------------*
*  | Target excited !!!!
  250 CONTINUE
         LPRDIF = .TRUE.
C=======================================================================
C     PROJECTILE MOVING WITH XHAD, TARGET EXCITED
C=======================================================================
         IIDIF = 1
         AMCH  = SQRT (XPIO) * UMO
         BITBIT = 0.5D+00
*  |  +----------------------------------------------------------------*
*  |  | The following condition roughly guarantees 1 GeV for the total
*  |  | energy of the pseudo-pion!
         IF ( AMCH .LE. AMTAR + BITBIT ) THEN
            AMCH = AMTAR + BITBIT
            XPIO = ( AMCH / UMO )**2
            XHAD = 1.D+00 - XPIO
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  | The following instructions make the division of the invariant
*  | mass of the system between the two particles, the two resulting
*  | energies being Eh1s and Eh2ex and the momentum Ph1s: the two
*  | particles are the original projectile and the excited target,
*  | ("mass" = amch)
         EH1S = (UMO**2 + AMPROJ**2 - AMCH**2) / (2.D0*UMO)
         IF (EH1S .LE. AMPROJ) GO TO 211
         EH2EX = UMO - EH1S
         PH1S  = SQRT (EH1S**2 - AMPROJ**2)
C*** AND  INT. CHAIN TRANSVERSE MOMENTA
         B3SAVE = B3BAMJ
         CALL GRNDM(RNDM,2)
         ES  = -2.D0/(B3BAMJ**2)*LOG(RNDM(1)*RNDM(2))
         HPS = SQRT(ES*ES+2.D0*ES*AMTAR)
         CALL SFECFE(SFE,CFE)
         PTXCH1 = HPS * CFE
         PTYCH1 = HPS * SFE
6171  CONTINUE
         PX1 = PTXCH1
         PY1 = PTYCH1
         ACH = PH1S**2 - PX1**2 - PY1**2
         IF (ACH .LE. 0.D+00) THEN
            PTXCH1 = 0.75D+00 * PTXCH1
            PTYCH1 = 0.75D+00 * PTYCH1
            GO TO 6171
         END IF
         PZ1  = SQRT (ACH)
*  | Now transform back the excited target to the lab system
         ECHCK  = GAMCM * EH2EX - BGCM * PZ1
         PXCHCK = - PX1
         PYCHCK = - PY1
         PZCHCK = - GAMCM * PZ1 + BGCM * EH2EX
*  | Now ..chck are the kinematical variables of the excited target
*  | in the lab frame, the invariant mass is always Amch
         CALL GRNDM(RNDM,1)
         IF (RNDM (1) .LE. 0.5D+00 ) THEN
            KPIO = 26
         ELSE
            KPIO = 23
         END IF
         AMPIO = AM (KPIO)
         EPIOL = 0.5D+00 * ( AMCH**2 - AMPIO**2 - AMTAR**2 ) / AMTAR
         PPIOL = SQRT ( EPIOL**2 - AMPIO**2 )
         ETOTX  = EPIOL + AMTAR
         AAFACT = ECHCK + ETOTX
         BBFACT = PPIOL - PZCHCK
         DDENOM = ETOTX * AAFACT - PPIOL * BBFACT
         GAM1 = ( ECHCK * AAFACT + PPIOL * BBFACT ) / DDENOM
         BGZ1 = - BBFACT * AAFACT / DDENOM
         BGX1 = PXCHCK * ( GAM1 + 1.D+00 ) / AAFACT
         BGY1 = PYCHCK * ( GAM1 + 1.D+00 ) / AAFACT
         CALL HADEVV ( NHAD, KPIO, KTARG, PPIOL, EPIOL, AMCH )
*   Restore the original b3bamj parameter
         B3BAMJ = B3SAVE
C     PRINT 888,PX1,PY1,PZ1,EH2EX
C        PRINT 888,PXX1,PYY1,PZZ1,EXXX
C 888    FORMAT(4F12.4)
*  | The following to go back to the original (lab) frame
*  |  +----------------------------------------------------------------*
*  |  |                      Looping over the produced particles
         DO 800 I=1,NHAD
            CALL ALTRA ( GAM1, BGX1, BGY1, BGZ1, PXH(I), PYH(I), PZH(I),
     &                   HEPH(I), PLR, PLRX, PLRY, PLRZ, ELR )
            PXH(I) = PLRX
            PYH(I) = PLRY
            PZH(I) = PLRZ
            HEPH(I) = ELR
*  |  |  Updating conservation counters
 800     CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  |  Add the original projectile to the final particles, transforming
*  |  it back to the lab system
         NHAD = NHAD+1
         PXH  (NHAD) = PX1
         PYH  (NHAD) = PY1
         PZH  (NHAD) = GAMCM * PZ1 + BGCM * EH1S
         HEPH (NHAD) = GAMCM * EH1S + BGCM * PZ1
         AMH  (NHAD) = AMPROJ
         ICHH (NHAD) = ICH  (KPROJ)
         IBARH(NHAD) = IBAR (KPROJ)
         ANH  (NHAD) = ANAME(KPROJ)
         NREH (NHAD) = KPROJ
      GO TO 600
*  |  end of excited target treatment !!
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |  Excited projectile !!!!!!
 260  CONTINUE
         LPRDIF = .FALSE.
C=======================================================================
C   THE TARGET PARTICLE GETS XHAD , THE PROJECTILE BECOMES EXCITED
C   WE GO TO THE PROJECTILE REST FRAME
C=======================================================================
         IIDIF = 2
         AMCH  = SQRT (XPIO) * UMO
         MK = 0
         DO 270 IQ = 1, 3
            MK = MK + ABS ( IQSCHR ( MQUARK ( IQ, KPTOIP(KPROJ) ) ) )
 270     CONTINUE
         BITBIT = 0.5D+00 + 0.2D+00 * MK
*  |  +----------------------------------------------------------------*
*  |  | The following condition roughly guarantees 1 GeV for the total
*  |  | energy of the pseudo-pion if the projectile has no strangeness
*  |  | a bit more if it has
         IF ( AMCH .LE. AMPROJ + BITBIT ) THEN
            AMCH = AMPROJ + BITBIT
            XPIO = ( AMCH / UMO )**2
            XHAD = 1.D+00 - XPIO
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  | The following instructions make the division of the invariant
*  | mass of the system between the two particles, the two resulting
*  | energies being Eh1s and Eh2ex and the momentum Ph1s: the two
*  | particles are the target nucleon and the excited projectile,
*  | ("mass" = amch)
         EH1S = (UMO**2 + AMTAR**2 -  AMCH**2) / (2.D+00*UMO)
         IF (EH1S .LE. AMTAR) GO TO 211
         EH2EX = UMO-EH1S
         PH1S  = SQRT (EH1S**2 - AMTAR**2)
C*** AND  INT. CHAIN TRANSVERSE MOMENTA
         B3SAVE = B3BAMJ
         CALL GRNDM(RNDM,2)
         ES  = -2.D0/(B3BAMJ**2)*LOG(RNDM(1)*RNDM(2))
         HPS = SQRT(ES*ES+2.D0*ES*AMPROJ)
         CALL SFECFE(SFE,CFE)
         PTXCH1 = HPS * CFE
         PTYCH1 = HPS * SFE
6181  CONTINUE
         PX1 = PTXCH1
         PY1 = PTYCH1
         ACH = PH1S**2 - PX1**2 - PY1**2
         IF (ACH .LE. 0.D+00) THEN
            PTXCH1 = 0.75D+00 * PTXCH1
            PTYCH1 = 0.75D+00 * PTYCH1
            GO TO 6181
         END IF
         PZ1  = SQRT (ACH)
*  | Now transform back the excited projectile to the lab system
         ECHCK  = GAMCM * EH2EX + BGCM * PZ1
         PXCHCK = PX1
         PYCHCK = PY1
         PZCHCK = GAMCM * PZ1   + BGCM * EH2EX
*  | Now ..chck are the kinematical variables of the excited projectile
*  | in the lab frame, the invariant mass is always Amch
         CALL GRNDM(RNDM,1)
         IF (RNDM (1) .LE. 0.5D+00 ) THEN
            KPIO = 26
         ELSE
            KPIO = 23
         END IF
         AMPIO = AM (KPIO)
         EPIOL = 0.5D+00 * ( AMCH**2 - AMPIO**2 - AMPROJ**2 ) / AMPROJ
         PPIOL = SQRT ( EPIOL**2 - AMPIO**2 )
         ETOTX  = EPIOL + AMPROJ
         AAFACT = ECHCK + ETOTX
         BBFACT = PPIOL - PZCHCK
         DDENOM = ETOTX * AAFACT - PPIOL * BBFACT
         GAM1 = ( ECHCK * AAFACT + PPIOL * BBFACT ) / DDENOM
         BGZ1 = - BBFACT * AAFACT / DDENOM
         BGX1 = PXCHCK * ( GAM1 + 1.D+00 ) / AAFACT
         BGY1 = PYCHCK * ( GAM1 + 1.D+00 ) / AAFACT
         CALL HADEVV ( NHAD, KPIO, KPROJ, PPIOL, EPIOL, AMCH )
*   Restore the original b3bamj parameter
         B3BAMJ = B3SAVE
C     PRINT 888,PX1,PY1,PZ1,EH2EX
C        PRINT 888,PXX1,PYY1,PZZ1,EXXX
*  | The following to go back to the original (lab) frame
*  |  +----------------------------------------------------------------*
*  |  |                      Looping over the produced particles
         DO 900 I=1,NHAD
            CALL ALTRA ( GAM1, BGX1, BGY1, BGZ1, PXH(I), PYH(I), PZH(I),
     &                   HEPH(I), PLR, PLRX, PLRY, PLRZ, ELR )
            PXH(I) = PLRX
            PYH(I) = PLRY
            PZH(I) = PLRZ
            HEPH(I) = ELR
*  |  |  Updating conservation counters
 900     CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  |  Add the target nucleon to the final particles, transforming
*  |  it back to the lab system
         NHAD = NHAD + 1
         PXH  (NHAD) = - PX1
         PYH  (NHAD) = - PY1
         PZH  (NHAD) = - GAMCM * PZ1  + BGCM * EH1S
         HEPH (NHAD) =   GAMCM * EH1S - BGCM * PZ1
         AMH  (NHAD) = AMTAR
         ICHH (NHAD) = ICH  (KTARG)
         IBARH(NHAD) = IBAR (KTARG)
         ANH  (NHAD) = ANAME(KTARG)
         NREH (NHAD) = KTARG
      GO TO 600
*  |  end of excited projectile !!!
*  +-------------------------------------------------------------------*
  600 CONTINUE
C
C********************************************************************
C
C*** PRINT AND TEST ENERGY CONSERVATION
C
C********************************************************************
C
      PUZZ = 0.D+00
      EUZZ = 0.D+00
      PUXX = 0.D+00
      PUYY = 0.D+00
      ICUU = 0
      IBUU = 0
      DO 82 I=1,NHAD
         PUXX = PUXX + PXH(I)
         PUYY = PUYY + PYH(I)
         PUZZ = PUZZ + PZH(I)
         EUZZ = EUZZ + HEPH(I)
         ICUU = ICUU + ICHH(I)
         IBUU = IBUU + IBARH(I)
  82  CONTINUE
      ICHTOT=ICH(KPROJ)+ICH(KTARG)
      IBTOT =IBAR(KPROJ)+IBAR(KTARG)
      PCHMIN = 1.D-10 * PPROJ
      IF ((ABS(PUXX) .GE. PCHMIN) .OR. (ABS(PUYY) .GE. PCHMIN) .OR.
     &    (ABS(PUZZ-PPROJ) .GE. PCHMIN) .OR. (ABS(EPROJ+AMTAR-EUZZ) .GE.
     &     1.D-10*EUZZ) .OR. (ICHTOT .NE. ICUU) .OR. (IBTOT .NE. IBUU))
     &   THEN
         WRITE(LUNERR,*)
     &               ' Difevt: failure!!!: NHAD, KPROJ, KTARG, IIDIF: ',
     &                 NHAD, KPROJ, KTARG, IIDIF
         WRITE(LUNERR,*)'                     ',
     &                  'ICHTOT, ICUU, IBTOT, IBUU: ',
     &                   ICHTOT, ICUU, IBTOT, IBUU
         WRITE(LUNERR,*)'                     ',
     &                  'PPROJ, PUXX, PUYY, PUZZ: ',
     &                   PPROJ, PUXX, PUYY, PUZZ
         WRITE(LUNERR,*)'                     EPROJ, EUZZ: ',
     &                   EPROJ, EUZZ
      END IF
      IF (IPRI .NE. 1) GO TO 90
      DO 84 I=1,NHAD
         WRITE(LUNERR,85)I,NREH(I),ICHH(I),IBARH(I),ANH(I),
     &                   PXH(I),PYH(I),PZH(I),HEPH(I),AMH(I)
  85     FORMAT (4I5,A8,5F12.6)
  84  CONTINUE
  90  CONTINUE
      RETURN
      END
