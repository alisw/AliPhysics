*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:34  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERTRAK (X1, P1, X2, P2, IPA, CHOPT)
*
************************************************************************
*                                                                      *
*          Perform the tracking of the track from point X1 to          *
*                    point X2                                          *
*          (Before calling this routine the user should also provide   *
*                    the input informations in /EROPTS/ and /ERTRIO/   *
*                    using subroutine EUFIL(L/P/V)                     *
*                 X1       - Starting coordinates (Cartesian)          *
*                 P1       - Starting 3-momentum  (Cartesian)          *
*                 X2       - Final coordinates    (Cartesian)          *
*                 P2       - Final 3-momentum     (Cartesian)          *
*                 IPA      - Particle code (a la GEANT) of the track   *
*                                                                      *
*                 CHOPT                                                *
*                     'B'   'Backward tracking' - i.e. energy loss     *
*                                        added to the current energy   *
*                     'E'   'Exact' calculation of errors assuming     *
*                                        helix (i.e. pathlength not    *
*                                        assumed as infinitesimal)     *
*                     'L'   Tracking upto prescribed Lengths reached   *
*                     'M'   'Mixed' prediction (not yet coded)         *
*                     'O'   Tracking 'Only' without calculating errors *
*                     'P'   Tracking upto prescribed Planes reached    *
*                     'V'   Tracking upto prescribed Volumes reached   *
*                     'X'   Tracking upto prescribed Point approached  *
*                                                                      *
*                Interface with GEANT :                                *
*             Track parameters are in /CGKINE/ and /GCTRAK/            *
*                                                                      *
*          ==>Called by : USER                                         *
*             Authors   M.Maire, E.Nagy  *********                     *
*                                                                      *
************************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gckine.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "geant321/ertrio.inc"
#include "geant321/erwork.inc"
#include "geant321/trcom3.inc"
*
      DIMENSION      P1(3), P2(3), X1(3), X2(3), DUM(15), IOPT(30)
      EQUIVALENCE    (IOPT(1),IOPTB), (IOPT(2),IOPTE), (IOPT(3),IOPTL),
     ,               (IOPT(4),IOPTM), (IOPT(5),IOPTO), (IOPT(6),IOPTP),
     ,               (IOPT(7),IOPTV), (IOPT(8),IOPTX)
      CHARACTER      CHOPT*(*)
*
*
* *** Decode character option
*
      CHOPTI = CHOPT
      CALL UOPTC (CHOPT, 'BELMOPVX', IOPT)
*
      IF (IOPTB.EQ.0) THEN
         BACKTR = 1.
      ELSE
         BACKTR = -1.
      ENDIF
*
      LEEXAC = IOPTE.NE.0
      LELENG = IOPTL.NE.0
      LEONLY = IOPTO.NE.0
      LEPLAN = IOPTP.NE.0
      LEVOLU = IOPTV.NE.0
      LEPOIN = IOPTX.NE.0
*
* *** Check consistency of the Ch-options
*
         IF ((LELENG .AND. LEVOLU) .OR. (LELENG .AND. LEPLAN) .OR.
     +       (LEVOLU .AND. LEPLAN)) THEN
            WRITE (LOUT, 779)
            GO TO 99
         ENDIF
*
* *** Initialization
*
      IF (NEPRED.LE.0) THEN
         WRITE (LOUT, 780)
         GO TO 99
      ENDIF
      ILPRED = 0
      TLGCM2 = 0.
      TLRAD  = 0.
      CALL VZERO (IEPRED ,   MXPRED)
      CALL VZERO (ERXOUT , 3*MXPRED)
      CALL VZERO (ERPOUT , 3*MXPRED)
      CALL VZERO (ERROUT ,15*MXPRED)
*
      DO 10 I = 1, 3
         ERXIN(I) = X1(I)
   10 CONTINUE
*
      PMOM2 = P1(1)**2 + P1(2)**2 + P1(3)**2
      IF(PMOM2.LE.1.E-20) THEN
         WRITE (LOUT, 778)
         GO TO 99
      ENDIF
      PABS = SQRT (PMOM2)
      ERPIN(1) = 1./PABS
      ERPIN(2) = ASIN (P1(3)*ERPIN(1))
      IF (ABS (P1(1)) .LT. 1.E-30) P1(1) = 1.E-30
      ERPIN(3) = ATAN2 (P1(2), P1(1))
*
* *** Initialize GCKINE common
*
      IF((IPA.LE.0).OR.(IPA.GT.NPART)) THEN
          WRITE (LOUT, 777) IPA
          GO TO 99
      ENDIF
*
      ITRA = 1
      ISTAK = 0
      IPART = IPA
      JPA = LQ(JPART-IPART)
      DO 26 I=1,5
         NAPART(I) = IQ(JPA+I)
   26 CONTINUE
      ITRTYP = Q(JPA+6)
      AMASS  = Q(JPA+7)
      CHARGE = Q(JPA+8)
      CHTR   = CHARGE*BACKTR
      TLIFE  = Q(JPA+9)
*
* *** Starting field
*
      CALL VZERO (HI, 9)
      CALL VZERO (HF, 9)
*
* *** Error matrix into SC System
*
      IF (LEPLAN) THEN
         IF (IFIELD.EQ.3) THEN
            HI(3) = FIELDM
         ELSEIF (IFIELD.NE.0) THEN
            CALL GUFLD (X1, HI)
         ENDIF
         CALL VZERO (DUM,15)
         CALL TRSCSD (ERPIN(1), DUM(1),   ERPIN(1), DUM(1),   HI(1),
     +                CHARGE, IERR, SPU, ERPLI(1,1), ERPLI(1,2))
         IF (LEONLY) GOTO 35
         CALL TRSDSC (ERPIN(1), ERRIN(1), DUM(1),   ERRIN(1), HI(1),
     +                CHARGE, IERR, SPU, ERPLI(1,1), ERPLI(1,2))
         DO 29 I = 1, 5
            DO 28 J = 1, 5
               ASDSC(I,J) = A(I,J)
   28       CONTINUE
   29    CONTINUE
      ENDIF
      IF (LEONLY) GOTO 35
*
* *** Error matrix into direction of tracking
*
      IF (BACKTR .LT. 0.) CALL ERBCER (ERRIN(1))
*
* *** Error matrix into double precision
*
      DO 30 I = 1,15
         EI(I) = ERRIN(I)
   30 CONTINUE
*
   35 CONTINUE
      DO 41 I = 1, 3
         VERT(I)  = X1(I)
         PVERT(I) = P1(I)*BACKTR
   41 CONTINUE
*
* *** Initialize GCTRAK common
*
      PVERT(4) = SQRT (PMOM2 + AMASS**2)
      VECT(7) = PABS
      DO 51 I=1,3
         VECT(I) = VERT(I)
         VECT(I+3) = PVERT(I)/VECT(7)
   51 CONTINUE
      GETOT = PVERT(4)
      GEKIN = GETOT - AMASS
      IF(GEKIN.LT.0.) GEKIN = 0.
      CALL GEKBIN
*
* *** Additional EMC initialization
*
      IF (.NOT.LEONLY) CALL ERPINI
*
* *** Ready for tracking
*
      CALL ERTRGO
*
* *** Copy the final point and momentum into the output buffer
*
      DO 60 I = 1, 3
         X2(I) = VECT(I)
         P2(I) = VECT(7)*BACKTR*VECT(I+3)
   60 CONTINUE
*
  777 FORMAT(/,4X,'Error in ERTRAK : particle type ', I4,
     *        '  unknown in GEANT' )
  778 FORMAT(/,4X,'Error in ERTRAK : Nul Momentum. Tracking stops now')
  779 FORMAT(/,4X,'Error in ERTRAK : Inconsistent character options',
     +                               '. Tracking stops now')
  780 FORMAT(/,4X,'Error in ERTRAK : No prediction. Tracking stops now')
*                                                           END ERTRAK
  99  END
