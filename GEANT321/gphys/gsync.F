*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:34  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 11/10/94  19.56.01  by  S.Giani
*-- Author :
      SUBROUTINE GSYNC
C
C.    ******************************************************************
C.    *                                                                *
C.    *   Simulation of synchrotron radiation for e- and e+.           *
C.    *                                                                *
C.    *   ==> Called by : GTELEC                                       *
C.    *       Authors   : Igor Gavrilenko                              *
C.    *                                                                *
C.    ******************************************************************
C
#include "geant321/gcphys.inc"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctmed.inc"
#include "geant321/gccuts.inc"
C...
      DIMENSION ETR(100),STR(100),VECTS(7),VECTN(7),FLD(3)
C...
      IF(ABS(CHARGE).NE.1..OR.STEP.LE.0.) GO TO 999
      GAMMA  = GETOT/AMASS
      IF(GAMMA.LE.1000.                 ) GO TO 999
      KCASE = NAMEC1(8)
      IF (IFIELD.EQ.3) THEN
          FLD(1) = 0.
          FLD(2) = 0.
          FLD(3) = FIELDM
      ELSE
          CALL GUFLD(VECT,FLD)
      END IF
      FLDS2  = FLD(1)**2+FLD(2)**2+FLD(3)**2
      IF(FLDS2.LE.1E-10) GO TO 999
      FLDT   = SQRT(FLDS2-
     +  (VECT(4)*FLD(1)+VECT(5)*FLD(2)+VECT(6)*FLD(3))**2)
C...
      CALL GSYNGE(GAMMA,AMASS,FLDT,STEP,NTR,ETR,STR)
      IF(NTR.EQ.0) GO TO 999
C...
      ELOSS          = 0.
C
C     When SYNC=1 the generation point of the photons is the current
C     particle point
C     When SYNC=3, the generation point of the photon is distributed
C     along the real trajectory of the particle (not along the step
C     taken) computed by GUSWIM.
      DO 10 N = 1,NTR
         IF (ISYNC.NE.2.AND.ETR(N).GT.CUTGAM) THEN
            NGKINE = NGKINE+1
            IF(NGKINE.GT.MXGKIN) GO TO 20
            IF (ISYNC.EQ.1) THEN
                GKIN(1,NGKINE) = ETR(N)*VECT(4)
                GKIN(2,NGKINE) = ETR(N)*VECT(5)
                GKIN(3,NGKINE) = ETR(N)*VECT(6)
                GKIN(4,NGKINE) = ETR(N)
                GKIN(5,NGKINE) = 1.
                GPOS(1,NGKINE) = VECT(1)
                GPOS(2,NGKINE) = VECT(2)
                GPOS(3,NGKINE) = VECT(3)
                TOFD (NGKINE) = 0.
            ELSE IF (ISYNC.EQ.3) THEN
                VECTS(1) = VECT(1)
                VECTS(2) = VECT(2)
                VECTS(3) = VECT(3)
                VECTS(4) = -VECT(4)
                VECTS(5) = -VECT(5)
                VECTS(6) = -VECT(6)
                VECTS(7) = VECT(7)
                CALL GUSWIM(1.,STR(N),VECTS,VECTN)
                GKIN(1,NGKINE) = -ETR(N)*VECTN(4)
                GKIN(2,NGKINE) = -ETR(N)*VECTN(5)
                GKIN(3,NGKINE) = -ETR(N)*VECTN(6)
                GKIN(4,NGKINE) = ETR(N)
                GKIN(5,NGKINE) = 1.
                GPOS(1,NGKINE) = VECTN(1)
                GPOS(2,NGKINE) = VECTN(2)
                GPOS(3,NGKINE) = VECTN(3)
                TOFD (NGKINE) = 0.
             END IF
         ELSE
             DESTEP = DESTEP +ETR(N)
         END IF
         ELOSS = ELOSS+ETR(N)
   10 CONTINUE
C
C..Update the energy of the primary particle
C
   20 GETOT = GETOT-ELOSS
      GEKIN = GEKIN-ELOSS
      VECT(7) = SQRT(GEKIN*(GEKIN+2.*AMASS))
      CALL GEKBIN
C...
 999  END
