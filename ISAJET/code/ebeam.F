#include "isajet/pilot.h"
      FUNCTION EBEAM(X,E)
C***********************************************************************
C* Modified from contributed subroutine by M. Drees (1/8/99)
C* Computes the effective single elctrn spectrum from beamstrahlung at *
C* e+e- colliders, using Chen's approximate expressions, for a given   *
C* beamstrahlung  parameter Y; is supposed to work for Y <= 10 or so.  *
C* The quantities in the COMMON block are the beamstrahlungs parameter *
C* Y, the bunch length XL in GeV, the number of photons NGAM, and the  *
C* parameters NUCL, NUGAM, W, XKAPPA defined by Chen, as well as the   *
C* pre-factor FAC. Y, E and XLMM are read in by BEAM when it is called *
C* for the first time, with INIT=1; in this first run the other para-  *
C* meters are then computed, and simply used in later calls with       *
C* INIT = 0. This COMMON block should be present in the main program   *
C* in order to guarantee the survival of these parameters. Finally, X  *
C* is the electron energy in units of the nominal beam energy. Notice  *
C* that BEAMEL is only the part which is NOT proportional to           *
C* delta(1-X); the coefficient of the delta-function is simply         *
C* (1-exp(-NGAM))/NGAM.                                                *
C***********************************************************************
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/eepar.inc"
C
      REAL XLMM,XL,GAM,RE,XKAPPA,NUCL,NUGAM,NGAM,X,NUBAR,
     $ETA,EPS,HFAC,RAT,MTERM,TERM,HBAR,XN,EBEAM,Y,E,GAMMA
      DOUBLE PRECISION DINCGM
C
      Y=UPSLON
      XLMM=SIGZ
      XLMM = 2.E0*SQRT(3.E0)*XLMM
      XL = XLMM*1.E12/.197327E0
      GAM = E/5.11E-4
      RE = 1.E0/(137.E0*5.11E-4)
      XKAPPA = 2.E0/(3.E0*Y)
      NUCL = 2.5E0*Y/(SQRT(3.E0)*137.E0**2*GAM*RE)
      NUGAM = NUCL/SQRT(1.E0+Y**.6666666E0)
      NGAM = .5E0*NUGAM*XL
      IF(X.LT.1.E-5) X=1.E-5
      IF(X.GT..99999) X=.99999
      NUBAR = X*NUCL + (1.0-X)*NUGAM
      ETA = XKAPPA*(1.0/X-1.0)
      IF(ETA.GT.5.E1) THEN
        EBEAM = 1.E-20
        RETURN
      ENDIF
      EPS = 1.E-4
      HFAC = EXP(-ETA)/(NGAM*(1.E0-X))
      IF(HFAC.LT.1.E0) EPS = EPS/SQRT(HFAC)
      RAT = NUBAR/NUGAM*(ETA**.33333333E0)
      MTERM = RAT
      TERM = MTERM/GAMMA(.33333E0)*DINCGM(2.D0,DBLE(NGAM),DBLE(EPS))
      HBAR = TERM
      XN = 1.0
1     XN = XN+1.0
      MTERM = MTERM*RAT/XN
      TERM = MTERM/GAMMA(XN/3.)*
     $DINCGM(DBLE(XN)+1.D0,DBLE(NGAM),DBLE(EPS))
      HBAR = HBAR+TERM
      IF(ABS(TERM/HBAR).GT.EPS) GO TO 1
      EBEAM = HFAC * HBAR
      IF(EBEAM.LT.0.) EBEAM = 1.E-20
      RETURN
      END
