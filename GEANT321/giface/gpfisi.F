*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE GPFISI
C
C *** CALCULATION OF PHOTO-FISSION AND PHOTO-ABSORBTION         ***
C *** CROSS SECTIONS IN CURRENT MATERIAL                        ***
C *** HMF 25-AUG-1989 RWTH AACHEN / NVE 11-MAY-1990 CERN GENEVA ***
C
C CALLED BY : GPHYSI
C ORIGIN    : F.CARMINATI, H.FESEFELDT
C
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
C
      DIMENSION CSGN(19),CSGF(19),EKFISS(19)
      SAVE CSGN,CSGF,EKFISS
C
C     MEASURED CROSS SECTIONS ON U(238) IN MILLBARN, ENERGIES IN GEV.
C
      DATA CSGN  /   0.,  12.,  20.,  46., 115., 260., 450., 520.,
     +             400., 270., 200., 150., 120.,  80.,  55.,  27.,
     +              15.,   1.,   0./
      DATA CSGF  /   0.,   0.,   0.,  14.,  30.,  65., 100., 130.,
     +             150., 180., 180., 150., 120.,  80.,  55.,  27.,
     +              15.,   1.,   0./
      DATA EKFISS/0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,
     +            0.012,0.013,0.014,0.015,0.016,0.018,0.020,0.025,
     +            0.030,0.035,0.040/
C
      SIG=0.
      SIGFIS=0.
      IF(JPFIS.LE.0)                  GOTO 99
      IF(ELOW(IEKBIN).LE.EKFISS(1))   GOTO 90
      IF(ELOW(IEKBIN).GE.EKFISS(19))  GOTO 90
      IF(A.LT.230..OR.A.GT.240.)      GOTO 90
C
C     A-DEPENDENCE, MEASURED AT THE MAXIMUM OF THE RESONANCE
C
      Z43BA=-67.0+38.7*Z**(4./3.)/A
      IF(Z43BA.LE.0.)                 GOTO 90
C
C     ENERGY DEPENDENCE, AS TABULATED IN THE DATA STATEMENTS
C
      DO 10 I=2,19
        IF(ELOW(IEKBIN).LT.EKFISS(I)) GOTO 20
   10 CONTINUE
      I=19
   20 RAT = (EKFISS(I)-ELOW(IEKBIN))/(EKFISS(I)-EKFISS(I-1))
C
      SIGRES = CSGN(I-1)*RAT + CSGN(I)*(1.-RAT)
      SIGRES = SIGRES*Z43BA/0.53
      SIGRES = SIGRES*0.0006022045*DENS/A
C
      SIGFIS = CSGF(I-1)*RAT + CSGF(I)*(1.-RAT)
      SIGFIS = SIGFIS*Z43BA/0.53
      SIGFIS = SIGFIS*0.0006022045*DENS/A
C
      SIG    = SIGRES+SIGFIS
C
   90 IF(SIG.GT.0.) THEN
        Q(JPFIS+IEKBIN) = 1./SIG
      ELSE
        Q(JPFIS+IEKBIN) = BIG
      ENDIF
C
      IF(SIGFIS.GT.0.) THEN
        Q(JPFIS+NEK1+IEKBIN) = 1./SIGFIS
      ELSE
        Q(JPFIS+NEK1+IEKBIN) = BIG
      ENDIF
C
  99  END
