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
      SUBROUTINE GPHMFI
C
C *** CALCULATION OF PHOTO-FISSION AND PHOTOABSORBTION          ***
C *** CROSS SECTIONS IN CURRENT MATERIAL                        ***
C *** THE ENERGY BINS OF THE CROSS SECTION TABLES ARE NOT FINE  ***
C *** ENOUGH IN ORDER TO ACCURATELY DESCRIBE THE RESONANCE      ***
C *** BEHAVIOUR OF PHOTOFISSION AND PHOTOABSORBTION.            ***
C *** THE PRESENT ROUTINE CALCULATES THE CROSS SECTIONS AT      ***
C *** TRACKING TIME                                             ***
C *** HMF 25-AUG-1989 RWTH AACHEN / NVE 11-MAY-1990 CERN GENEVA ***
C
C CALLED BY : GTGAMA
C ORIGIN    : H.FESEFELDT
C
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gsecti.inc"
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
      IF(JPFIS.LE.0)                  GOTO 90
      IF(GEKIN.LE.EKFISS(1))          GOTO 90
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
        IF(GEKIN.LT.EKFISS(I))        GOTO 20
   10 CONTINUE
      I=19
   20 RAT     = (EKFISS(I)-GEKIN)/(EKFISS(I)-EKFISS(I-1))
      AIIN(1) = CSGN(I-1)*RAT + CSGN(I)*(1.-RAT)
      AIFI(1) = CSGF(I-1)*RAT + CSGF(I)*(1.-RAT)
      SIG     = AIIN(1) + AIFI(1)
C
      SIG = SIG*Z43BA/0.53
      SIG = SIG*0.0006022045*DENS/A
C
   90 IF(SIG.GT.0.) THEN
        ALAM = 1./SIG
      ELSE
        ALAM = BIG
      ENDIF
C
      END
