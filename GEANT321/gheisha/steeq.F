*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE STEEQ(XXH,IPV)
C
C *** CORRECTIONS FOR SINGLE PARTICLE SPECTRA (SHOWER PARTICLES) ***
C *** NVE 16-MAR-1988 CERN GENEVA
C
C ORIGIN : H.FESEFELDT (06-SEP-1985)
C
#include "geant321/mxgkgh.inc"
#include "geant321/s_curpar.inc"
#include "geant321/s_blank.inc"
#include "geant321/s_result.inc"
C
      DIMENSION ALEM(7),VAL0(7)
      DIMENSION RNDM(1)
C**   DATA   EM/ 4.0 , 10.  , 15.  , 20.  ,  30. , 100. , 1000./
      DATA ALEM/ 1.40, 2.30 , 2.70 , 3.00 , 3.40 , 4.60 , 7.00 /
      DATA VAL0/ 0.00, 0.40 , 0.48 , 0.51 , 0.54 , 0.60 , 0.65 /
C
      XXH=1.
C
      IF ((IPART .NE. 7) .AND. (IPART .NE. 9)) GO TO 9999
      IF (ABS(IPA(IPV)) .NE. 8) GO TO 9999
      CALL GRNDM(RNDM,1)
      IF (RNDM(1) .GT. LOG(ATNO2)) GO TO 9999
      EKW=PV(4,MXGKPV)-ABS(PV(5,MXGKPV))
      ALEKW=LOG(EKW)
      IF (ALEKW .LE. ALEM(1)) GO TO 9999
C
C --- GET ENERGY BIN ---
      DO 1 I=2,7
      IF (ALEKW .LT. ALEM(I)) GO TO 2
 1    CONTINUE
      XXH=VAL0(7)
      GO TO 3
C
C *** USE LINEAR INTERPOLATION OR EXTRAPOLATION BY Y=RC*X+B ***
 2    CONTINUE
      I1=I-1
      I2=I
      DXNVE=ALEM(I2)-ALEM(I1)
      DYNVE=VAL0(I2)-VAL0(I1)
      RCNVE=DYNVE/DXNVE
      BNVE=VAL0(I1)-RCNVE*ALEM(I1)
      XXH=RCNVE*ALEKW+BNVE
C
 3    CONTINUE
      XXH=1.-XXH
C
 9999 CONTINUE
      RETURN
      END
