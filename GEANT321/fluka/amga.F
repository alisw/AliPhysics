*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE AMGA.FOR
*COPY AMGA
*
*=== amga =============================================================*
*
      FUNCTION AMGA(IT)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C*** RANDOM SELECTION OF MASSES OF DECAYING PARTICLES
      COMMON/FKABLT/AMM(110),GA(110),TAU(110),ICH(110)
     *,IBAR(110),K1(110),K2(110)
      COMMON /FKGAMR/ REDU,AMO,AM (15 )
      DIMENSION GASUNI(14)
      REAL RNDM(1)
      SAVE GASUNI, GAUNO, GAUNON, IO, NSTAB
      DATA GASUNI/
     *-1.D0,-.98D0,-.95D0,-.87D0,-.72D0,-.48D0,-.17D0,.17D0,.48D0,
     *.72D0,.87D0,.95D0,.98D0,1.D0/
      DATA GAUNO/2.352D0/
      DATA GAUNON/2.4D0/
      DATA IO/14/
      DATA NSTAB/23/
      I=1
      IF (IT.LE.0) GO TO 3
      IF (IT.LE.NSTAB) GO TO 2
      DGAUNI=GAUNO*GAUNON/(IO-1.D0)
      CALL GRNDM(RNDM,1)
      VV = RNDM(1)
      VV = VV*2.D0-1.D0+1.D-16
    1 CONTINUE
      VO=GASUNI(I)
      I=I+1
      V1=GASUNI(I)
      IF (VV.GT.V1) GO TO 1
      UNIGA=DGAUNI*(I-2.D0+(VV-VO+1.D-16)/(V1-VO)-(IO-1.D0)*.5D0)
      DAM=GA(IT)*UNIGA/GAUNO
      AAM=AMM(IT)+DAM
      AMGA=AAM
      RETURN
    2 CONTINUE
      AMGA=AMM(IT)
      RETURN
    3 CONTINUE
      AMGA=0.D0
      RETURN
      END
