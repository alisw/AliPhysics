*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      FUNCTION PMLTPC(NP,NM,NZ,N,B,C)
C
C *** NVE 03-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (14-SEP-1987)
C
#include "geant321/limits.inc"
      RLNNPF=0.
      IF(NP.LE.1) GOTO 2
      DO 1 I=2,NP
    1 RLNNPF=RLNNPF+LOG(I*1.)
    2 RLNNMF=0.
      IF(NM.LE.1) GOTO 4
      DO 3 I=2,NM
    3 RLNNMF=RLNNMF+LOG(I*1.)
    4 RLNNZF=0.
      IF(NZ.LE.1) GOTO 6
      DO 5 I=2,NZ
    5 RLNNZF=RLNNZF+LOG(I*1.)
    6 PMLTPC=-(NP-NM+NZ+B)**2/(2*(C*N)**2)-RLNNPF-RLNNMF-RLNNZF
      IF(PMLTPC.LT.EXPXL) PMLTPC=EXPXL
      PMLTPC=EXP(PMLTPC)
      RETURN
      END
