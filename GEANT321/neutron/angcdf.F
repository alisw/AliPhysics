*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE ANGCDF(D,LD,LZ)
C       THIS ROUTINE READS THE INPUT ANGULAR DISTRIBUTION FILES
C       AND CONVERTS THEM TO A NORMALIZED CDF
      DIMENSION D(*),LD(*)
      IPP=1
      NR=LD(IPP)
      NE=LD(IPP+1)
      NR2=2*NR
      II=2+NR2
   10 CONTINUE
      E=D(II+1)
      NP=LD(II+2)
      A1=-1.0
      PL=D(II+4)
      D(II+4)=0.0
      PROB=0.0
      DO 20 I=2,NP
         N=II+2*I+2
         A2=D(N-1)
         PH=D(N)
         PROB=PROB+(PH+PL)*(A2-A1)/2.0
         PL=PH
         D(N)=PROB
         A1=A2
   20 CONTINUE
      DO 30 I=1,NP
         N=II+2*I+2
         D(N)=D(N)/PROB
   30 CONTINUE
      II=II+2*NP+2
      IF(II.GE.LZ)GO TO 40
      GO TO 10
   40 RETURN
      END
