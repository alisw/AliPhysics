*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:30:27  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:45  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE SFMSFG(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
c  This program returns the fitted parton distributions for MRS(G)   c
c  ... it is an analytic approximation to the "grid version" and     c
c  gives much smoother distributions. It is, however, much slower.   c
c  Full details of the fit are contained in the preprint "Pinning    c
c  Down the Gluon in the Proton", by A.D. Martin, R.G. Roberts and   c
c  W.J. Stirling, Rutherford preprint RAL-95-021 (1995), to which    c
c  reference should be made if appropriate.                          c
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
#include "pdf/impdp.inc"
C
      DIMENSION Z(18,8),ZZ(18),F(8)
      DATA ((Z(I,K),I=1,18),K=1,4) /
     +1.665297, -0.590753,  0.086807,  0.676836, -0.027729,  0.006186,
     +3.340307,  0.260426,  0.024214, -0.787023,  1.101313, -0.302934,
     +6.098057, -2.645910,  0.175731,  0.557278,  0.363394, -0.001549,

     +0.171268, -0.086429,  0.012191,  0.388998,  0.139088, -0.014565,
     +4.049580, -0.495801,  0.071001,  8.829433,  0.246082, -0.319317,
     +5.310913, -6.033850,  1.335430,  0.391294,  0.974633, -0.026356,

     +1.121933, -0.980891,  0.250001, -0.279071, -0.118949,  0.001370,
     +5.733053,  0.962694, -0.537345, -4.338513,  4.055084, -1.249222,
     +11.591153, -9.385871,  2.080363, 0.233747,  0.649513, -0.033386,

     +0.811684, -0.732591,  0.188398, -0.146012, -0.450745,  0.162278,
     +10.716037,  0.423859, -0.332372, -3.347546, 1.476767, -0.164336,
     +10.221619, -2.432161,  0.693669, -0.372843, -0.451730, 0.548307 /

      DATA ((Z(I,K),I=1,18),K=5,8) /
     +0.018819, -0.019057,  0.014581,  0.473385, -0.026701, -0.191575,
     +8.942582, -0.562553,  2.098186, -1.302327, 29.018376,-11.868145,
     +105.320622,-53.953297,  1.329035, 0.878000, 0.499894, -0.971650,

     +0.093684, -0.074091,  0.033660, -0.107765, -0.408746,  0.100063,
     +9.275198,  0.692364,  0.110208, -2.787409,  0.041578,  0.488094,
     +10.637883, 3.815478, -4.848576,  0.184162,  0.185983, -0.186796,

     +0.001026,  0.015587, -0.001414, -0.111311, -0.346509,  0.094542,
     +5.764679,  4.958334, -2.102611, -1.341127, -2.951861,  1.855254,
     +12.079860, -0.297841, -3.809543, 1.480855, -0.719267,  0.175987,

     +-0.001533,  0.004429, -0.000718, -0.097267, -0.221771, 0.056117,
     +3.630148,  2.941006, -0.512863, -0.483095,  1.599856, -1.197675,
     +8.366443, -4.004534,  1.022508,  2.526478, -0.713666,  0.277209/

      AL=0.254
      Q2=SCALE*SCALE
      Q0=2.
      ARG=( LOG(SCALE/AL)/ LOG(Q0/AL))
      S= LOG(ARG)

      DO 10 J=1,8

      A=Z(1,J)+Z(2,J)*S+Z(3,J)*S*S
      B=Z(4,J)+Z(5,J)*S+Z(6,J)*S*S
      C=Z(7,J)+Z(8,J)*S+Z(9,J)*S*S
      D=Z(10,J)+Z(11,J)*S+Z(12,J)*S*S
      E=Z(13,J)+Z(14,J)*S+Z(15,J)*S*S
      H=Z(16,J)+Z(17,J)*S+Z(18,J)*S*S
   10 F(J)=A*X**B*(1.-X)**C*(1.+D*SQRT(X)+E*X)*( LOG(1./X))**H
      UPV=F(1)
      DNV=F(2)
      GLU=F(3)
      DUSUM=F(4)
      DUDIF=F(5)
      STR=F(6)
      CHM=F(7)
      BOT=F(8)
      DSEA=0.5D0*(DUSUM+DUDIF)
      USEA=0.5D0*(DUSUM-DUDIF)
      IF(Q2.LT.2.7D0) CHM=0.D0
      IF(Q2.LT.30.D0) BOT=0.D0
      IF(CHM.LT.0.D0) CHM=0.D0
      IF(BOT.LT.0.D0) BOT=0.D0
      RETURN
      END
