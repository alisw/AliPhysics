#include "pdf/pilot.h"
      SUBROUTINE SFMFR1(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
c  This program returns the fitted parton distributions for MRS(R1)  c
c                                                                    c
c  Note::  The gluon fit below Q2=2 GeV2 is not reliable             c
c                                                                    c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
#include "pdf/impdp.inc"
C
#include "pdf/expdp.inc"
     +       F(8)
      REAL   Z(18,8)
      DATA Z/
     .1.851372, -0.519807,  0.012441,  0.633560, -0.037308,  0.015912,
     .3.375911,  0.342401, -0.051574, -1.120042,  1.261935, -0.260027,
     .7.269788, -3.120148,  0.333912,  0.267482,  0.294300,  0.057360,

     .0.118684, -0.026694,  0.001418,  0.309073,  0.036590,  0.020734,
     .4.069154,  0.290830, -0.125778,  4.041269,  7.626688, -2.807994,
     .29.897693,-22.979920,  4.180347, 0.334020,  0.330986,  0.114984,

     .9.053009, -8.482741,  2.830186,  0.053474, -0.589212,  0.144551,
     .5.812161,  2.209058,  0.424067, -4.014258,  0.658552, -0.367525,
     .5.807276,  0.224226,  0.154323, -0.520160,  0.182673, -0.104011,

     .0.120167,  0.167327, -0.053467, -0.054880, -0.148988, -0.048956,
     .7.759331,  2.553270, -0.438661,  1.178532, -5.710315,  2.380271,
     .11.895434,  0.861384, -2.528848, 0.493662, -0.073144, -0.207566,

     .0.177E-03,-0.159E-03, 0.953E-04,-0.863E+00,-0.352E+00, 0.338E+00,
     .0.126E+02, 0.887E+00,-0.929E+00, 0.394E+03,-0.978E+03, 0.542E+03,
     .0.236E+05, 0.358E+04,-0.137E+05,-0.326E+01,-0.248E-01, 0.973E+00,

     .0.049569,  0.044433, -0.016279, -0.105288, -0.087260, -0.050832,
     .8.515981, -0.203986,  0.963570,  0.588399, -4.727176,  1.568474,
     .9.610173, -1.946171,  0.998092,  0.108582,  0.643896, -0.358571,

     .0.001324,  0.020995,  0.033861, -0.116455, -0.101398, -0.054482,
     .7.220589,  1.135903,  0.961005, -5.142936,  2.318274, -1.057351,
     .14.319705,-10.025316,  3.619585, -0.692507, 2.084232, -1.105329,

     .-0.011228,  0.014937, -0.000135, -0.109754, -0.019479, -0.064273,
     . 2.973424,  3.819107, -0.168367, -0.662015, -2.241463,  0.521234,
     . 0.936045,  6.588813, -2.031513,  1.507347,  0.341242, -0.313753/
      DATA ZEROD/0.D0/,PONED/0.1D0/,ONED/1.D0/,ONEDO/1.1D0/,TWOD/2.D0/
C
      UPV= ZEROD
      DNV= ZEROD
      USEA=ZEROD
      DSEA=ZEROD
      STR= ZEROD
      CHM= ZEROD
      BOT= ZEROD
      GLU= ZEROD
      IF(X.GE.ONED) RETURN
C
      AL=0.241D0
      AL2=AL*AL
      Q02=ONED
      ARG=( LOG(Q2/AL2)/ LOG(Q02/AL2))
      S= LOG(ARG)
C
      DO 10 J=1,8
C
      A=Z(1,J)+Z(2,J)*S+Z(3,J)*S*S
      B=Z(4,J)+Z(5,J)*S+Z(6,J)*S*S
      C=Z(7,J)+Z(8,J)*S+Z(9,J)*S*S
      D=Z(10,J)+Z(11,J)*S+Z(12,J)*S*S
      E=Z(13,J)+Z(14,J)*S+Z(15,J)*S*S
      H=Z(16,J)+Z(17,J)*S+Z(18,J)*S*S 
   10 F(J)=A*X**B*(ONED-X)**C*(ONED+D*SQRT(X)+E*X)*( LOG(ONED/X))**H
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
      IF(Q2.LT.2.7D0) CHM=ZEROD
      IF(Q2.LT.30.D0) BOT=ZEROD
      IF(CHM.LT.ZEROD) CHM=ZEROD
      IF(BOT.LT.ZEROD) BOT=ZEROD
      RETURN
      END
