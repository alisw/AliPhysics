#include "pdf/pilot.h"
      SUBROUTINE SFMFR2(X,Q2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
c  This program returns the fitted parton distributions for MRS(R2)  c
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
     .1.840938, -0.613912,  0.069060,  0.664453, -0.047312,  0.008339,
     .3.199642,  0.299148, -0.005252, -1.176342,  1.452607, -0.242739,
     .7.607104, -3.009800,  0.110586,  0.305600,  0.324526,  0.004361,

     .0.116204, -0.036926,  0.005132,  0.347559,  0.044734,  0.025428,
     .3.776571,  0.129517, -0.106180,  4.583171,  6.631731, -2.496486,
     .29.786677,-19.577891,  3.120062, 0.410322,  0.467855,  0.101884,

     .9.333184, -8.070131,  2.146747, -0.089195, -0.436328,  0.122834,
     .6.948321, -0.174262,  0.921887, -4.043493,  2.187837, -1.088141,
     .7.196254, -4.912953,  2.464480, -1.359569,  0.705893, -0.012063,

     .0.083723,  0.102184, -0.032540, -0.069562, -0.107065, -0.042555,
     .6.848496,  2.149966, -0.330963,  2.365100, -5.254575,  1.846434,
     .11.934685, -0.271707, -1.758501, 0.554814,  0.050133, -0.160565,

     .-0.221E-04, 0.120E-04,-0.234E-05,-0.603E+00, 0.251E+00,-0.146E-01,
     . 0.108E+02,-0.722E+00, 0.876E-02,-0.774E+04,-0.783E+03, 0.152E+04,
     .-0.125E+06, 0.519E+05,-0.507E+04,-0.244E+01, 0.146E+01, 0.237E-01,

     .0.285E-01, 0.481E-01,-0.128E-01,-0.686E-01,-0.139E+00,-0.255E-01,
     .0.694E+01, 0.204E+01,-0.290E+00, 0.105E+01,-0.445E+01, 0.154E+01,
     .0.988E+01, 0.519E+00,-0.159E+01, 0.447E+00, 0.123E+00,-0.158E+00,

     .0.536E-02,-0.149E-01, 0.337E-01,-0.136E+00, 0.397E-01,-0.114E+00,
     .0.459E+01, 0.333E+01,-0.115E+00,-0.309E+01, 0.189E+00,-0.109E+00,
     .0.182E+02,-0.122E+02, 0.286E+01,-0.323E+00, 0.207E+01,-0.108E+01,

     .-0.111E-01, 0.134E-01,-0.146E-02,-0.130E+00,-0.157E-01,-0.418E-01,
     . 0.276E+01, 0.277E+01,-0.697E-01, 0.124E+01,-0.255E+01, 0.466E+00,
     .-0.125E+01, 0.581E+01,-0.154E+01, 0.148E+01, 0.353E+00,-0.196E+00/
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
      AL=0.344D0
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
