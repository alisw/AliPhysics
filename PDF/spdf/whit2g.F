*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:00  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:15  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION function WHIT2G(x,Q2)
c               input: x,Q2
c               output: clg
c                        (gluon dist.)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c arg
      DOUBLE PRECISION Q2,x
c const
      DOUBLE PRECISION q42it,q52it,lam42,lam52
      DOUBLE PRECISION alinv
c local
      DOUBLE PRECISION A0g,B0g,C0g,A1g,AA1g,B1g,C1g
      DOUBLE PRECISION s,s2,s3,s4,prsccf
      DOUBLE PRECISION x1
c parameters
      parameter(lam42=0.16d0, lam52=0.091411319d0)
      parameter(Q42IT=4.0d0, Q52IT=100.0d0)
      parameter(alinv=137.036d0)
      common /scale/ s,s2,s3,s4,prsccf
c
c begin
      x1=1.0d0-x
c
      if(Q2.lt.100.0d0) then
c  under 100 GeV^2
c
cccccc   WHIT2 gluon (U100)
c
      A0g = 5.000000d+00+s*(-1.499000d+01)+s2*( 2.617000d+01)
     $          +s3*(-2.530000d+01)+s4*( 1.012000d+01)
      B0g =              s*(-9.370000d-01)+s2*( 4.100000d-01)
     $          +s3*( 3.390000d-02)
      C0g = 9.000000d+00+s*( 7.090000d-01)+s2*( 3.118000d+00)
     $          +s3*(-5.820000d-04)
      A1g =              s*( 4.610000d-01)+s2*( 1.041000d-01)
     $          +s3*(-1.753000d-02)+s4*(-2.717000d-01)
      AA1g=              s*( 9.680000d-03)+s2*(-4.170000d-01)
     $          +s3*(-3.950000d-01)+s4*( 8.430000d-01)
      B1g =-4.140000d-01+s*(-6.060000d-02)+s2*( 2.847000d-01)
     $          +s3*(-5.070000d-01)
      C1g = 1.244000d+00+s*( 5.880000d-01)+s2*(-1.228000d+00)
     $          +s3*( 8.090000d-01)
      else
c over 100 GeV^2
c
cccccc   WHIT2 gluon (O100)
c
      A0g = 1.095000d+00+s*(-2.388000d+00)+s2*( 9.190000d+00)
     $          +s3*(-3.032000d+01)+s4*( 3.480000d+01)
      B0g =-4.410000d-01+s*(-9.070000d-01)+s2*( 4.680000d+00)
     $          +s3*(-1.866000d+01)+s4*( 2.717000d+01)
      C0g = 1.099000d+01+s*( 4.710000d+00)+s2*( 2.801000d+01)
     $          +s3*(-1.279000d+02)+s4*( 1.640000d+02)
      A1g = 3.010000d-01+s*( 1.275000d+00)+s2*(-1.563000d+00)
     $          +s3*( 4.100000d+00)+s4*(-1.337000d+01)
      AA1g=-1.305000d-01+s*(-1.245000d+00)+s2*( 2.438000d+00)
     $          +s3*(-2.539000d+00)+s4*( 1.273000d+01)
      B1g =-4.890000d-01+s*( 9.550000d-01)+s2*(-4.400000d+00)
     $          +s3*( 1.022000d+01)+s4*(-1.713000d+01)
      C1g = 1.331000d+00+s*(-2.481000d-01)+s2*( 1.950000d+00)
     $          +s3*(-2.072000d+00)
      endif
c
         WHIT2G = prsccf/alinv/x*
     $            ( A0g * x**B0g * x1**C0g
     $             +(A1g+AA1g*x) * x**B1g * x1**C1g )
c
      return
      end
