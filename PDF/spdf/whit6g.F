*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:06  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:16  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION function WHIT6G(x,Q2)
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
cccccc   WHIT6 gluon (U100)
c
      A0g = 1.600000d+01+s*(-6.100000d+01)+s2*( 1.278000d+02)
     $          +s3*(-1.399000d+02)+s4*( 5.990000d+01)
      B0g =              s*(-1.109000d+00)+s2*( 8.450000d-01)
     $          +s3*(-3.510000d-01)
      C0g = 1.500000d+01+s*( 1.596000d-01)+s2*( 4.180000d+00)
     $          +s3*(-1.765000d-01)
      A1g =              s*( 6.020000d-01)+s2*( 5.090000d-01)
     $          +s3*(-2.054000d+00)+s4*( 1.392000d+00)
      AA1g=              s*(-9.220000d-02)+s2*(-1.899000d+00)
     $          +s3*( 4.180000d+00)+s4*(-2.494000d+00)
      B1g =-2.895000d-01+s*( 3.760000d-01)+s2*(-1.719000d+00)
     $          +s3*( 1.116000d+00)
      C1g = 1.439000d+00+s*(-5.570000d-01)+s2*( 3.660000d-01)
     $          +s3*( 7.330000d-01)+s4*(-7.620000d-01)
      else
c over 100 GeV^2
c
cccccc   WHIT6 gluon (O100)
c
      A0g = 2.378000d+00+s*(-4.380000d+00)+s2*( 5.850000d-01)
     $          +s3*( 8.340000d+00)+s4*(-9.920000d+00)
      B0g =-4.790000d-01+s*(-6.070000d-01)+s2*( 1.458000d+00)
     $          +s3*(-6.030000d+00)+s4*( 9.330000d+00)
      C0g = 1.706000d+01+s*( 4.960000d+00)+s2*( 2.497000d+01)
     $          +s3*(-1.582000d+02)+s4*( 2.954000d+02)
      A1g = 2.992000d-01+s*( 1.179000d+00)+s2*(-1.915000d+00)
     $          +s3*( 7.260000d+00)+s4*(-1.839000d+01)
      AA1g=-1.600000d-01+s*(-1.114000d+00)+s2*( 2.939000d+00)
     $          +s3*(-6.660000d+00)+s4*( 1.923000d+01)
      B1g =-4.830000d-01+s*( 7.550000d-01)+s2*(-3.800000d+00)
     $          +s3*( 1.075000d+01)+s4*(-1.993000d+01)
      C1g = 1.297000d+00+s*(-1.669000d-01)+s2*( 1.906000d+00)
     $          +s3*(-2.057000d+00)
      endif
c
         WHIT6G = prsccf/alinv/x*
     $            ( A0g * x**B0g * x1**C0g
     $             +(A1g+AA1g*x) * x**B1g * x1**C1g )
c
      return
      end
