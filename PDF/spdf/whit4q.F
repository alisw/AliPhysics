*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:04  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:16  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
c
cccccccccccccccccccccccccccccc
c   QPM calculation
      subroutine WHIT4Q(x,r,cv,cs)
ccc INPUTS : x,r=mc^2/Q^2
ccc OUTPUTS:   cv,cs (valence- and sea- charm quark dist)
ccc                  cv <-- cv / ( alpha / 2PI)
ccc                  cs <-- cs / ( alpha_s/2PI)
c
      implicit none
c arg
      DOUBLE PRECISION x,r
      DOUBLE PRECISION cv,cs
c CONST
      DOUBLE PRECISION ec,mc
      parameter(ec=2.0d0/3.0d0,mc=1.5d0)
c  N=15 Gauss int. weights and points
      integer GN,i
      parameter(GN=15)
      DOUBLE PRECISION XG(GN), XW(GN)
      DATA (XG(i),i=1,GN)/6.003741d-03,
     $ 3.136330d-02, 7.589671d-02, 1.377911d-01, 2.145139d-01,
     $ 3.029243d-01, 3.994030d-01, 5.000000d-01, 6.005970d-01,
     $ 6.970757d-01, 7.854861d-01, 8.622089d-01, 9.241033d-01,
     $ 9.686367d-01, 9.939963d-01/
      DATA (XW(i),i=1,GN)/1.537662d-02,
     $ 3.518302d-02, 5.357961d-02, 6.978534d-02, 8.313460d-02,
     $ 9.308050d-02, 9.921574d-02, 1.012891d-01, 9.921574d-02,
     $ 9.308050d-02, 8.313460d-02, 6.978534d-02, 5.357961d-02,
     $ 3.518302d-02, 1.537662d-02/
c local
      DOUBLE PRECISION sum, y,z,beta,w,WHIT4G,L
      parameter(L=4.0d0)
      DOUBLE PRECISION x1,rx,z1,rz
c
c begin
      x1=1.0d0-x
      rx=4.0d0*r*x
c
c direct
         beta=dsqrt(1.0d0-rx/x1)
         w=x*( beta*(-1.0d0+8.0d0*x*x1-rx*x1)
     $        +(x**2+x1**2+rx*(1.0d0-3.0d0*x)-0.5d0*rx**2)
     $         *  log( (1.0d0+beta)/(1.0d0-beta) ))
         cv = 3.0d0*ec**2 * w / x
c
c resolved
      sum=0.0d0
      do 10 i=1,GN
         y= x+rx + (x1-rx)*XG(i)**L
         z=x/y
         z1=1.0d0-z
         rz=4.0d0*r*z
         beta=dsqrt(1.0d0-rz/z1)
         w=z*( beta*(-1.0d0+8.0d0*z*z1-rz*z1)
     $        +(z**2+z1**2+rz*(1.0d0-3.0d0*z)-0.5d0*rz**2)
     $         *  log( (1.0d0+beta)/(1.0d0-beta) ))
c
         sum= sum + w * WHIT4G(y,mc**2/r)* L*XG(i)**(L-1.0d0)*XW(i)
c
   10 continue
c
         cs = 0.5d0/x * (x1-rx) * sum
c
      return
      end
