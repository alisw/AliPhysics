C**********************************************************************
C*                                                                    *
C* Exact calculation of the total differential e+ e- -pair production *
C* in Relativistic Heavy Ion Collisions for a point particle in an    *
C* external field approach.                                           *
C*                                                                    *
C* For details see the publication:                                   *
C* "Multiple electromagnetic electron positron pair production in     *
C*  relativistic heavy ion collisions"                                *
C*  Adrian Alscher, Kai Hencken, Dirk Trautmann, and Gerhard Baur     *
C*  Phys. Rev. A55 (1997) 396.                                        *
C*                                                                    *
C* Copyright (c) 1997, 1998, 2002, Adrian Alscher and Kai Hencken     *
C*                                                                    *
C* Permission to use, copy and distribute this software and its       *
C* documentation strictly for non-commercial purposes is hereby       *
C* granted without fee, provided that the above copyright notice      *
C* appears in all copies, that both the copyright notice and this     *
C* permission notice appear in the supporting documentation and that  *
C* the use of this program is acknowledged in scientific              *
C* publications (see reference above). The authors make no claims     *
C* about the suitability of this software for any purpose. It is      *
C* provided "as is" without express or implied warranty. Any change   *
C* of the code should be submitted to the authors                     *
C*                                                                    *
C* To use this program at LHC energies, please make sure that         *
C* "double precision" variables should better be real*16              *
C**********************************************************************

C======================================================================
C
C call this routine first to initialize some parameter needed in the
C function. 
C
C gm = Gamma_cm, that is, gm of each ion (~100 for RHIC ~3000 for LHC)
C mass = mass of the produced particle in MeV (~0.511 for e,~100 for mu)
C======================================================================

      SUBROUTINE Initdiffcross (gm,mass)
      IMPLICIT NONE
      DOUBLE PRECISION gm,mass

      DOUBLE PRECISION gamma,beta,m,w1xw1,w1xw2,w2xw2,wl,wy
      COMMON/PHYSPARAM/gamma,beta,m,w1xw1,w1xw2,w2xw2,wl,wy

      DOUBLE PRECISION ARCOSH,x
      ARCOSH(x)=LOG(x+SQRT((x-1D0)*(x+1D0)))

      gamma=gm
      beta=sqrt((1D0-1D0/gamma)*(1D0+1D0/gamma))
      m=mass
      
      w1xw1 = 1D0/gamma**2
      w2xw2 = 1D0/gamma**2
      w1xw2 = 2D0-1D0/gamma**2
      wl=SQRT (w1xw1)
      wy=ARCOSH (gamma)
      
      RETURN 
      END

C ============================================================================
C
C Diffcross calculates the fivefold differential cross section
C
C dsigma/ dp+t dp-t dy+ dy- ddeltaphi
C
C the trivial integration over the total phi-dependence is not included
C
C ppvt= absolute value of the positron transverse momentum (MeV)
C pmvt=  - " -                electron - " -
C dphi= phi-angle between the electron and positron transverse momentum
C yp = rapidity of the positron
C ym = rapidity of the electron
C      defined to make E = sqrt(pt^2 + m^2) cosh (y)
C      and             Pz= sqrt(pt^2 + m^2) sinh (y)
C
C dsigma = differential cross section in kbarn/MeV^4/ (Zalpha)**4
C
C to get the total cross section, you have to integrate over
C
C Integral dsigma dyp dym ddphi 2 pi ppt dppt pmt dpmt
C
C (see also source for sigma.f)
C
C======================================================================

      SUBROUTINE Diffcross(ppvt,yp,pmvt,ym,dphi,dsigma)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION ppvt,pmvt,dphi,dsigma,yp,ym
      DOUBLE PRECISION ppvl,pmvl
      
      DOUBLE PRECISION k1(2)
      DOUBLE PRECISION kd(2)
      DOUBLE PRECISION kx(2)
      DOUBLE PRECISION pmt(2)
      DOUBLE PRECISION ppt(2)
      DOUBLE PRECISION qb
      DOUBLE PRECISION mk1(2),mk1d(2),mk1x(2)
      DOUBLE PRECISION mkd(2),mkd1(2)
      DOUBLE PRECISION mkx(2),mkx1(2)
      
      DOUBLE PRECISION Iz2,Id1,Id2,Id3,Iv1,Iv2
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793238462643D0)
      
      DOUBLE PRECISION gamma,beta,m,w1xw1,w1xw2,w2xw2,wl,wy
      COMMON/PHYSPARAM/gamma,beta,m,w1xw1,w1xw2,w2xw2,wl,wy
      
      DOUBLE PRECISION m0,m1,md,mx
      DOUBLE PRECISION N1,N2,N3,N4,N5
      DOUBLE PRECISION N6,N7,N8,N9,N10
      DOUBLE PRECISION N11,N12,N13,N14,N15
      DOUBLE PRECISION N16,N17,N18,NT
      
      DOUBLE PRECISION pmlxpml,pmlxppl,pmlxql
      DOUBLE PRECISION pmtxpmt,pmtxppt,pptxppt
      DOUBLE PRECISION pplxppl,pplxql,qlxql
      DOUBLE PRECISION w1xpml,w1xppl,w1xql
      DOUBLE PRECISION w2xpml,w2xppl,w2xql
      
      DOUBLE PRECISION r1,rd,rx

      INTEGER setzero
      COMMON /SETZPARAM/ setzero
      INTEGER badcount
      COMMON/badpar/badcount

      setzero=0
      
      ppt(1)=ppvt
      ppt(2)=0D0
      ppvl=SQRT(ppvt**2+m**2)
      
      pmt(1)=pmvt*cos(dphi)
      pmt(2)=pmvt*sin(dphi)
      pmvl=SQRT(pmvt**2+m**2)
      
      pmlxpml = pmvl**2
      pplxppl = ppvl**2
      pmlxppl = pmvl*ppvl*COSH(yp-ym)
      w2xpml = wl*pmvl*COSH(ym+wy)
      w2xppl = wl*ppvl*COSH(yp+wy)
      w1xpml = wl*pmvl*COSH(ym-wy)
      w1xppl = wl*ppvl*COSH(yp-wy)
      w1xql=0D0
      w2xql=w2xpml+w2xppl
      qb=gamma*(w2xppl+w2xpml)/SINH(2D0*wy)
      qlxql=-qb**2
      pmlxql = qb*pmvl*SINH(wy-ym)
      pplxql = qb*ppvl*SINH(wy-yp)
      
      pmtxpmt =-pmvt**2
      pmtxppt =-pmvt*ppvt*cos(dphi)
      pptxppt =-ppvt**2
      
C======================================================================
C
C Definition of propagatorterms
C
      m0 = -qlxql
      
      k1(1) = -ppt(1) -pmt(1)
      k1(2) = -pmt(2)
      m1=-qlxql-pplxppl-pmlxpml+2D0*(pplxql+pmlxql-pmlxppl)
      
      kd(1)= -pmt(1)
      kd(2)= -pmt(2)
      md= m**2 - (qlxql+pmlxpml-2D0*pmlxql)
      
      kx(1)= -ppt(1)	
      kx(2)=0D0
      mx= m**2 - (qlxql+pplxppl-2D0*pplxql)
      
      r1 = m1 - m0 + k1(1)*k1(1)+k1(2)*k1(2)
      rd = md - m0 + kd(1)*kd(1)+kd(2)*kd(2)
      rx = mx - m0 + kx(1)*kx(1)+kx(2)*kx(2)
      
      mk1(1)= -k1(1)
      mk1(2)= -k1(2)
      mkd(1)= -kd(1)
      mkd(2)= -kd(2)
      mkx(1)= -kx(1)
      mkx(2)= -kx(2)
      mk1d(1)= k1(1)-kd(1)
      mk1d(2)= k1(2)-kd(2)
      mk1x(1)= k1(1)-kx(1)
      mk1x(2)= k1(2)-kx(2)
      mkd1(1)= kd(1)-k1(1)
      mkd1(2)= kd(2)-k1(2)
      mkx1(1)= kx(1)-k1(1)
      mkx1(2)= kx(2)-k1(2)

C
C calculate all different integrals
C

      N1 = Id1(mkd,mk1d,md,m0,m1) * (  - 2*w1xw1*w2xw2 )
      
      N2 = Id1(mkx,mk1x,mx,m0,m1) * (  - 2*w1xw1*w2xw2 )
      
      N3 = Iv1(mk1,mkd1,mkx1,m1,m0,md,mx) * ( 4*w1xw1*w2xw2*
     &  pmlxpml + 8*w1xw1*w2xw2*pmlxppl - 8*w1xw1*w2xw2*pmlxql + 4*
     &  w1xw1*w2xw2*pmtxpmt + 8*w1xw1*w2xw2*pmtxppt + 4*w1xw1*w2xw2*
     &  pplxppl - 8*w1xw1*w2xw2*pplxql + 4*w1xw1*w2xw2*pptxppt + 4*
     &  w1xw1*w2xw2*rd + 4*w1xw1*w2xw2*rx - 16*w1xw1*w2xpml*w2xppl + 
     &  8*w1xw1*w2xpml*w2xql + 8*w1xw1*w2xppl*w2xql - 8*w1xw2**2*
     &  pmlxpml - 16*w1xw2**2*pmlxppl + 16*w1xw2**2*pmlxql - 8*
     &  w1xw2**2*pmtxpmt - 16*w1xw2**2*pmtxppt - 8*w1xw2**2*pplxppl
     &  + 16*w1xw2**2*pplxql - 8*w1xw2**2*pptxppt - 8*w1xw2**2*rd - 
     &  8*w1xw2**2*rx + 16*w1xw2*w1xpml*w2xpml + 16*w1xw2*w1xpml*
     &  w2xppl - 16*w1xw2*w1xpml*w2xql + 16*w1xw2*w1xppl*w2xpml + 16*
     &  w1xw2*w1xppl*w2xppl - 16*w1xw2*w1xppl*w2xql - 16*w1xw2*w1xql*
     &  w2xpml - 16*w1xw2*w1xql*w2xppl - 8*w1xpml**2*w2xw2 + 8*w1xpml
     &  *w1xql*w2xw2 - 8*w1xppl**2*w2xw2 + 8*w1xppl*w1xql*w2xw2 )
      
      N4 = Id1(mk1,mkd1,m1,m0,md) * (  - 2*w1xw1*w2xw2 + 8*
     &  w1xw2**2 )
      
      N5 = Id2(mk1d,mkd,md,m1,m0) * ( 4*w1xw1*w2xw2*pmlxppl + 4*
     &  w1xw1*w2xw2*pmtxppt - 4*w1xw1*w2xw2*pplxql + 4*w1xw1*w2xw2*
     &  m**2 + 2*w1xw1*w2xw2*r1 - 2*w1xw1*w2xw2*rd - 8*w1xw1*w2xpml*
     &  w2xppl + 8*w1xw1*w2xppl*w2xql )
      
      N6 = Id1(mk1,mkx1,m1,m0,mx) * (  - 2*w1xw1*w2xw2 + 8*
     &  w1xw2**2 )
      
      N7 = Id2(mk1x,mkx,mx,m1,m0) * ( 4*w1xw1*w2xw2*pmlxppl - 4*
     &  w1xw1*w2xw2*pmlxql + 4*w1xw1*w2xw2*pmtxppt + 4*w1xw1*w2xw2*
     &  m**2 + 2*w1xw1*w2xw2*r1 - 2*w1xw1*w2xw2*rx - 8*w1xw1*w2xpml*
     &  w2xppl + 8*w1xw1*w2xpml*w2xql )
      
      N8 = Id1(k1,kd,m0,m1,md) * ( 2*w1xw1*w2xw2 )
      
      N9 = Id2(kd,k1,m0,md,m1) * (  - 4*w1xw1*w2xw2*pmlxpml
     &  + 4*w1xw1*w2xw2*pmlxql - 4*w1xw1*w2xw2*pmtxpmt + 4*w1xw1*
     &  w2xw2*m**2 - 2*w1xw1*w2xw2*rd + 8*w1xpml**2*w2xw2 - 8*w1xpml*
     &  w1xql*w2xw2 )
      
      N10 = Id1(k1,kx,m0,m1,mx) * ( 2*w1xw1*w2xw2 )
      
      N11 = Id2(kx,k1,m0,mx,m1) * (  - 4*w1xw1*w2xw2*pplxppl
     &  + 4*w1xw1*w2xw2*pplxql - 4*w1xw1*w2xw2*pptxppt + 4*w1xw1*
     &  w2xw2*m**2 - 2*w1xw1*w2xw2*rx + 8*w1xppl**2*w2xw2 - 8*w1xppl*
     &  w1xql*w2xw2 )
      
      N12 = Iv2(k1,kd,kx,m0,m1,md,mx) * ( 8*w1xw1*w2xw2*
     &  pmlxpml*pplxppl - 8*w1xw1*w2xw2*pmlxpml*pplxql + 8*w1xw1*
     &  w2xw2*pmlxpml*pptxppt - 8*w1xw1*w2xw2*pmlxpml*m**2 + 4*w1xw1*
     &  w2xw2*pmlxpml*rx - 8*w1xw1*w2xw2*pmlxppl*qlxql - 8*w1xw1*
     &  w2xw2*pmlxppl*m0 - 8*w1xw1*w2xw2*pmlxql*pplxppl + 16*w1xw1
     &  *w2xw2*pmlxql*pplxql - 8*w1xw1*w2xw2*pmlxql*pptxppt + 8*w1xw1
     &  *w2xw2*pmlxql*m**2 - 8*w1xw1*w2xw2*pmlxql*rx + 8*w1xw1*w2xw2*
     &  pmtxpmt*pplxppl - 8*w1xw1*w2xw2*pmtxpmt*pplxql + 8*w1xw1*
     &  w2xw2*pmtxpmt*pptxppt - 8*w1xw1*w2xw2*pmtxpmt*m**2 + 4*w1xw1*
     &  w2xw2*pmtxpmt*rx - 8*w1xw1*w2xw2*pmtxppt*qlxql - 8*w1xw1*
     &  w2xw2*pmtxppt*m0 - 8*w1xw1*w2xw2*pplxppl*m**2 + 4*w1xw1*
     &  w2xw2*pplxppl*rd + 8*w1xw1*w2xw2*pplxql*m**2 - 8*w1xw1*w2xw2*
     &  pplxql*rd - 8*w1xw1*w2xw2*pptxppt*m**2 + 4*w1xw1*w2xw2*
     &  pptxppt*rd - 8*w1xw1*w2xw2*qlxql*m**2 + 8*w1xw1*w2xw2*m**4 - 
     &  8*w1xw1*w2xw2*m**2*m0 - 4*w1xw1*w2xw2*m**2*rd - 4*w1xw1*
     &  w2xw2*m**2*rx + 4*w1xw1*w2xw2*rd*rx + 16*w1xw1*w2xpml*w2xppl*
     &  qlxql +  16*w1xw1*w2xpml*
     &  w2xppl*m0 - 16*w1xw1*w2xpml*w2xql*pplxql + 8*w1xw1*w2xpml*
     &  w2xql*rx - 16*w1xw1*w2xppl*w2xql*pmlxql + 8*w1xw1*w2xppl*
     &  w2xql*rd + 16*w1xw1*w2xql**2*pmlxppl + 16*w1xw1*w2xql**2*
     &  pmtxppt + 16*w1xw1*w2xql**2*m**2 - 16*w1xw2**2*pmlxpml*
     &  pplxppl + 16*w1xw2**2*pmlxpml*pplxql - 16*w1xw2**2*pmlxpml*
     &  pptxppt + 16*w1xw2**2*pmlxpml*m**2 - 8*w1xw2**2*pmlxpml*rx + 
     &  16*w1xw2**2*pmlxppl*qlxql + 16*w1xw2**2*pmlxppl*m0 + 16*
     &  w1xw2**2*pmlxql*pplxppl - 32*w1xw2**2*pmlxql*pplxql + 16*
     &  w1xw2**2*pmlxql*pptxppt - 16*w1xw2**2*pmlxql*m**2 + 16*
     &  w1xw2**2*pmlxql*rx - 16*w1xw2**2*pmtxpmt*pplxppl + 16*
     &  w1xw2**2*pmtxpmt*pplxql - 16*w1xw2**2*pmtxpmt*pptxppt + 16*
     &  w1xw2**2*pmtxpmt*m**2 - 8*w1xw2**2*pmtxpmt*rx + 16*w1xw2**2*
     &  pmtxppt*qlxql + 16*w1xw2**2*pmtxppt*m0 + 16*w1xw2**2*
     &  pplxppl*m**2 - 8*w1xw2**2*pplxppl*rd - 16*w1xw2**2*pplxql*
     &  m**2 + 16*w1xw2**2*pplxql*rd + 16*w1xw2**2*pptxppt*m**2 - 8*
     &  w1xw2**2*pptxppt*rd + 16*w1xw2**2*qlxql
     &  *m**2 - 16*w1xw2**2*m**4 + 16*w1xw2**2*m**2*m0 + 8*
     &  w1xw2**2*m**2*rd + 8*w1xw2**2*m**2*rx - 8*w1xw2**2*rd*rx + 32
     &  *w1xw2*w1xpml*w2xpml*pplxppl - 32*w1xw2*w1xpml*w2xpml*pplxql
     &  + 32*w1xw2*w1xpml*w2xpml*pptxppt - 32*w1xw2*w1xpml*w2xpml*
     &  m**2 + 16*w1xw2*w1xpml*w2xpml*rx - 16*w1xw2*w1xpml*w2xppl*
     &  qlxql - 16*w1xw2*w1xpml*w2xppl*m0 - 16*w1xw2*w1xpml*w2xql*
     &  pplxppl + 32*w1xw2*w1xpml*w2xql*pplxql - 16*w1xw2*w1xpml*
     &  w2xql*pptxppt + 16*w1xw2*w1xpml*w2xql*m**2 - 16*w1xw2*w1xpml*
     &  w2xql*rx - 16*w1xw2*w1xppl*w2xpml*qlxql - 16*w1xw2*w1xppl*
     &  w2xpml*m0 + 32*w1xw2*w1xppl*w2xppl*pmlxpml - 32*w1xw2*
     &  w1xppl*w2xppl*pmlxql + 32*w1xw2*w1xppl*w2xppl*pmtxpmt - 32*
     &  w1xw2*w1xppl*w2xppl*m**2 + 16*w1xw2*w1xppl*w2xppl*rd - 16*
     &  w1xw2*w1xppl*w2xql*pmlxpml + 32*w1xw2*w1xppl*w2xql*pmlxql - 
     &  16*w1xw2*w1xppl*w2xql*pmtxpmt + 16*w1xw2*w1xppl*w2xql*m**2 - 
     &  16*w1xw2*w1xppl*w2xql*rd - 16*w1xw2*w1xql*w2xpml*pplxppl + 32
     &  *w1xw2*w1xql*w2xpml*pplxql - 16*w1xw2*w1xql
     &  *w2xpml*pptxppt + 16*w1xw2*w1xql*w2xpml*m**2 - 16*w1xw2*w1xql
     &  *w2xpml*rx - 16*w1xw2*w1xql*w2xppl*pmlxpml + 32*w1xw2*w1xql*
     &  w2xppl*pmlxql - 16*w1xw2*w1xql*w2xppl*pmtxpmt + 16*w1xw2*
     &  w1xql*w2xppl*m**2 - 16*w1xw2*w1xql*w2xppl*rd - 32*w1xw2*w1xql
     &  *w2xql*pmlxppl - 32*w1xw2*w1xql*w2xql*pmtxppt - 32*w1xw2*
     &  w1xql*w2xql*m**2 - 16*w1xpml**2*w2xw2*pplxppl + 16*w1xpml**2*
     &  w2xw2*pplxql - 16*w1xpml**2*w2xw2*pptxppt + 16*w1xpml**2*
     &  w2xw2*m**2 - 8*w1xpml**2*w2xw2*rx + 32*w1xpml*w1xppl*w2xw2*
     &  pmlxppl - 16*w1xpml*w1xppl*w2xw2*pmlxql + 32*w1xpml*w1xppl*
     &  w2xw2*pmtxppt - 16*w1xpml*w1xppl*w2xw2*pplxql + 16*w1xpml*
     &  w1xppl*w2xw2*qlxql + 32*w1xpml*w1xppl*w2xw2*m**2 + 16*w1xpml*
     &  w1xppl*w2xw2*m0 + 8*w1xpml*w1xppl*w2xw2*rd + 8*w1xpml*
     &  w1xppl*w2xw2*rx - 64*w1xpml*w1xppl*w2xpml*w2xppl + 32*w1xpml*
     &  w1xppl*w2xpml*w2xql + 32*w1xpml*w1xppl*w2xppl*w2xql - 32*
     &  w1xpml*w1xppl*w2xql**2 - 16*w1xpml*w1xql*w2xw2*pmlxppl - 16*
     &  w1xpml*w1xql*w2xw2*pmtxppt + 16*w1xpml*w1xql*
     &  w2xw2*pplxppl - 16*w1xpml*w1xql*w2xw2*pplxql + 16*w1xpml*
     &  w1xql*w2xw2*pptxppt - 32*w1xpml*w1xql*w2xw2*m**2 + 8*w1xpml*
     &  w1xql*w2xw2*rx + 32*w1xpml*w1xql*w2xpml*w2xppl - 16*w1xppl**2
     &  *w2xw2*pmlxpml + 16*w1xppl**2*w2xw2*pmlxql - 16*w1xppl**2*
     &  w2xw2*pmtxpmt + 16*w1xppl**2*w2xw2*m**2 - 8*w1xppl**2*w2xw2*
     &  rd + 16*w1xppl*w1xql*w2xw2*pmlxpml - 16*w1xppl*w1xql*w2xw2*
     &  pmlxppl - 16*w1xppl*w1xql*w2xw2*pmlxql + 16*w1xppl*w1xql*
     &  w2xw2*pmtxpmt - 16*w1xppl*w1xql*w2xw2*pmtxppt - 32*w1xppl*
     &  w1xql*w2xw2*m**2 + 8*w1xppl*w1xql*w2xw2*rd + 32*w1xppl*w1xql*
     &  w2xpml*w2xppl + 16*w1xql**2*w2xw2*pmlxppl + 16*w1xql**2*w2xw2
     &  *pmtxppt + 16*w1xql**2*w2xw2*m**2 - 32*w1xql**2*w2xpml*w2xppl
     &  )
      
      N13 = Id2(k1,kd,m0,m1,md) * ( 4*w1xw1*w2xw2*pmlxql + 4*
     &  w1xw1*w2xw2*pplxql - 2*w1xw1*w2xw2*r1 - 8*w1xw1*w2xpml*w2xql
     &  - 8*w1xw1*w2xppl*w2xql + 8*w1xw2**2*pmlxpml - 16*w1xw2**2*
     &  pmlxql + 8*w1xw2**2*pmtxpmt - 8*w1xw2**2*m**2 + 8*w1xw2**2*rd
     &  - 16*w1xw2*w1xpml*w2xpml + 16*w1xw2*w1xpml*w2xppl + 16*w1xw2
     &  *w1xpml*w2xql + 16*w1xw2*w1xql*w2xpml - 16*w1xpml*w1xppl*
     &  w2xw2 )
      
      N14 = Id3(k1,kd,m0,m1,md) * ( 4*w1xw1*w2xw2*pmlxpml*
     &  pmlxppl + 4*w1xw1*w2xw2*pmlxpml*pmtxppt - 8*w1xw1*w2xw2*
     &  pmlxpml*pplxql + 4*w1xw1*w2xw2*pmlxpml*m**2 + 4*w1xw1*w2xw2*
     &  pmlxpml*r1 - 4*w1xw1*w2xw2*pmlxpml*rd + 4*w1xw1*w2xw2*pmlxppl
     &  *pmtxpmt - 4*w1xw1*w2xw2*pmlxppl*qlxql - 4*w1xw1*w2xw2*
     &  pmlxppl*m**2 - 4*w1xw1*w2xw2*pmlxppl*m0 + 8*w1xw1*w2xw2*
     &  pmlxql*pplxql - 4*w1xw1*w2xw2*pmlxql*r1 + 4*w1xw1*w2xw2*
     &  pmlxql*rd + 4*w1xw1*w2xw2*pmtxpmt*pmtxppt - 8*w1xw1*w2xw2*
     &  pmtxpmt*pplxql + 4*w1xw1*w2xw2*pmtxpmt*m**2 + 4*w1xw1*w2xw2*
     &  pmtxpmt*r1 - 4*w1xw1*w2xw2*pmtxpmt*rd - 4*w1xw1*w2xw2*pmtxppt
     &  *qlxql - 4*w1xw1*w2xw2*pmtxppt*m**2 - 4*w1xw1*w2xw2*pmtxppt*
     &  m0 + 8*w1xw1*w2xw2*pplxql*m**2 - 4*w1xw1*w2xw2*pplxql*rd
     &  - 4*w1xw1*w2xw2*qlxql*m**2 - 4*w1xw1*w2xw2*m**4 - 4*w1xw1*
     &  w2xw2*m**2*m0 - 4*w1xw1*w2xw2*m**2*r1 + 4*w1xw1*w2xw2*m**2
     &  *rd + 2*w1xw1*w2xw2*r1*rd - 2*w1xw1*w2xw2*rd**2 - 8*w1xw1*
     &  w2xpml*w2xppl*pmlxpml - 8*w1xw1*w2xpml*w2xppl*pmtxpmt + 8*
     &  w1xw1*w2xpml*w2xppl*qlxql + 8*w1xw1*w2xpml*w2xppl*m**2
     &  + 8*w1xw1*w2xpml*w2xppl*m0 + 16*w1xw1*w2xppl*w2xql*
     &  pmlxpml - 16*w1xw1*w2xppl*w2xql*pmlxql + 16*w1xw1*w2xppl*
     &  w2xql*pmtxpmt - 16*w1xw1*w2xppl*w2xql*m**2 + 8*w1xw1*w2xppl*
     &  w2xql*rd - 16*w1xw2*w1xpml*w2xppl*pmlxpml + 32*w1xw2*w1xpml*
     &  w2xppl*pmlxql - 16*w1xw2*w1xpml*w2xppl*pmtxpmt - 16*w1xw2*
     &  w1xpml*w2xppl*qlxql + 16*w1xw2*w1xpml*w2xppl*m**2 - 16*w1xw2*
     &  w1xpml*w2xppl*m0 - 16*w1xw2*w1xpml*w2xppl*rd - 16*
     &  w1xpml**2*w2xw2*pmlxppl - 16*w1xpml**2*w2xw2*pmtxppt + 16*
     &  w1xpml**2*w2xw2*pplxql - 16*w1xpml**2*w2xw2*m**2 - 8*
     &  w1xpml**2*w2xw2*r1 + 8*w1xpml**2*w2xw2*rd + 32*w1xpml**2*
     &  w2xpml*w2xppl - 32*w1xpml**2*w2xppl*w2xql + 8*w1xpml*w1xppl*
     &  w2xw2*pmlxpml - 16*w1xpml*w1xppl*w2xw2*pmlxql + 8*w1xpml*
     &  w1xppl*w2xw2*pmtxpmt + 8*w1xpml*w1xppl*w2xw2*qlxql - 8*w1xpml
     &  *w1xppl*w2xw2*m**2 + 8*w1xpml*w1xppl*w2xw2*m0 + 8*w1xpml*
     &  w1xppl*w2xw2*rd + 16*w1xpml*w1xql*w2xw2*pmlxppl + 16*w1xpml*
     &  w1xql*w2xw2*pmtxppt - 16*w1xpml*w1xql*w2xw2*
     &  pplxql + 16*w1xpml*w1xql*w2xw2*m**2 + 8*w1xpml*w1xql*w2xw2*r1
     &  - 8*w1xpml*w1xql*w2xw2*rd - 32*w1xpml*w1xql*w2xpml*w2xppl + 
     &  32*w1xpml*w1xql*w2xppl*w2xql )
      
      N15 = Id2(k1,kx,m0,m1,mx) * ( 4*w1xw1*w2xw2*pmlxql + 4*
     &  w1xw1*w2xw2*pplxql - 2*w1xw1*w2xw2*r1 - 8*w1xw1*w2xpml*w2xql
     &  - 8*w1xw1*w2xppl*w2xql + 8*w1xw2**2*pplxppl - 16*w1xw2**2*
     &  pplxql + 8*w1xw2**2*pptxppt - 8*w1xw2**2*m**2 + 8*w1xw2**2*rx
     &  + 16*w1xw2*w1xppl*w2xpml - 16*w1xw2*w1xppl*w2xppl + 16*w1xw2
     &  *w1xppl*w2xql + 16*w1xw2*w1xql*w2xppl - 16*w1xpml*w1xppl*
     &  w2xw2 )
      
      N16 = Id3(k1,kx,m0,m1,mx) * ( 4*w1xw1*w2xw2*pmlxppl*
     &  pplxppl + 4*w1xw1*w2xw2*pmlxppl*pptxppt - 4*w1xw1*w2xw2*
     &  pmlxppl*qlxql - 4*w1xw1*w2xw2*pmlxppl*m**2 - 4*w1xw1*w2xw2*
     &  pmlxppl*m0 - 8*w1xw1*w2xw2*pmlxql*pplxppl + 8*w1xw1*w2xw2*
     &  pmlxql*pplxql - 8*w1xw1*w2xw2*pmlxql*pptxppt + 8*w1xw1*w2xw2*
     &  pmlxql*m**2 - 4*w1xw1*w2xw2*pmlxql*rx + 4*w1xw1*w2xw2*pmtxppt
     &  *pplxppl + 4*w1xw1*w2xw2*pmtxppt*pptxppt - 4*w1xw1*w2xw2*
     &  pmtxppt*qlxql - 4*w1xw1*w2xw2*pmtxppt*m**2 - 4*w1xw1*w2xw2*
     &  pmtxppt*m0 + 4*w1xw1*w2xw2*pplxppl*m**2 + 4*w1xw1*w2xw2*
     &  pplxppl*r1 - 4*w1xw1*w2xw2*pplxppl*rx - 4*w1xw1*w2xw2*pplxql*
     &  r1 + 4*w1xw1*w2xw2*pplxql*rx + 4*w1xw1*w2xw2*pptxppt*m**2 + 4
     &  *w1xw1*w2xw2*pptxppt*r1 - 4*w1xw1*w2xw2*pptxppt*rx - 4*w1xw1*
     &  w2xw2*qlxql*m**2 - 4*w1xw1*w2xw2*m**4 - 4*w1xw1*w2xw2*m**2*
     &  m0 - 4*w1xw1*w2xw2*m**2*r1 + 4*w1xw1*w2xw2*m**2*rx + 2*
     &  w1xw1*w2xw2*r1*rx - 2*w1xw1*w2xw2*rx**2 - 8*w1xw1*w2xpml*
     &  w2xppl*pplxppl - 8*w1xw1*w2xpml*w2xppl*pptxppt + 8*w1xw1*
     &  w2xpml*w2xppl*qlxql + 8*w1xw1*w2xpml*w2xppl*m**2
     &  + 8*w1xw1*w2xpml*w2xppl*m0 + 16*w1xw1*w2xpml*w2xql*
     &  pplxppl - 16*w1xw1*w2xpml*w2xql*pplxql + 16*w1xw1*w2xpml*
     &  w2xql*pptxppt - 16*w1xw1*w2xpml*w2xql*m**2 + 8*w1xw1*w2xpml*
     &  w2xql*rx - 16*w1xw2*w1xppl*w2xpml*pplxppl + 32*w1xw2*w1xppl*
     &  w2xpml*pplxql - 16*w1xw2*w1xppl*w2xpml*pptxppt - 16*w1xw2*
     &  w1xppl*w2xpml*qlxql + 16*w1xw2*w1xppl*w2xpml*m**2 - 16*w1xw2*
     &  w1xppl*w2xpml*m0 - 16*w1xw2*w1xppl*w2xpml*rx + 8*w1xpml*
     &  w1xppl*w2xw2*pplxppl - 16*w1xpml*w1xppl*w2xw2*pplxql + 8*
     &  w1xpml*w1xppl*w2xw2*pptxppt + 8*w1xpml*w1xppl*w2xw2*qlxql - 8
     &  *w1xpml*w1xppl*w2xw2*m**2 + 8*w1xpml*w1xppl*w2xw2*m0 + 8*
     &  w1xpml*w1xppl*w2xw2*rx - 16*w1xppl**2*w2xw2*pmlxppl + 16*
     &  w1xppl**2*w2xw2*pmlxql - 16*w1xppl**2*w2xw2*pmtxppt - 16*
     &  w1xppl**2*w2xw2*m**2 - 8*w1xppl**2*w2xw2*r1 + 8*w1xppl**2*
     &  w2xw2*rx + 32*w1xppl**2*w2xpml*w2xppl - 32*w1xppl**2*w2xpml*
     &  w2xql + 16*w1xppl*w1xql*w2xw2*pmlxppl - 16*w1xppl*w1xql*w2xw2
     &  *pmlxql + 16*w1xppl*w1xql*w2xw2*
     &  pmtxppt + 16*w1xppl*w1xql*w2xw2*m**2 + 8*w1xppl*w1xql*w2xw2*
     &  r1 - 8*w1xppl*w1xql*w2xw2*rx - 32*w1xppl*w1xql*w2xpml*w2xppl
     &  + 32*w1xppl*w1xql*w2xpml*w2xql )
      
      N17 = Iz2(k1,m0,m1) * (  - 8*w1xw2**2 )
      
      N18 = Id1(mkd1,mkx1,m1,md,mx) * ( 4*w1xw1*w2xw2 - 8*w1xw2**2 )
      
C
C dsigma is summ of all terms
C
      NT=N1+N2+N3+N4+N5+N6+N7+N8+N9+N10+N11+N12+N13+N14+
     &  N15+N16+N17+N18

C
C correction from w/u
      NT=NT*4D0/beta**2
C 1/(2pi)**6 from d3p, (2*pi)**2 from F.T.
      NT=NT/(2*pi)**6*(2*pi)**2
C from 1/2E+ 1/2E-
      NT=NT/4D0
C transform from MeV^-2 to kbarn
      NT=NT*(1.9733D0)**2/10D0
      
      dsigma=NT
      IF((setzero.EQ.1).OR.(dsigma.LT.0)) THEN
      dsigma=0D0
      badcount=badcount+1
      ENDIF 
      END
      
C========================================================================
C All the differential integral forms are calculated here
C
      DOUBLE PRECISION FUNCTION Iz0(x,u,v)
      IMPLICIT NONE
      
      DOUBLE PRECISION x(2)
      DOUBLE PRECISION u,v
      DOUBLE PRECISION s,arg
      
      INTEGER setzero
      COMMON /SETZPARAM/ setzero
      
      DOUBLE PRECISION tepxx
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793238462643D0)
      
      tepxx=x(1)*x(1)+x(2)*x(2)
      s=SQRT((tepxx+u+v)**2 - 4*u*v)
      arg=(tepxx+u+v+s)**2/(4*u*v)
      IF(arg.lt.0D0) THEN
      setzero=1
      Iz0=0D0
      ELSE 
      Iz0=pi * LOG((tepxx+u+v+s)**2/(4*u*v)) / s        
      ENDIF 
      END
      
C     ---------------------------------------------------------------     

      DOUBLE PRECISION FUNCTION Iz1(x,u,v)
      IMPLICIT NONE
      DOUBLE PRECISION Iz0
      EXTERNAL         Iz0

      DOUBLE PRECISION x(2)
      DOUBLE PRECISION u,v
      DOUBLE PRECISION a,b,s
      
      DOUBLE PRECISION tepxx
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793238462643D0)
      
      tepxx=x(1)*x(1)+x(2)*x(2)
      
      s=SQRT((tepxx+u+v)**2 - 4*u*v)
      A= (2*(tepxx+u-v+s)) / (s**2*(tepxx+u+v+s)) 
     &  - 1/(s*u)
      B= -Iz0(x,u,v)/pi * (tepxx+u-v)/s**2     
      Iz1=(-1)*pi*(A + B)
      
      END
      
C     -----------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION Iz2(x,u,v)
      IMPLICIT NONE
      
      DOUBLE PRECISION Iz0,Iz1
      EXTERNAL         Iz0,Iz1

      DOUBLE PRECISION x(2)
      DOUBLE PRECISION u,v
      DOUBLE PRECISION s,A,B,C,D,E
      
      DOUBLE PRECISION tepxx
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793238462643D0)

      tepxx=x(1)*x(1)+x(2)*x(2)
      
      s=SQRT((tepxx+u+v)**2 - 4*u*v)
      A= (2*(tepxx-u+v-s))/(s**3*(tepxx+u+v+s))
      B= -2*(tepxx+u-v+s)*( 
     &  2*(tepxx-u+v)*(tepxx+u+v+s) +
     &  s*(tepxx-u+v+s)  ) / 
     &  (s**4*(tepxx+u+v+s)**2)
      C= (tepxx-u+v)/(u*s**3)
      D= (Iz1(x,v,u)/pi)*(tepxx+u-v)/(s**2)
      E= (Iz0(x,u,v)/pi)*( 
     &  1/s**2 + (tepxx+u-v)*2*(tepxx-u+v) /
     &  (s**4) )
      Iz2= pi*(A + B + C +D + E)
      
      END

C     -------------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION Id0(x,y,u,v,w)
      IMPLICIT NONE
      DOUBLE PRECISION Iz0,Iz1
      EXTERNAL         Iz0,Iz1
      DOUBLE PRECISION x(2),y(2),dyx(2)
      DOUBLE PRECISION u,v,w
      DOUBLE PRECISION A,AXY,B,C,D,E,RX,RY
      DOUBLE PRECISION tepxx,tepxy,tepyy

      tepxx=x(1)*x(1)+x(2)*x(2)
      tepyy=y(1)*y(1)+y(2)*y(2)
      tepxy=x(1)*y(1)+x(2)*y(2)
      
      rx= v-u+tepxx
      ry= w-u+tepyy
      axy= x(1)*y(2)-x(2)*y(1)
      
      A= rx**2*tepyy - 2*rx*ry*tepxy + ry**2*tepxx
      B= 4*u*axy**2 + A
      C= 2*axy**2 + (rx+ry)*tepxy - rx*tepyy - ry*tepxx
      D= rx*tepyy - ry*tepxy
      E= ry*tepxx - rx*tepxy
      dyx(1)= y(1)-x(1)
      dyx(2)= y(2)-x(2)
      Id0=1/B*(C*Iz0(dyx,v,w) + D*Iz0(y,u,w) + E*Iz0(x,u,v))
      
      END
      
C     -------------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION Id1(x,y,u,v,w)
      IMPLICIT NONE
      DOUBLE PRECISION Id0,Iz0,Iz1
      EXTERNAL         Id0,Iz0,Iz1
      
      DOUBLE PRECISION x(2),y(2),dyx(2)
      DOUBLE PRECISION u,v,w
      
      DOUBLE PRECISION A,AXY,B,C,D,E,F,G,H,RX,RY
      DOUBLE PRECISION tepxx,tepxy,tepyy
      
      tepxx=x(1)*x(1)+x(2)*x(2)
      tepxy=x(1)*y(1)+x(2)*y(2)
      tepyy=y(1)*y(1)+y(2)*y(2)
      
      rx= v-u+tepxx
      ry= w-u+tepyy
      axy= x(1)*y(2)-x(2)*y(1)
      A= rx**2*tepyy - 2*rx*ry*tepxy + ry**2*tepxx
      B= 4*u*axy**2 + A
      C= -4*axy**2 -2*(-rx*tepyy+(rx+ry)*tepxy-ry*tepxx)
      D= -2*tepxy + tepxx + tepyy
      E= -tepyy + tepxy
      F= rx*tepyy - ry*tepxy
      G= -tepxx + tepxy
      H= ry*tepxx - rx*tepxy
      
      dyx(1)= y(1)-x(1)
      dyx(2)= y(2)-x(2)
      
      Id1= -( (C/B)*Id0(x,y,u,v,w) + (1/B)*(
     &  D*Iz0(dyx,v,w) +
     &  E*Iz0(y,u,w) - F*Iz1(y,u,w) +
     &  G*Iz0(x,u,v) - H*Iz1(x,u,v)   ) )
      
      END
      
C     --------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION Id2(x,y,u,v,w)
      IMPLICIT NONE
      DOUBLE PRECISION Id1,Id0,Iz0,Iz1,Iz2
      EXTERNAL         Id1,Id0,Iz0,Iz1,Iz2
      
      DOUBLE PRECISION x(2),y(2),mx(2),dyx(2)
      DOUBLE PRECISION u,v,w
      
      DOUBLE PRECISION A0,A1,A2,A3,A4,AXY
      DOUBLE PRECISION B,C1,C2,C3,D1,D2,D3
      DOUBLE PRECISION E1,E2,E3,F1,F2,F3
      DOUBLE PRECISION G1,G2,G3,RX,RY
      DOUBLE PRECISION tepxx,tepxy,tepyy
      
      tepxx=x(1)*x(1)+x(2)*x(2)
      tepxy=x(1)*y(1)+x(2)*y(2)
      tepyy=y(1)*y(1)+y(2)*y(2)
      
      mx(1)=-x(1)
      mx(2)=-x(2)
      dyx(1)= y(1)-x(1)
      dyx(2)= y(2)-x(2)
      rx= v-u+tepxx
      ry= w-u+tepyy
      axy= x(1)*y(2)-x(2)*y(1)
      
      A0= rx**2*tepyy - 2*rx*ry*tepxy + ry**2*tepxx
      B= 4*u*axy**2 + A0
      A1= 2*(tepyy-tepxy)*B
      A2= (-4*axy**2 +2*(rx*tepyy-(rx+ry)*tepxy+ry*tepxx))*
     &  2*2*(rx*tepyy-ry*tepxy)
      A3= -4*axy**2 +2*(rx*tepyy-(rx+ry)*tepxy+ry*tepxx)
      A4= -2*(rx*tepyy-ry*tepxy)
      C1= tepxy-tepyy
      D1= 2*axy**2+(rx+ry)*tepxy-rx*tepyy-ry*tepxx
      E1= tepyy
      F1= -tepxy
      G1= ry*tepxx-rx*tepxy
      C2= -2*tepxy+tepxx+tepyy
      D2= -tepyy + tepxy
      E2= rx*tepyy - ry*tepxy
      F2= -tepxx + tepxy
      G2= ry*tepxx - rx*tepxy
      C3= -2*tepxy+tepxx+tepyy
      D3= -tepyy
      E3= -tepxx+tepxy
      F3= tepxy
      G3= -(ry*tepxx-rx*tepxy)
      
      Id2= ((A1-A2)/B**2)*Id0(x,y,u,v,w) + 
     &  (A3/B**2)*( 
     &  C1*Iz0(dyx,v,w) - D1*Iz1(dyx,v,w) + 
     &  E1*Iz0(y,u,w) +
     &  F1*Iz0(x,u,v) - G1*Iz1(mx,v,u) )  +      
     &  (A4/B**2)*(                                        
     &  C2*Iz0(dyx,v,w) +
     &  D2*Iz0(y,u,w) - E2*Iz1(y,u,w) +
     &  F2*Iz0(x,u,v) - G2*Iz1(x,u,v)  ) +
     &  (1/B)*(
     &  - C3*Iz1(dyx,v,w) 
     &  + D3*Iz1(y,u,w)
     &  - E3*Iz1(mx,v,u) + F3*Iz1(x,u,v) - G3*Iz2(mx,v,u) )
      
      END
      
C     --------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION Id3(x,y,u,v,w)
      IMPLICIT NONE
      DOUBLE PRECISION Id2,Id1,Id0,Iz0,Iz1,Iz2
      EXTERNAL         Id2,Id1,Id0,Iz0,Iz1,Iz2
      
      DOUBLE PRECISION x(2),y(2),dyx(2),mdyx(2),mx(2),my(2)
      DOUBLE PRECISION u,v,w
      
      DOUBLE PRECISION A0,A1,A2,A3A,A3B,A3C,A3D
      DOUBLE PRECISION A4,A5A,A5B,A5C,A5D,A5E,A5F
      DOUBLE PRECISION A6,A8A,A8B,A8C,A9,A10,AXY
      DOUBLE PRECISION B,C2,C5,C6,C7,C8,C9,C10,C11
      DOUBLE PRECISION D2,D5,D6,D7,D8,D9,D10,D11
      DOUBLE PRECISION E2,E5,E6,E7,E8,E9,E10,E11
      DOUBLE PRECISION F2,F5,F6,F7,F8,F9,F10
      DOUBLE PRECISION G2,G5,G6,G7,G8,G9,G10
      DOUBLE PRECISION RX,RY
      DOUBLE PRECISION tepxx,tepxy,tepyy

      tepxx=x(1)*x(1)+x(2)*x(2)
      tepxy=x(1)*y(1)+x(2)*y(2)
      tepyy=y(1)*y(1)+y(2)*y(2)
      
      dyx(1)= y(1)-x(1)
      dyx(2)= y(2)-x(2)
      mdyx(1)=-dyx(1)
      mdyx(2)=-dyx(2)
      mx(1)=-x(1)
      mx(2)=-x(2)
      my(1)=-y(1)
      my(2)=-y(2)
      
      rx= v-u+tepxx
      ry= w-u+tepyy
      axy= x(1)*y(2)-x(2)*y(1)
      
      A0= rx**2*tepyy - 2*rx*ry*tepxy + ry**2*tepxx
      B= 4*u*axy**2 + A0
      
      A1= 2*(tepyy-tepxy)*(-2)*2*(-rx*tepxy+ry*tepxx)
      A2= 2*(tepyy-tepxy)
      C2= tepxy-tepxx
      D2= 2*axy**2+(rx+ry)*tepxy-rx*tepyy-ry*tepxx
      E2= -tepxy
      F2= rx*tepyy-ry*tepxy
      G2= tepxx
      
      A3a= 2*(-tepxy+tepxx) * 2*(rx*tepyy-ry*tepxy)     
      A3b= (-4*axy**2 +2*(rx*tepyy-(rx+ry)*tepxy+ry*tepxx))
     +  * 2*(-tepxy)
      A3c= (-4*axy**2 +2*(rx*tepyy-(rx+ry)*tepxy+ry*tepxx))
     &  * 2*(rx*tepyy-ry*tepxy)                  
      A3d= 3*2*(-rx*tepxy+ry*tepxx)
      
      A4= A3c
      A5a= 2*(-tepxy+tepxx)
      A5b= -4*axy**2 +2*(rx*tepyy-(rx+ry)*tepxy+ry*tepxx)
      A5c= 2*2*(-rx*tepxy+ ry*tepxx)
      A5d= -2*(rx*tepyy-ry*tepxy)
      A5e= -2*(-tepxy)
      A5f= -2*(rx*tepyy-ry*tepxy)*(2)*(-rx*tepxy+ry*tepxx)
      C5= tepxy-tepyy
      D5= 2*axy**2+(rx+ry)*tepxy-rx*tepyy-ry*tepxx
      E5= tepyy
      F5= -tepxy
      G5= ry*tepxx-rx*tepxy
      
      A6= -2*(rx*tepyy- ry*tepxy)
      C6= tepxy-tepxx
      D6= 2*axy**2+(rx+ry)*tepxy-rx*tepyy-ry*tepxx
      E6= -tepxy
      F6= rx*tepyy-ry*tepxy
      G6= tepxx
      
      C7= tepxy-tepyy
      D7= -(tepxy-tepxx)
      E7= 2*axy**2+(rx+ry)*tepxy-rx*tepyy-ry*tepxx
      F7= tepyy
      G7= -tepxx
      
      A8a= -2*(-tepxy)
      A8b= -2*(rx*tepyy-ry*tepxy)
      A8c=  2*2*(-rx*tepxy+ry*tepxx)
      C8= -2*tepxy+tepxx+tepyy
      D8= -tepyy + tepxy
      E8= rx*tepyy - ry*tepxy
      F8= -tepxx + tepxy
      G8= ry*tepxx - rx*tepxy
      
      A9= -2*(rx*tepyy-ry*tepxy)
      C9= -2*tepxy+tepxx+tepyy
      D9= -tepyy+tepxy  
      E9= tepxy
      F9= -(rx*tepyy-ry*tepxy)
      G9= -tepxx
      
      A10= -2*(-rx*tepxy+ry*tepxx)
      C10= -2*tepxy+tepxx+tepyy
      D10= -tepyy
      E10= -tepxx+tepxy
      F10= tepxy
      G10= -(ry*tepxx-rx*tepxy)
      
      C11= -(-2*tepxy+tepyy+tepxx)
      D11= -tepyy
      E11= tepxx
      
     	Id3= -( (A1/B**2)*Id0(x,y,u,v,w) + (A2/B**2)*(
     &  C2*Iz0(dyx,v,w) - D2*Iz1(mdyx,w,v)+
     &  E2*Iz0(y,u,w) - F2*Iz1(my,w,u) + 
     &  G2*Iz0(x,u,v) )      
     &  -( ((A3a+A3b)*B - A3c*A3d)/B**3 * Id0(x,y,u,v,w) +
     &  (A4/B**3)*(
     &  C2*Iz0(dyx,v,w) - D2*Iz1(mdyx,w,v)+
     &  E2*Iz0(y,u,w) - F2*Iz1(my,w,u)+
     &  G2*Iz0(x,u,v) ) ) +
     &  (A5a*B - A5b*A5c)/B**3*( A5d*Id0(x,y,u,v,w) + 
     &  C5*Iz0(dyx,v,w) - D5*Iz1(dyx,v,w) + 
     &  E5*Iz0(y,u,w) +
     &  F5*Iz0(x,u,v) - G5*Iz1(mx,v,u) ) +
     &  A5b/B**2 * ( (A5e - A5f/B) * Id0(x,y,u,v,w) + 
     &  (A6/B)*(  
     &  C6*Iz0(dyx,v,w) - D6*Iz1(mdyx,w,v) + 
     &  E6*Iz0(y,u,w) - F6*Iz1(my,w,u)+
     &  G6*Iz0(x,u,v)  ) 
     &  - C7*Iz1(mdyx,w,v) + D7*Iz1(dyx,v,w) 
     &  + E7*Iz2(mdyx,w,v)
     &  - F7*Iz1(my,w,u) 
     &  + G7*Iz1(mx,v,u) )+
     &  (A8a*B - A8b*A8c)/B**3*(
     &  C8*Iz0(dyx,v,w) +
     &  D8*Iz0(y,u,w) - E8*Iz1(y,u,w) + 
     &  F8*Iz0(x,u,v) - G8*Iz1(x,u,v)  ) +
     &  A9/B**2*(  
     &  - C9*Iz1(mdyx,w,v) 
     &  - D9*Iz1(my,w,u) + E9*Iz1(y,u,w) - F9*Iz2(y,u,w) 
     &  + G9*Iz1(x,u,v)  ) +
     &  A10/B**2*(                  
     &  - C10*Iz1(dyx,v,w) 
     &  + D10*Iz1(y,u,w)
     &  - E10*Iz1(mx,v,u) + F10*Iz1(x,u,v) - G10*Iz2(mx,v,u) )+
     &  1/B*(
     &  - C11*Iz2(dyx,v,w) 
     &  - D11*Iz2(y,u,w)
     &  + E11*Iz2(x,u,v)  ))
      
      END
      
C     --------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION Iv0(x,y,z,u,v,w,w2)
      IMPLICIT NONE
      DOUBLE PRECISION Id0
      EXTERNAL         Id0
      DOUBLE PRECISION x(2),y(2),z(2),dyx(2),dzx(2)
      DOUBLE PRECISION u,v,w,w2
      DOUBLE PRECISION rx,ry,rz
      DOUBLE PRECISION axy,ayz,azx
      DOUBLE PRECISION A

      rx= v-u+x(1)*x(1)+x(2)*x(2)
      ry= w-u+y(1)*y(1)+y(2)*y(2)
      rz= w2-u+z(1)*z(1)+z(2)*z(2)
      
      axy= x(1)*y(2)-x(2)*y(1)
      ayz= y(1)*z(2)-y(2)*z(1)
      azx= z(1)*x(2)-z(2)*x(1)
      
      dyx(1)= y(1)-x(1)
      dyx(2)= y(2)-x(2)
      dzx(1)= z(1)-x(1)
      dzx(2)= z(2)-x(2)
      
      A= 1/(axy*rz + ayz*rx + azx*ry)
      
      Iv0= A*(-(axy+ayz+azx)*Id0(dyx,dzx,v,w,w2)
     &  +axy*Id0(x,y,u,v,w)
     &  +ayz*Id0(y,z,u,w,w2)
     &  +azx*Id0(x,z,u,v,w2) )
      
      END
     
C     --------------------------------------------------------

      DOUBLE PRECISION FUNCTION Iv1(x,y,z,u,v,w,w2)
      IMPLICIT NONE
      DOUBLE PRECISION Iv0,Id0,Id1
      EXTERNAL         Iv0,Id0,Id1
      
      DOUBLE PRECISION x(2),y(2),z(2),dyx(2),dzx(2)
      DOUBLE PRECISION u,v,w,w2
      DOUBLE PRECISION A,B,RX,RY,RZ
      DOUBLE PRECISION AXY,AYZ,AZX
      
      dyx(1)= y(1)-x(1)
      dyx(2)= y(2)-x(2)
      dzx(1)= z(1)-x(1)
      dzx(2)= z(2)-x(2)
      
      rx= v-u+x(1)*x(1)+x(2)*x(2)
      ry= w-u+y(1)*y(1)+y(2)*y(2)
      rz= w2-u+z(1)*z(1)+z(2)*z(2)
      
      axy= x(1)*y(2)-x(2)*y(1)
      ayz= y(1)*z(2)-y(2)*z(1)
      azx= z(1)*x(2)-z(2)*x(1)
      
      A= 1/(axy*rz + ayz*rx + azx*ry)
      B= axy+ayz+azx
      Iv1=-((B*A**2)*(-(axy+ayz+azx)*Id0(dyx,dzx,v,w,w2)
     &  +axy*Id0(x,y,u,v,w)
     &  +ayz*Id0(y,z,u,w,w2)
     &  +azx*Id0(x,z,u,v,w2) )+
     &  A*( -axy*Id1(x,y,u,v,w)
     &  -ayz*Id1(y,z,u,w,w2)
     &  -azx*Id1(x,z,u,v,w2) ))
      
      END

C     --------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION Iv2(x,y,z,u,v,w,w2)
      IMPLICIT NONE
      DOUBLE PRECISION Iv1,Iv0,Id0,Id1,Id2
      EXTERNAL         Iv1,Iv0,Id0,Id1,Id2
      DOUBLE PRECISION x(2),mx(2),y(2),z(2),dyx(2)
      DOUBLE PRECISION dzx(2)
      DOUBLE PRECISION u,v,w,w2
      DOUBLE PRECISION axy,azx,ayz
      DOUBLE PRECISION rx,ry,rz,A,B
      
      mx(1)=-x(1)
      mx(2)=-x(2)
      dyx(1)= y(1)-x(1)
      dyx(2)= y(2)-x(2)
      dzx(1)= z(1)-x(1)
      dzx(2)= z(2)-x(2)
      
      rx= v-u+x(1)*x(1)+x(2)*x(2)
      ry= w-u+y(1)*y(1)+y(2)*y(2)
      rz= w2-u+z(1)*z(1)+z(2)*z(2)
      
      axy= x(1)*y(2)-x(2)*y(1)
      azx= z(1)*x(2)-z(2)*x(1)
      ayz= y(1)*z(2)-y(2)*z(1)
      
      A= 1/(axy*rz + ayz*rx + azx*ry)
      B= axy+ayz+azx
      
      Iv2=(-B*2D0*ayz*A**3)*(
     &  -B*Id0(dyx,dzx,v,w,w2)
     &  +axy*Id0(x,y,u,v,w)
     &  +ayz*Id0(y,z,u,w,w2)
     &  +azx*Id0(x,z,u,v,w2) )+
     &  (B*A**2)*(B*Id1(dyx,dzx,v,w,w2)
     &  -axy*Id1(mx,dyx,v,u,w)
     &  -azx*Id1(mx,dzx,v,u,w2) )+
     &  (-ayz*A**2)*(-axy*Id1(x,y,u,v,w)
     &  -ayz*Id1(y,z,u,w,w2)
     &  -azx*Id1(x,z,u,v,w2) )+
     &  A*(axy*Id2(x,y,u,v,w)
     &  +azx*Id2(x,z,u,v,w2) )
      
      END
