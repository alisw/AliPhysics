      subroutine WHITGevolvep(xin,qin,p2in,ip2in,pdf)
      include 'parmsetup.inc'
      real*8 xin,qin,q2in,p2in,pdf(-6:6),xval(45),qcdl4,qcdl5
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer nset
      
      save 
      call getnset(iset)
      call getnmem(iset,imem)
      
      if(imem.eq.1.or.imem.eq.0) then
        call SFWHI1(xin,qin,upv,dnv,usea,dsea,str,chm,glu)

      elseif(imem.eq.2) then
        call SFWHI2(xin,qin,upv,dnv,usea,dsea,str,chm,glu)

      elseif(imem.eq.3) then
        call SFWHI3(xin,qin,upv,dnv,usea,dsea,str,chm,glu)

      elseif(imem.eq.4) then
        call SFWHI4(xin,qin,upv,dnv,usea,dsea,str,chm,glu)

      elseif(imem.eq.5) then
        call SFWHI5(xin,qin,upv,dnv,usea,dsea,str,chm,glu)

      elseif(imem.eq.6) then
        call SFWHI6(xin,qin,upv,dnv,usea,dsea,str,chm,glu)

      else
        CONTINUE
      endif     

      pdf(-6)= 0.0d0
      pdf(6)= 0.0d0
      pdf(-5)= 0.0d0
      pdf(5 )= 0.0d0
      pdf(-4)= chm
      pdf(4 )= chm
      pdf(-3)= str
      pdf(3 )= str
      pdf(-2)= usea
      pdf(2 )= upv
      pdf(-1)= dsea
      pdf(1 )= dnv
      pdf(0 )= glu
      
      return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry WHITGread(nset)
      read(1,*)nmem(nset),ndef(nset)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry WHITGalfa(alfas,qalfa)
        call getnset(iset)
	call getnmem(iset,imem)
	call GetOrderAsM(iset,iord)
        call Getlam4M(iset,imem,qcdl4)
        call Getlam5M(iset,imem,qcdl5)
        call aspdflib(alfas,Qalfa,iord,qcdl5)

      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry WHITGinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry WHITGpdf(mem)
      call getnset(iset)
      call setnmem(iset,mem)
c      imem = mem
      return
c
 1000 format(5e13.5)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-------------------------------------------------------
      subroutine SFWHI1(ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL)
c-------------------------------------------------------
c     WHIT1 parton distribution in the photon
c
c     INPUT:  integer ic  : if ic=0 then qc=0
c                           else qc is calculated
c             DOUBLE PRECISION  Q2  : energy scale Q^2 (GeV^2)
c             DOUBLE PRECISION  x   : energy fraction
c
c     OUTPUT: DOUBLE PRECISION  qu  : up-quark dist.
c             DOUBLE PRECISION  qd  : down- or strange-quark dist.
c             DOUBLE PRECISION  qc  : charm-quark dist.
c             DOUBLE PRECISION  g   : gluon dist.
c-------------------------------------------------------
c     Modified by M.Tanaka on July 22, 1994.
c     The bug pointed out by M.Drees is fixed.
c-------------------------------------------------------
c     Modified by I.Watanabe on July 22, 1994.
c-------------------------------------------------------
      implicit none
      external WHIT1G
      double precision
     +       ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL
c arg
      integer ic
      DOUBLE PRECISION Q2,x
      DOUBLE PRECISION qu,qd,qc,g
c const
      DOUBLE PRECISION q42it,q52it,lam42,lam52
      DOUBLE PRECISION alinv,mc,PI
c local
      DOUBLE PRECISION qv,qsea,cv,cs,dcv,dcs
      DOUBLE PRECISION A0val,A1val,A2val,Bval,Cval,
     $       A0sea,B0sea,BB0sea,C0sea
      DOUBLE PRECISION A0dcv,A1dcv,A2dcv,A3dcv,Bdcv,Cdcv
      DOUBLE PRECISION Adcs, B0dcs, B1dcs, Cdcs
      DOUBLE PRECISION x1,x2,mc2q2
      DOUBLE PRECISION s,s2,s3,s4,prsccf,alstpi
      DOUBLE PRECISION WHIT1G
c parameters
      parameter(lam42=0.16d0, lam52=0.091411319d0)
      parameter(Q42IT=4.0d0, Q52IT=100.0d0)
      parameter(alinv=137.036d0, mc=1.5d0)
      parameter(pi=3.14159265358979323846d0)
      common /scale/ s,s2,s3,s4,prsccf
c
c begin
      x=ZX
      Q2=ZQ*ZQ
      ic=1
c
      x1=1.0d0-x
      x2=x**2
      mc2q2=mc**2/Q2
c
      if(Q2.lt.100.0d0) then
c  under 100 GeV^2
c
c  set  scale s
         if(Q2.lt.4.0d0) then
cccc  for under 4GeV^2 prescription
            s=  0.0d0
            prsccf =  log(Q2/LAM42)/ log(Q42IT/LAM42)
            alstpi = 6.0d0/25.0d0/ log(Q42IT/LAM42)
         else
            s=   log(  log(Q2/LAM42)/ log(Q42IT/LAM42))
            prsccf = 1.0d0
            alstpi = 6.0d0/25.0d0/ log(Q2/LAM42)
         endif
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT1 quark (U100)
c
      A0val= 1.882000d+00+s*( 1.213000d+00)+s2*( 6.970000d-01)
      A1val=              s*(-2.361000d+00)+s2*(-1.136000d+00)
      A2val=              s*( 5.280000d-01)+s2*( 2.406000d+00)
      Bval = 5.000000d-01+s*( 2.107000d-02)+s2*( 4.130000d-03)
      Cval = 2.500000d-01+s*(-2.376000d-01)+s2*( 2.018000d-01)
     $           +s3*(-5.040000d-02)
      A0sea= 6.510000d-01+s*( 1.291000d+00)+s2*(-4.470000d+00)
     $           +s3*( 5.140000d+00)+s4*(-2.091000d+00)
      B0sea=-3.820000d-02+s*( 9.010000d-02)+s2*(-1.356000d+00)
     $           +s3*( 1.582000d+00)+s4*(-6.440000d-01)
      BB0sea=2.084000d+00+s*( 7.740000d+00)+s2*(-2.970000d+01)
     $           +s3*( 3.860000d+01)+s4*(-1.705000d+01)
      C0sea= 7.000000d+00+s*(-1.608000d+01)+s2*( 4.670000d+01)
     $           +s3*(-5.710000d+01)+s4*( 2.386000d+01)
c
         qv  = prsccf/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= prsccf/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
            call WHIT1Q(x,mc2q2,cv,cs)
            qc = cv/alinv/2.0d0/PI + cs*alstpi
            qc  = qc*x
            ZCB=qc
         else
            qc = 0.0d0
            ZCB=qc
         endif
c
         g   = WHIT1G(x,Q2)
         g   = g*x
         ZGL=g
c
      else
c over 100 GeV^2
c
c set scale s
         s=   log(  log(Q2/LAM52)/ log(Q52IT/LAM52))
         prsccf = 1.0d0
         alstpi = 6.0d0/23.0d0/ log(Q2/LAM52)
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT1 quark (O100)
c
      A0val= 3.058000d+00+s*( 2.474000d+00)+s2*( 1.002000d+00)
      A1val=-2.182000d+00+s*(-4.480000d+00)+s2*(-2.251000d-01)
      A2val= 1.522000d+00+s*( 4.310000d+00)+s2*( 1.314000d+00)
      Bval = 5.170000d-01+s*( 4.040000d-02)+s2*(-2.100000d-02)
      Cval = 1.655000d-01+s*(-2.062000d-02)+s2*( 5.360000d-02)
      A0sea= 6.250000d-01+s*(-5.890000d-01)+s2*( 4.180000d+00)
     $           +s3*(-1.206000d+01)+s4*( 1.257000d+01)
      B0sea=-2.492000d-01+s*(-4.110000d-01)+s2*( 9.660000d-01)
     $           +s3*(-2.584000d+00)+s4*( 2.670000d+00)
      BB0sea=2.100000d+00+s*(-5.750000d+00)+s2*( 4.780000d+01)
     $           +s3*(-1.407000d+02)+s4*( 1.476000d+02)
      C0sea= 4.780000d+00+s*( 4.860000d+00)+s2*(-4.890000d+01)
     $           +s3*( 1.477000d+02)+s4*(-1.602000d+02)
c
         qv  = 1.0d0/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= 1.0d0/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
         g   = WHIT1G(x,Q2)
         g   = g*x
         ZGL=g
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
      A0dcv=              s*( 1.219000d-01)+s2*( 6.200000d+00)
     $           +s3*(-2.504000d+01)+s4*( 3.098000d+01)
      A1dcv=              s*( 1.913000d+00)+s2*(-7.690000d+01)
     $           +s3*( 3.180000d+02)+s4*(-3.920000d+02)
      A2dcv=              s*(-7.160000d+00)+s2*( 2.503000d+02)
     $           +s3*(-1.062000d+03)+s4*( 1.308000d+03)
      A3dcv=              s*( 3.190000d+00)+s2*(-2.301000d+02)
     $           +s3*( 1.012000d+03)+s4*(-1.250000d+03)
      Bdcv = 4.990000d-01+s*( 3.470000d+00)+s2*(-1.526000d+01)
     $           +s3*( 1.967000d+01)
      Cdcv = 3.290000d-01+s*( 8.240000d+00)+s2*(-3.800000d+01)
     $           +s3*( 4.630000d+01)
      Adcs =              s*(-1.815000d-02)+s2*( 2.043000d-03)
     $           +s3*(-4.130000d-03)
      B0dcs=-3.086000d-01+s*(-2.565000d-01)+s2*( 9.840000d-02)
      B1dcs= 1.376000d+00+s*(-4.630000d-01)+s2*( 1.232000d+00)
      Cdcs = 3.650000d+00+s*( 7.290000d-01)+s2*(-7.570000d+00)
     $           +s3*( 7.790000d+00)
c
         dcv = 1.0d0/alinv/x*
     $         (A0dcv+x*A1dcv+x2*A2dcv+x2*x*A3dcv) * x**Bdcv * x1**Cdcv
         dcs = 1.0d0/alinv/x*
     $         Adcs * x**(B0dcs+B1dcs*x) * x1**Cdcs
c
           call WHIT1Q(x,mc*mc/Q2,cv,cs)
           qc = cv/alinv/2.0d0/PI + cs*alstpi + dcs + dcv
           qc  = qc*x
           ZCB=qc
         else
           qc = 0.0d0
           ZCB=qc
         endif
      endif
c
      return
      end
c-------------------------------------------------------
      subroutine SFWHI2(ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL)
c-------------------------------------------------------
c     WHIT2 parton distribution in the photon
c
c     INPUT:  integer ic  : if ic=0 then qc=0
c                           else qc is calculated
c             DOUBLE PRECISION  Q2  : energy scale Q^2 (GeV^2)
c             DOUBLE PRECISION  x   : energy fraction
c
c     OUTPUT: DOUBLE PRECISION  qu  : up-quark dist.
c             DOUBLE PRECISION  qd  : down- or strange-quark dist.
c             DOUBLE PRECISION  qc  : charm-quark dist.
c             DOUBLE PRECISION  g   : gluon dist.
c-------------------------------------------------------
c     Modified by M.Tanaka on July 22, 1994.
c     The bug pointed out by M.Drees is fixed.
c-------------------------------------------------------
c     Modified by I.Watanabe on July 22, 1994.
c-------------------------------------------------------
      implicit none
      external WHIT2G
      double precision
     +       ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL
c arg
      integer ic
      DOUBLE PRECISION Q2,x
      DOUBLE PRECISION qu,qd,qc,g
c const
      DOUBLE PRECISION q42it,q52it,lam42,lam52
      DOUBLE PRECISION alinv,mc,PI
c local
      DOUBLE PRECISION qv,qsea,cv,cs,dcv,dcs
      DOUBLE PRECISION A0val,A1val,A2val,Bval,Cval,
     $       A0sea,B0sea,BB0sea,C0sea
      DOUBLE PRECISION A0dcv,A1dcv,A2dcv,A3dcv,Bdcv,Cdcv
      DOUBLE PRECISION Adcs, B0dcs, B1dcs, Cdcs
      DOUBLE PRECISION x1,x2,mc2q2
      DOUBLE PRECISION s,s2,s3,s4,prsccf,alstpi
      DOUBLE PRECISION WHIT2G
c parameters
      parameter(lam42=0.16d0, lam52=0.091411319d0)
      parameter(Q42IT=4.0d0, Q52IT=100.0d0)
      parameter(alinv=137.036d0, mc=1.5d0)
      parameter(pi=3.14159265358979323846d0)
      common /scale/ s,s2,s3,s4,prsccf
c
c begin
      x=ZX
      Q2=ZQ*ZQ
      ic=1
c
      x1=1.0d0-x
      x2=x**2
      mc2q2=mc**2/Q2
c
      if(Q2.lt.100.0d0) then
c  under 100 GeV^2
c
c  set  scale s
         if(Q2.lt.4.0d0) then
cccc  for under 4GeV^2 prescription
            s=  0.0d0
            prsccf =  log(Q2/LAM42)/ log(Q42IT/LAM42)
            alstpi = 6.0d0/25.0d0/ log(Q42IT/LAM42)
         else
            s=   log(  log(Q2/LAM42)/ log(Q42IT/LAM42))
            prsccf = 1.0d0
            alstpi = 6.0d0/25.0d0/ log(Q2/LAM42)
         endif
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT2 quark (U100)
c
      A0val= 1.882000d+00+s*( 1.213000d+00)+s2*( 6.970000d-01)
      A1val=              s*(-2.361000d+00)+s2*(-1.136000d+00)
      A2val=              s*( 5.280000d-01)+s2*( 2.406000d+00)
      Bval=  5.000000d-01+s*( 2.107000d-02)+s2*( 4.130000d-03)
      Cval=  2.500000d-01+s*(-2.376000d-01)+s2*( 2.018000d-01)
     $           +s3*(-5.040000d-02)
      A0sea= 1.237000d+00+s*( 3.390000d+00)+s2*(-1.075000d+01)
     $           +s3*( 1.246000d+01)+s4*(-5.580000d+00)
      B0sea=-7.270000d-02+s*( 1.748000d-01)+s2*(-1.392000d+00)
     $           +s3*( 1.711000d+00)+s4*(-7.960000d-01)
      BB0sea=4.290000d+00+s*( 1.787000d+01)+s2*(-5.810000d+01)
     $           +s3*( 8.190000d+01)+s4*(-4.140000d+01)
      C0sea= 1.434000d+01+s*(-4.490000d+01)+s2*( 1.197000d+02)
     $           +s3*(-1.585000d+02)+s4*( 7.530000d+01)
c
         qv  = prsccf/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= prsccf/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
            call WHIT2Q(x,mc2q2,cv,cs)
            qc = cv/alinv/2.0d0/PI + cs*alstpi
            qc  = qc*x
            ZCB=qc
         else
            qc = 0.0d0
            ZCB=qc
         endif
c
         g   = WHIT2G(x,Q2)
         g   = g*x
         ZGL=g
c
      else
c over 100 GeV^2
c
c set scale s
         s=   log(  log(Q2/LAM52)/ log(Q52IT/LAM52))
         prsccf = 1.0d0
         alstpi = 6.0d0/23.0d0/ log(Q2/LAM52)
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT2 quark (O100)
c
      A0val= 3.058000d+00+s*( 2.474000d+00)+s2*( 1.002000d+00)
      A1val=-2.182000d+00+s*(-4.480000d+00)+s2*(-2.259000d-01)
      A2val= 1.522000d+00+s*( 4.300000d+00)+s2*( 1.315000d+00)
      Bval = 5.170000d-01+s*( 4.030000d-02)+s2*(-2.098000d-02)
      Cval = 1.655000d-01+s*(-2.063000d-02)+s2*( 5.370000d-02)
      A0sea= 1.287000d+00+s*(-2.069000d+00)+s2*( 1.157000d+01)
     $           +s3*(-3.570000d+01)+s4*( 3.740000d+01)
      B0sea=-2.340000d-01+s*(-4.430000d-01)+s2*( 1.235000d+00)
     $           +s3*(-3.720000d+00)+s4*( 3.840000d+00)
      BB0sea=6.460000d+00+s*(-1.048000d+01)+s2*( 8.980000d+01)
     $           +s3*(-2.847000d+02)+s4*( 2.998000d+02)
      C0sea= 5.350000d+00+s*( 1.011000d+01)+s2*(-1.337000d+02)
     $           +s3*( 4.270000d+02)+s4*(-4.570000d+02)
c
         qv  = 1.0d0/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= 1.0d0/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
         g   = WHIT2G(x,Q2)
         g   = g*x
         ZGL=g
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
      A0dcv=              s*( 1.219000d-01)+s2*( 6.200000d+00)
     $           +s3*(-2.504000d+01)+s4*( 3.098000d+01)
      A1dcv=              s*( 1.913000d+00)+s2*(-7.690000d+01)
     $           +s3*( 3.180000d+02)+s4*(-3.920000d+02)
      A2dcv=              s*(-7.160000d+00)+s2*( 2.503000d+02)
     $           +s3*(-1.062000d+03)+s4*( 1.308000d+03)
      A3dcv=              s*( 3.190000d+00)+s2*(-2.301000d+02)
     $           +s3*( 1.012000d+03)+s4*(-1.250000d+03)
      Bdcv = 4.990000d-01+s*( 3.470000d+00)+s2*(-1.526000d+01)
     $           +s3*( 1.967000d+01)
      Cdcv = 3.290000d-01+s*( 8.240000d+00)+s2*(-3.800000d+01)
     $           +s3*( 4.630000d+01)
      Adcs =              s*(-2.786000d-02)+s2*( 3.490000d-02)
     $           +s3*(-2.223000d-02)
      B0dcs=-3.141000d-01+s*(-4.250000d-01)+s2*( 1.564000d-01)
      B1dcs= 4.720000d+00+s*(-5.480000d+00)+s2*( 2.686000d+00)
      Cdcs = 2.961000d+00+s*( 7.760000d-01)+s2*(-8.280000d+00)
     $           +s3*( 9.780000d+00)
c
         dcv = 1.0d0/alinv/x*
     $         (A0dcv+x*A1dcv+x2*A2dcv+x2*x*A3dcv) * x**Bdcv * x1**Cdcv
         dcs = 1.0d0/alinv/x*
     $         Adcs * x**(B0dcs+B1dcs*x) * x1**Cdcs
c
           call WHIT2Q(x,mc*mc/Q2,cv,cs)
           qc = cv/alinv/2.0d0/PI + cs*alstpi + dcs + dcv
           qc  = qc*x
           ZCB=qc
         else
           qc = 0.0d0
           ZCB=qc
         endif
      endif
c
      return
      end
c-------------------------------------------------------
      subroutine SFWHI3(ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL)
c-------------------------------------------------------
c     WHIT3 parton distribution in the photon
c
c     INPUT:  integer ic  : if ic=0 then qc=0
c                           else qc is calculated
c             DOUBLE PRECISION  Q2  : energy scale Q^2 (GeV^2)
c             DOUBLE PRECISION  x   : energy fraction
c
c     OUTPUT: DOUBLE PRECISION  qu  : up-quark dist.
c             DOUBLE PRECISION  qd  : down- or strange-quark dist.
c             DOUBLE PRECISION  qc  : charm-quark dist.
c             DOUBLE PRECISION  g   : gluon dist.
c-------------------------------------------------------
c     Modified by M.Tanaka on July 22, 1994.
c     The bug pointed out by M.Drees is fixed.
c-------------------------------------------------------
c     Modified by I.Watanabe on July 22, 1994.
c-------------------------------------------------------
      implicit none
      external whit3g
      double precision
     +       ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL
c arg
      integer ic
      DOUBLE PRECISION Q2,x
      DOUBLE PRECISION qu,qd,qc,g
c const
      DOUBLE PRECISION q42it,q52it,lam42,lam52
      DOUBLE PRECISION alinv,mc,PI
c local
      DOUBLE PRECISION qv,qsea,cv,cs,dcv,dcs
      DOUBLE PRECISION A0val,A1val,A2val,Bval,Cval,
     $       A0sea,B0sea,BB0sea,C0sea
      DOUBLE PRECISION A0dcv,A1dcv,A2dcv,A3dcv,Bdcv,Cdcv
      DOUBLE PRECISION Adcs, B0dcs, B1dcs, Cdcs
      DOUBLE PRECISION x1,x2,mc2q2
      DOUBLE PRECISION s,s2,s3,s4,prsccf,alstpi
      DOUBLE PRECISION WHIT3g
c parameters
      parameter(lam42=0.16d0, lam52=0.091411319d0)
      parameter(Q42IT=4.0d0, Q52IT=100.0d0)
      parameter(alinv=137.036d0, mc=1.5d0)
      parameter(pi=3.14159265358979323846d0)
      common /scale/ s,s2,s3,s4,prsccf
c
c begin
      x=ZX
      Q2=ZQ*ZQ
      ic=1
c
      x1=1.0d0-x
      x2=x**2
      mc2q2=mc**2/Q2
c
      if(Q2.lt.100.0d0) then
c  under 100 GeV^2
c
c  set  scale s
         if(Q2.lt.4.0d0) then
cccc  for under 4GeV^2 prescription
            s=  0.0d0
            prsccf =  log(Q2/LAM42)/ log(Q42IT/LAM42)
            alstpi = 6.0d0/25.0d0/ log(Q42IT/LAM42)
         else
            s=   log(  log(Q2/LAM42)/ log(Q42IT/LAM42))
            prsccf = 1.0d0
            alstpi = 6.0d0/25.0d0/ log(Q2/LAM42)
         endif
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT3 quark (U100)
c
      A0val= 1.882000d+00+s*( 1.213000d+00)+s2*( 6.970000d-01)
      A1val=              s*(-2.361000d+00)+s2*(-1.136000d+00)
      A2val=              s*( 5.280000d-01)+s2*( 2.406000d+00)
      Bval = 5.000000d-01+s*( 2.107000d-02)+s2*( 4.130000d-03)
      Cval = 2.500000d-01+s*(-2.376000d-01)+s2*( 2.018000d-01)
     $           +s3*(-5.040000d-02)
      A0sea= 1.587000d+00+s*( 5.050000d+00)+s2*(-1.126000d+01)
     $           +s3*( 7.560000d+00)+s4*(-1.471000d+00)
      B0sea=-1.006000d-01+s*( 2.259000d-01)+s2*(-1.195000d+00)
     $           +s3*( 1.175000d+00)+s4*(-4.460000d-01)
      BB0sea=5.730000d+00+s*( 2.564000d+01)+s2*(-5.870000d+01)
     $           +s3*( 6.320000d+01)+s4*(-2.577000d+01)
      C0sea= 2.136000d+01+s*(-7.290000d+01)+s2*( 1.532000d+02)
     $           +s3*(-1.679000d+02)+s4*( 6.740000d+01)
c
         qv  = prsccf/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= prsccf/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
            call WHIT3Q(x,mc2q2,cv,cs)
            qc = cv/alinv/2.0d0/PI + cs*alstpi
            qc  = qc*x
            ZCB=qc
         else
            qc = 0.0d0
            ZCB=qc
         endif
c
         g   = WHIT3G(x,Q2)
         g   = g*x
         ZGL=g
c
      else
c over 100 GeV^2
c
c set scale s
         s=   log(  log(Q2/LAM52)/ log(Q52IT/LAM52))
         prsccf = 1.0d0
         alstpi = 6.0d0/23.0d0/ log(Q2/LAM52)
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT3 quark (O100)
c
      A0val= 3.058000d+00+s*( 2.474000d+00)+s2*( 1.002000d+00)
      A1val=-2.182000d+00+s*(-4.480000d+00)+s2*(-2.264000d-01)
      A2val= 1.522000d+00+s*( 4.300000d+00)+s2*( 1.315000d+00)
      Bval = 5.170000d-01+s*( 4.030000d-02)+s2*(-2.097000d-02)
      Cval = 1.655000d-01+s*(-2.064000d-02)+s2*( 5.370000d-02)
      A0sea= 1.850000d+00+s*(-3.670000d+00)+s2*( 2.714000d+01)
     $           +s3*(-1.066000d+02)+s4*( 1.309000d+02)
      B0sea=-2.299000d-01+s*(-4.970000d-01)+s2*( 2.464000d+00)
     $           +s3*(-9.950000d+00)+s4*( 1.232000d+01)
      BB0sea=1.042000d+01+s*(-1.074000d+01)+s2*( 1.327000d+02)
     $           +s3*(-5.390000d+02)+s4*( 6.560000d+02)
      C0sea= 4.070000d+00+s*( 4.110000d+00)+s2*(-1.719000d+02)
     $           +s3*( 7.070000d+02)+s4*(-8.590000d+02)
c
         qv  = 1.0d0/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= 1.0d0/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
         g   = WHIT3G(x,Q2)
         g   = g*x
         ZGL=g
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
      A0dcv=              s*( 1.219000d-01)+s2*( 6.200000d+00)
     $           +s3*(-2.504000d+01)+s4*( 3.098000d+01)
      A1dcv=              s*( 1.913000d+00)+s2*(-7.690000d+01)
     $           +s3*( 3.180000d+02)+s4*(-3.920000d+02)
      A2dcv=              s*(-7.160000d+00)+s2*( 2.503000d+02)
     $           +s3*(-1.062000d+03)+s4*( 1.308000d+03)
      A3dcv=              s*( 3.190000d+00)+s2*(-2.301000d+02)
     $           +s3*( 1.012000d+03)+s4*(-1.250000d+03)
      Bdcv = 4.990000d-01+s*( 3.470000d+00)+s2*(-1.526000d+01)
     $           +s3*( 1.967000d+01)
      Cdcv = 3.290000d-01+s*( 8.240000d+00)+s2*(-3.800000d+01)
     $           +s3*( 4.630000d+01)
      Adcs =              s*(-1.948000d-02)+s2*( 2.861000d-02)
     $           +s3*(-2.036000d-02)
      B0dcs=-4.130000d-01+s*(-4.390000d-01)+s2*( 1.810000d-01)
      B1dcs= 5.190000d+00+s*(-7.400000d+00)+s2*( 3.400000d+00)
      Cdcs = 2.359000d+00+s*( 9.770000d-01)+s2*(-7.730000d+00)
     $           +s3*( 9.480000d+00)
c
         dcv = 1.0d0/alinv/x*
     $         (A0dcv+x*A1dcv+x2*A2dcv+x2*x*A3dcv) * x**Bdcv * x1**Cdcv
         dcs = 1.0d0/alinv/x*
     $         Adcs * x**(B0dcs+B1dcs*x) * x1**Cdcs
c
           call WHIT3Q(x,mc*mc/Q2,cv,cs)
           qc = cv/alinv/2.0d0/PI + cs*alstpi + dcs + dcv
           qc  = qc*x
           ZCB=qc
         else
           qc = 0.0d0
           ZCB=qc
         endif
      endif
c
      return
      end
c-------------------------------------------------------
      subroutine SFWHI4(ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL)
c-------------------------------------------------------
c     WHIT4 parton distribution in the photon
c
c     INPUT:  integer ic  : if ic=0 then qc=0
c                           else qc is calculated
c             DOUBLE PRECISION  Q2  : energy scale Q^2 (GeV^2)
c             DOUBLE PRECISION  x   : energy fraction
c
c     OUTPUT: DOUBLE PRECISION  qu  : up-quark dist.
c             DOUBLE PRECISION  qd  : down- or strange-quark dist.
c             DOUBLE PRECISION  qc  : charm-quark dist.
c             DOUBLE PRECISION  g   : gluon dist.
c-------------------------------------------------------
c     Modified by M.Tanaka on July 22, 1994.
c     The bug pointed out by M.Drees is fixed.
c-------------------------------------------------------
c     Modified by I.Watanabe on July 22, 1994.
c-------------------------------------------------------
      implicit none
      external WHIT4G
      double precision
     +       ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL
c arg
      integer ic
      DOUBLE PRECISION Q2,x
      DOUBLE PRECISION qu,qd,qc,g
c const
      DOUBLE PRECISION q42it,q52it,lam42,lam52
      DOUBLE PRECISION alinv,mc,PI
c local
      DOUBLE PRECISION qv,qsea,cv,cs,dcv,dcs
      DOUBLE PRECISION A0val,A1val,A2val,Bval,Cval,
     $       A0sea,B0sea,BB0sea,C0sea
      DOUBLE PRECISION A0dcv,A1dcv,A2dcv,A3dcv,Bdcv,Cdcv
      DOUBLE PRECISION Adcs, B0dcs, B1dcs, Cdcs
      DOUBLE PRECISION x1,x2,mc2q2
      DOUBLE PRECISION s,s2,s3,s4,prsccf,alstpi
      DOUBLE PRECISION WHIT4G
c parameters
      parameter(lam42=0.16d0, lam52=0.091411319d0)
      parameter(Q42IT=4.0d0, Q52IT=100.0d0)
      parameter(alinv=137.036d0, mc=1.5d0)
      parameter(pi=3.14159265358979323846d0)
      common /scale/ s,s2,s3,s4,prsccf
c
c begin
      x=ZX
      Q2=ZQ*ZQ
      ic=1
c
      x1=1.0d0-x
      x2=x**2
      mc2q2=mc**2/Q2
c
      if(Q2.lt.100.0d0) then
c  under 100 GeV^2
c
c  set  scale s
         if(Q2.lt.4.0d0) then
cccc  for under 4GeV^2 prescription
            s=  0.0d0
            prsccf =  log(Q2/LAM42)/ log(Q42IT/LAM42)
            alstpi = 6.0d0/25.0d0/ log(Q42IT/LAM42)
         else
            s=   log(  log(Q2/LAM42)/ log(Q42IT/LAM42))
            prsccf = 1.0d0
            alstpi = 6.0d0/25.0d0/ log(Q2/LAM42)
         endif
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT4 quark (U100)
c
      A0val= 2.540000d+00+s*( 2.000000d+00)+s2*( 7.180000d-01)
      A1val= 6.230000d-02+s*(-7.010000d+00)+s2*( 1.251000d-01)
      A2val=-1.642000d-01+s*(-4.360000d-01)+s2*( 1.048000d+01)
     $           +s3*(-5.200000d+00)
      Bval = 6.990000d-01+s*(-2.796000d-02)+s2*(-3.650000d-03)
      Cval = 4.420000d-01+s*(-1.255000d+00)+s2*( 1.941000d+00)
     $           +s3*(-9.950000d-01)
      A0sea= 1.308000d+00+s*( 2.315000d+00)+s2*(-7.880000d+00)
     $           +s3*( 8.260000d+00)+s4*(-3.004000d+00)
      B0sea=-3.730000d-02+s*( 5.630000d-02)+s2*(-1.133000d+00)
     $           +s3*( 1.185000d+00)+s4*(-4.180000d-01)
      BB0sea=2.103000d+00+s*( 4.850000d+00)+s2*(-1.781000d+01)
     $           +s3*( 2.062000d+01)+s4*(-7.940000d+00)
      C0sea= 7.000000d+00+s*(-1.017000d+01)+s2*( 2.600000d+01)
     $           +s3*(-2.960000d+01)+s4*( 1.227000d+01)
c
         qv  = prsccf/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= prsccf/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
            call WHIT4Q(x,mc2q2,cv,cs)
            qc = cv/alinv/2.0d0/PI + cs*alstpi
            qc  = qc*x
            ZCB=qc
         else
            qc = 0.0d0
            ZCB=qc
         endif
c
         g   = WHIT4G(x,Q2)
         g   = g*x
         ZGL=g
c
      else
c over 100 GeV^2
c
c set scale s
         s=   log(  log(Q2/LAM52)/ log(Q52IT/LAM52))
         prsccf = 1.0d0
         alstpi = 6.0d0/23.0d0/ log(Q2/LAM52)
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT4 quark (O100)
c
      A0val= 4.270000d+00+s*( 3.096000d+00)+s2*( 1.619000d+00)
      A1val=-4.740000d+00+s*(-6.900000d+00)+s2*(-2.430000d+00)
      A2val= 2.837000d+00+s*( 6.470000d+00)+s2*( 4.090000d+00)
      Bval = 6.780000d-01+s*(-3.940000d-02)+s2*( 1.756000d-02)
      Cval = 1.728000d-01+s*(-2.479000d-02)+s2*( 1.446000d-01)
      A0sea= 1.188000d+00+s*(-1.396000d+00)+s2*( 8.710000d+00)
     $           +s3*(-2.542000d+01)+s4*( 2.492000d+01)
      B0sea=-2.448000d-01+s*(-4.190000d-01)+s2*( 1.007000d+00)
     $           +s3*(-2.689000d+00)+s4*( 2.517000d+00)
      BB0sea=1.942000d+00+s*(-6.040000d+00)+s2*( 5.030000d+01)
     $           +s3*(-1.478000d+02)+s4*( 1.481000d+02)
      C0sea= 5.420000d+00+s*( 6.110000d+00)+s2*(-5.380000d+01)
     $           +s3*( 1.632000d+02)+s4*(-1.716000d+02)
c
         qv  = 1.0d0/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= 1.0d0/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
         g   = WHIT4G(x,Q2)
         g   = g*x
         ZGL=g
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
      A0dcv=              s*( 1.219000d-01)+s2*( 6.200000d+00)
     $           +s3*(-2.504000d+01)+s4*( 3.098000d+01)
      A1dcv=              s*( 1.913000d+00)+s2*(-7.690000d+01)
     $           +s3*( 3.180000d+02)+s4*(-3.920000d+02)
      A2dcv=              s*(-7.160000d+00)+s2*( 2.503000d+02)
     $           +s3*(-1.062000d+03)+s4*( 1.308000d+03)
      A3dcv=              s*( 3.190000d+00)+s2*(-2.301000d+02)
     $           +s3*( 1.012000d+03)+s4*(-1.250000d+03)
      Bdcv = 4.990000d-01+s*( 3.470000d+00)+s2*(-1.526000d+01)
     $           +s3*( 1.967000d+01)
      Cdcv = 3.290000d-01+s*( 8.240000d+00)+s2*(-3.800000d+01)
     $           +s3*( 4.630000d+01)
      Adcs =              s*(-2.821000d-02)+s2*(-2.649000d-04)
     $           +s3*( 7.040000d-03)
      B0dcs=-3.270000d-01+s*(-2.298000d-01)+s2*( 3.500000d-02)
      B1dcs= 1.254000d+00+s*( 8.780000d-01)+s2*( 2.086000d-01)
      Cdcs = 4.170000d+00+s*( 6.400000d-01)+s2*(-7.630000d+00)
     $           +s3*( 7.170000d+00)
c
         dcv = 1.0d0/alinv/x*
     $         (A0dcv+x*A1dcv+x2*A2dcv+x2*x*A3dcv) * x**Bdcv * x1**Cdcv
         dcs = 1.0d0/alinv/x*
     $         Adcs * x**(B0dcs+B1dcs*x) * x1**Cdcs
c
           call WHIT4Q(x,mc*mc/Q2,cv,cs)
           qc = cv/alinv/2.0d0/PI + cs*alstpi + dcs + dcv
           qc  = qc*x
           ZCB=qc
         else
           qc = 0.0d0
           ZCB=qc
         endif
      endif
c
      return
      end
c-------------------------------------------------------
      subroutine SFWHI5(ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL)
c-------------------------------------------------------
c     WHIT5 parton distribution in the photon
c
c     INPUT:  integer ic  : if ic=0 then qc=0
c                           else qc is calculated
c             DOUBLE PRECISION  Q2  : energy scale Q^2 (GeV^2)
c             DOUBLE PRECISION  x   : energy fraction
c
c     OUTPUT: DOUBLE PRECISION  qu  : up-quark dist.
c             DOUBLE PRECISION  qd  : down- or strange-quark dist.
c             DOUBLE PRECISION  qc  : charm-quark dist.
c             DOUBLE PRECISION  g   : gluon dist.
c-------------------------------------------------------
c     Modified by M.Tanaka on July 22, 1994.
c     The bug pointed out by M.Drees is fixed.
c-------------------------------------------------------
c     Modified by I.Watanabe on July 22, 1994.
c-------------------------------------------------------
      implicit none
      external WHIT5G
      double precision
     +       ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL
c arg
      integer ic
      DOUBLE PRECISION Q2,x
      DOUBLE PRECISION qu,qd,qc,g
c const
      DOUBLE PRECISION q42it,q52it,lam42,lam52
      DOUBLE PRECISION alinv,mc,PI
c local
      DOUBLE PRECISION qv,qsea,cv,cs,dcv,dcs
      DOUBLE PRECISION A0val,A1val,A2val,Bval,Cval,
     $       A0sea,B0sea,BB0sea,C0sea
      DOUBLE PRECISION A0dcv,A1dcv,A2dcv,A3dcv,Bdcv,Cdcv
      DOUBLE PRECISION Adcs, B0dcs, B1dcs, Cdcs
      DOUBLE PRECISION x1,x2,mc2q2
      DOUBLE PRECISION s,s2,s3,s4,prsccf,alstpi
      DOUBLE PRECISION WHIT5G
c parameters
      parameter(lam42=0.16d0, lam52=0.091411319d0)
      parameter(Q42IT=4.0d0, Q52IT=100.0d0)
      parameter(alinv=137.036d0, mc=1.5d0)
      parameter(pi=3.14159265358979323846d0)
      common /scale/ s,s2,s3,s4,prsccf
c
c begin
      x=ZX
      Q2=ZQ*ZQ
      ic=1
c
      x1=1.0d0-x
      x2=x**2
      mc2q2=mc**2/Q2
c
      if(Q2.lt.100.0d0) then
c  under 100 GeV^2
c
c  set  scale s
         if(Q2.lt.4.0d0) then
cccc  for under 4GeV^2 prescription
            s=  0.0d0
            prsccf =  log(Q2/LAM42)/ log(Q42IT/LAM42)
            alstpi = 6.0d0/25.0d0/ log(Q42IT/LAM42)
         else
            s=   log(  log(Q2/LAM42)/ log(Q42IT/LAM42))
            prsccf = 1.0d0
            alstpi = 6.0d0/25.0d0/ log(Q2/LAM42)
         endif
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT5 quark (U100)
c
      A0val= 2.540000d+00+s*( 2.000000d+00)+s2*( 7.180000d-01)
      A1val= 6.230000d-02+s*(-7.010000d+00)+s2*( 1.251000d-01)
      A2val=-1.642000d-01+s*(-4.360000d-01)+s2*( 1.048000d+01)
     $           +s3*(-5.200000d+00)
      Bval = 6.990000d-01+s*(-2.796000d-02)+s2*(-3.650000d-03)
      Cval = 4.420000d-01+s*(-1.255000d+00)+s2*( 1.941000d+00)
     $           +s3*(-9.950000d-01)
      A0sea= 2.227000d+00+s*( 5.720000d+00)+s2*(-1.295000d+01)
     $           +s3*( 7.220000d+00)+s4*(-2.514000d-01)
      B0sea=-8.810000d-02+s*( 1.465000d-01)+s2*(-9.750000d-01)
     $           +s3*( 7.820000d-01)+s4*(-2.074000d-01)
      BB0sea=3.370000d+00+s*( 1.416000d+01)+s2*(-3.150000d+01)
     $           +s3*( 2.789000d+01)+s4*(-8.710000d+00)
      C0sea= 1.581000d+01+s*(-3.630000d+01)+s2*( 7.710000d+01)
     $           +s3*(-7.810000d+01)+s4*( 2.948000d+01)
c
         qv  = prsccf/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= prsccf/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
            call WHIT5Q(x,mc2q2,cv,cs)
            qc = cv/alinv/2.0d0/PI + cs*alstpi
            qc  = qc*x
            ZCB=qc
         else
            qc = 0.0d0
            ZCB=qc
         endif
c
         g   = WHIT5G(x,Q2)
         g   = g*x
         ZGL=g
c
      else
c over 100 GeV^2
c
c set scale s
         s=   log(  log(Q2/LAM52)/ log(Q52IT/LAM52))
         prsccf = 1.0d0
         alstpi = 6.0d0/23.0d0/ log(Q2/LAM52)
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT5 quark (O100)
c
      A0val= 4.270000d+00+s*( 3.096000d+00)+s2*( 1.617000d+00)
      A1val=-4.740000d+00+s*(-6.900000d+00)+s2*(-2.417000d+00)
      A2val= 2.837000d+00+s*( 6.470000d+00)+s2*( 4.070000d+00)
      Bval = 6.780000d-01+s*(-3.940000d-02)+s2*( 1.750000d-02)
      Cval = 1.728000d-01+s*(-2.457000d-02)+s2*( 1.440000d-01)
      A0sea= 2.318000d+00+s*(-3.760000d+00)+s2*( 2.026000d+01)
     $           +s3*(-5.950000d+01)+s4*( 5.900000d+01)
      B0sea=-2.425000d-01+s*(-4.360000d-01)+s2*( 1.241000d+00)
     $           +s3*(-3.510000d+00)+s4*( 3.360000d+00)
      BB0sea=5.330000d+00+s*(-8.680000d+00)+s2*( 7.420000d+01)
     $           +s3*(-2.070000d+02)+s4*( 1.967000d+02)
      C0sea= 8.480000d+00+s*( 9.310000d+00)+s2*(-1.041000d+02)
     $           +s3*( 2.801000d+02)+s4*(-2.663000d+02)
c
         qv  = 1.0d0/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= 1.0d0/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
         g   = WHIT5G(x,Q2)
         g   = g*x
         ZGL=g
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
      A0dcv=              s*( 1.219000d-01)+s2*( 6.200000d+00)
     $           +s3*(-2.504000d+01)+s4*( 3.098000d+01)
      A1dcv=              s*( 1.913000d+00)+s2*(-7.690000d+01)
     $           +s3*( 3.180000d+02)+s4*(-3.920000d+02)
      A2dcv=              s*(-7.160000d+00)+s2*( 2.503000d+02)
     $           +s3*(-1.062000d+03)+s4*( 1.308000d+03)
      A3dcv=              s*( 3.190000d+00)+s2*(-2.301000d+02)
     $           +s3*( 1.012000d+03)+s4*(-1.250000d+03)
      Bdcv = 4.990000d-01+s*( 3.470000d+00)+s2*(-1.526000d+01)
     $           +s3*( 1.967000d+01)
      Cdcv = 3.290000d-01+s*( 8.240000d+00)+s2*(-3.800000d+01)
     $           +s3*( 4.630000d+01)
      Adcs =              s*(-6.580000d-02)+s2*( 1.059000d-01)
     $           +s3*(-6.630000d-02)
      B0dcs=-2.750000d-01+s*(-4.760000d-01)+s2*( 1.191000d-01)
      B1dcs= 6.370000d+00+s*(-5.320000d+00)+s2*( 1.986000d+00)
      Cdcs = 3.400000d+00+s*( 3.750000d-01)+s2*(-8.790000d+00)
     $           +s3*( 1.001000d+01)
c
         dcv = 1.0d0/alinv/x*
     $         (A0dcv+x*A1dcv+x2*A2dcv+x2*x*A3dcv) * x**Bdcv * x1**Cdcv
         dcs = 1.0d0/alinv/x*
     $         Adcs * x**(B0dcs+B1dcs*x) * x1**Cdcs
c
           call WHIT5Q(x,mc*mc/Q2,cv,cs)
           qc = cv/alinv/2.0d0/PI + cs*alstpi + dcs + dcv
           qc  = qc*x
           ZCB=qc
         else
           qc = 0.0d0
           ZCB=qc
         endif
      endif
c
      return
      end
c-------------------------------------------------------
      subroutine SFWHI6(ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL)
c-------------------------------------------------------
c     WHIT6 parton distribution in the photon
c
c     INPUT:  integer ic  : if ic=0 then qc=0
c                           else qc is calculated
c             DOUBLE PRECISION  Q2  : energy scale Q^2 (GeV^2)
c             DOUBLE PRECISION  x   : energy fraction
c
c     OUTPUT: DOUBLE PRECISION  qu  : up-quark dist.
c             DOUBLE PRECISION  qd  : down- or strange-quark dist.
c             DOUBLE PRECISION  qc  : charm-quark dist.
c             DOUBLE PRECISION  g   : gluon dist.
c-------------------------------------------------------
c     Modified by M.Tanaka on July 22, 1994.
c     The bug pointed out by M.Drees is fixed.
c-------------------------------------------------------
c     Modified by I.Watanabe on July 22, 1994.
c-------------------------------------------------------
      implicit none
      external WHIT6G
      double precision
     +       ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZGL
c arg
      integer ic
      DOUBLE PRECISION Q2,x
      DOUBLE PRECISION qu,qd,qc,g
c const
      DOUBLE PRECISION q42it,q52it,lam42,lam52
      DOUBLE PRECISION alinv,mc,PI
c local
      DOUBLE PRECISION qv,qsea,cv,cs,dcv,dcs
      DOUBLE PRECISION A0val,A1val,A2val,Bval,Cval,
     $       A0sea,B0sea,BB0sea,C0sea
      DOUBLE PRECISION A0dcv,A1dcv,A2dcv,A3dcv,Bdcv,Cdcv
      DOUBLE PRECISION Adcs, B0dcs, B1dcs, Cdcs
      DOUBLE PRECISION x1,x2,mc2q2
      DOUBLE PRECISION s,s2,s3,s4,prsccf,alstpi
      DOUBLE PRECISION WHIT6G
c parameters
      parameter(lam42=0.16d0, lam52=0.091411319d0)
      parameter(Q42IT=4.0d0, Q52IT=100.0d0)
      parameter(alinv=137.036d0, mc=1.5d0)
      parameter(pi=3.14159265358979323846d0)
      common /scale/ s,s2,s3,s4,prsccf
c
c begin
      x=ZX
      Q2=ZQ*ZQ
      ic=1
c
      x1=1.0d0-x
      x2=x**2
      mc2q2=mc**2/Q2
c
      if(Q2.lt.100.0d0) then
c  under 100 GeV^2
c
c  set  scale s
         if(Q2.lt.4.0d0) then
cccc  for under 4GeV^2 prescription
            s=  0.0d0
            prsccf =  log(Q2/LAM42)/ log(Q42IT/LAM42)
            alstpi = 6.0d0/25.0d0/ log(Q42IT/LAM42)
         else
            s=   log(  log(Q2/LAM42)/ log(Q42IT/LAM42))
            prsccf = 1.0d0
            alstpi = 6.0d0/25.0d0/ log(Q2/LAM42)
         endif
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT6 quark (U100)
c
      A0val= 2.540000d+00+s*( 2.000000d+00)+s2*( 7.180000d-01)
      A1val= 6.230000d-02+s*(-7.010000d+00)+s2*( 1.251000d-01)
      A2val=-1.642000d-01+s*(-4.360000d-01)+s2*( 1.048000d+01)
     $           +s3*(-5.200000d+00)
      Bval = 6.990000d-01+s*(-2.796000d-02)+s2*(-3.650000d-03)
      Cval = 4.420000d-01+s*(-1.255000d+00)+s2*( 1.941000d+00)
     $           +s3*(-9.950000d-01)
      A0sea= 3.180000d+00+s*( 8.690000d+00)+s2*(-2.287000d+01)
     $           +s3*( 1.896000d+01)+s4*(-5.140000d+00)
      B0sea=-1.003000d-01+s*( 1.603000d-01)+s2*(-1.037000d+00)
     $           +s3*( 9.440000d-01)+s4*(-2.915000d-01)
      BB0sea=5.690000d+00+s*( 1.867000d+01)+s2*(-4.670000d+01)
     $           +s3*( 5.050000d+01)+s4*(-1.835000d+01)
      C0sea= 2.149000d+01+s*(-5.650000d+01)+s2*( 1.293000d+02)
     $           +s3*(-1.459000d+02)+s4*( 5.750000d+01)
c
         qv  = prsccf/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= prsccf/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
            call WHIT6Q(x,mc2q2,cv,cs)
            qc = cv/alinv/2.0d0/PI + cs*alstpi
            qc  = qc*x
            ZCB=qc
         else
            qc = 0.0d0
            ZCB=qc
         endif
c
         g   = WHIT6G(x,Q2)
         g   = g*x
         ZGL=g
c
      else
c over 100 GeV^2
c
c set scale s
         s=   log(  log(Q2/LAM52)/ log(Q52IT/LAM52))
         prsccf = 1.0d0
         alstpi = 6.0d0/23.0d0/ log(Q2/LAM52)
            s2=s**2
            s3=s2*s
            s4=s2**2
c
cccccc   WHIT6 quark (O100)
c
      A0val= 4.270000d+00+s*( 3.096000d+00)+s2*( 1.621000d+00)
      A1val=-4.740000d+00+s*(-6.900000d+00)+s2*(-2.439000d+00)
      A2val= 2.837000d+00+s*( 6.460000d+00)+s2*( 4.100000d+00)
      Bval = 6.780000d-01+s*(-3.940000d-02)+s2*( 1.758000d-02)
      Cval = 1.728000d-01+s*(-2.493000d-02)+s2*( 1.451000d-01)
      A0sea= 3.340000d+00+s*(-5.610000d+00)+s2*( 5.000000d+01)
     $           +s3*(-2.207000d+02)+s4*( 3.028000d+02)
      B0sea=-2.402000d-01+s*(-4.090000d-01)+s2*( 2.263000d+00)
     $           +s3*(-1.050000d+01)+s4*( 1.487000d+01)
      BB0sea=8.790000d+00+s*(-8.860000d+00)+s2*( 1.640000d+02)
     $           +s3*(-7.120000d+02)+s4*( 9.730000d+02)
      C0sea= 9.160000d+00+s*( 9.290000d+00)+s2*(-2.784000d+02)
     $           +s3*( 1.175000d+03)+s4*(-1.592000d+03)
c
         qv  = 1.0d0/alinv/x*
     $         (A0val+A1val*x+A2val*x2) * x**Bval * x1**Cval
         qsea= 1.0d0/alinv/x*
     $         A0sea * x**(B0sea+BB0sea*x) * x1**C0sea
c
         qu  = qv/3.0d0  + qsea/6.0d0
         qu  = qu*x
         ZUV=qu
         ZUB=qu
         qd  = qv/12.0d0 + qsea/6.0d0
         qd  = qd*x
         ZDV=qd
         ZDB=qd
         ZSB=qd
         g   = WHIT6G(x,Q2)
         g   = g*x
         ZGL=g
c
         if((ic.ne.0) .and. (x*(1.0d0+4.0d0*mc2q2).lt.1.0d0)) then
      A0dcv=              s*( 1.219000d-01)+s2*( 6.200000d+00)
     $           +s3*(-2.504000d+01)+s4*( 3.098000d+01)
      A1dcv=              s*( 1.913000d+00)+s2*(-7.690000d+01)
     $           +s3*( 3.180000d+02)+s4*(-3.920000d+02)
      A2dcv=              s*(-7.160000d+00)+s2*( 2.503000d+02)
     $           +s3*(-1.062000d+03)+s4*( 1.308000d+03)
      A3dcv=              s*( 3.190000d+00)+s2*(-2.301000d+02)
     $           +s3*( 1.012000d+03)+s4*(-1.250000d+03)
      Bdcv = 4.990000d-01+s*( 3.470000d+00)+s2*(-1.526000d+01)
     $           +s3*( 1.967000d+01)
      Cdcv = 3.290000d-01+s*( 8.240000d+00)+s2*(-3.800000d+01)
     $           +s3*( 4.630000d+01)
      Adcs =              s*(-4.990000d-02)+s2*( 1.026000d-01)
     $           +s3*(-7.870000d-02)
      B0dcs=-3.610000d-01+s*(-5.760000d-01)+s2*( 2.257000d-01)
      B1dcs= 7.680000d+00+s*(-8.830000d+00)+s2*( 3.880000d+00)
      Cdcs = 2.548000d+00+s*( 6.910000d-01)+s2*(-8.700000d+00)
     $           +s3*( 1.065000d+01)
c
         dcv = 1.0d0/alinv/x*
     $         (A0dcv+x*A1dcv+x2*A2dcv+x2*x*A3dcv) * x**Bdcv * x1**Cdcv
         dcs = 1.0d0/alinv/x*
     $         Adcs * x**(B0dcs+B1dcs*x) * x1**Cdcs
c
           call WHIT6Q(x,mc*mc/Q2,cv,cs)
           qc = cv/alinv/2.0d0/PI + cs*alstpi + dcs + dcv
           qc  = qc*x
           ZCB=qc
         else
           qc = 0.0d0
           ZCB=qc
         endif
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION function WHIT1G(x,Q2)
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
cccccc   WHIT1 gluon (U100)
c
      A0g = 2.000000d+00+s*(-3.280000d+00)+s2*( 2.894000d+00)
     $          +s3*(-1.561000d+00)+s4*( 8.180000d-01)
      B0g =              s*(-7.610000d-01)+s2*(-4.900000d-02)
     $          +s3*( 4.460000d-01)
      C0g = 3.000000d+00+s*( 1.586000d+00)+s2*(-9.490000d-01)
     $          +s3*( 2.425000d+00)
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
cccccc   WHIT1 gluon (O100)
c
      A0g = 7.840000d-01+s*(-2.238000d+00)+s2*( 1.617000d+01)
     $          +s3*(-6.250000d+01)+s4*( 8.390000d+01)
      B0g =-4.030000d-01+s*(-1.307000d+00)+s2*( 8.780000d+00)
     $          +s3*(-3.580000d+01)+s4*( 5.350000d+01)
      C0g = 4.450000d+00+s*( 1.027000d+00)+s2*( 4.460000d+01)
     $          +s3*(-1.600000d+02)+s4*( 1.816000d+02)
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
         WHIT1G = prsccf/alinv/x*
     $            ( A0g * x**B0g * x1**C0g
     $             +(A1g+AA1g*x) * x**B1g * x1**C1g )
c
      return
      end
c
cccccccccccccccccccccccccccccc
c   QPM calculation
      subroutine WHIT1Q(x,r,cv,cs)
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
      DOUBLE PRECISION sum, y,z,beta,w,WHIT1G,L
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
         sum= sum + w * WHIT1G(y,mc**2/r)* L*XG(i)**(L-1.0d0)*XW(i)
c
   10 continue
c
         cs = 0.5d0/x * (x1-rx) * sum
c
      return
      end
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
c
cccccccccccccccccccccccccccccc
c   QPM calculation
      subroutine WHIT2Q(x,r,cv,cs)
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
      DOUBLE PRECISION sum, y,z,beta,w,WHIT2G,L
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
         sum= sum + w * WHIT2G(y,mc**2/r)* L*XG(i)**(L-1.0d0)*XW(i)
c
   10 continue
c
         cs = 0.5d0/x * (x1-rx) * sum
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION function WHIT3G(x,Q2)
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
cccccc   WHIT3 gluon (U100)
c
      A0g = 8.000000d+00+s*(-2.864000d+01)+s2*( 5.590000d+01)
     $          +s3*(-5.760000d+01)+s4*( 2.366000d+01)
      B0g =              s*(-9.870000d-01)+s2*( 5.100000d-01)
     $          +s3*(-6.670000d-02)
      C0g = 1.500000d+01+s*( 3.310000d-01)+s2*( 3.500000d+00)
     $          +s3*( 8.920000d-01)
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
cccccc   WHIT3 gluon (O100)
c
      A0g = 1.270000d+00+s*(-2.817000d+00)+s2*( 5.740000d+00)
     $          +s3*(-1.327000d+01)+s4*( 1.268000d+01)
      B0g =-4.610000d-01+s*(-8.170000d-01)+s2*( 3.320000d+00)
     $          +s3*(-1.296000d+01)+s4*( 1.893000d+01)
      C0g = 1.721000d+01+s*( 1.257000d+00)+s2*( 5.050000d+01)
     $          +s3*(-2.761000d+02)+s4*( 4.900000d+02)
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
         WHIT3G = prsccf/alinv/x*
     $            ( A0g * x**B0g * x1**C0g
     $             +(A1g+AA1g*x) * x**B1g * x1**C1g )
c
      return
      end
c
cccccccccccccccccccccccccccccc
c   QPM calculation
      subroutine WHIT3Q(x,r,cv,cs)
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
      DOUBLE PRECISION sum, y,z,beta,w,WHIT3G,L
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
         sum= sum + w * WHIT3G(y,mc**2/r)* L*XG(i)**(L-1.0d0)*XW(i)
c
   10 continue
c
         cs = 0.5d0/x * (x1-rx) * sum
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION function WHIT4G(x,Q2)
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
cccccc   WHIT4 gluon (U100)
c
      A0g = 4.000000d+00+s*(-9.400000d+00)+s2*( 1.555000d+01)
     $          +s3*(-1.450000d+01)+s4*( 5.470000d+00)
      B0g =              s*(-1.142000d+00)+s2*( 1.034000d+00)
     $          +s3*(-4.410000d-01)
      C0g = 3.000000d+00+s*( 8.720000d-01)+s2*( 1.006000d+00)
     $          +s3*( 3.560000d-01)
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
cccccc   WHIT4 gluon (O100)
c
      A0g = 1.384000d+00+s*(-2.455000d+00)+s2*( 8.940000d+00)
     $          +s3*(-2.906000d+01)+s4*( 3.710000d+01)
      B0g =-4.420000d-01+s*(-7.190000d-01)+s2*( 2.961000d+00)
     $          +s3*(-1.209000d+01)+s4*( 1.916000d+01)
      C0g = 4.210000d+00+s*( 2.524000d+00)+s2*( 1.003000d+01)
     $          +s3*(-1.827000d+01)+s4*( 2.162000d+00)
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
         WHIT4G = prsccf/alinv/x*
     $            ( A0g * x**B0g * x1**C0g
     $             +(A1g+AA1g*x) * x**B1g * x1**C1g )
c
      return
      end
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION function WHIT5G(x,Q2)
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
cccccc   WHIT5 gluon (U100)
c
      A0g = 1.000000d+01+s*(-3.400000d+01)+s2*( 6.900000d+01)
     $          +s3*(-7.530000d+01)+s4*( 3.230000d+01)
      B0g =              s*(-1.126000d+00)+s2*( 9.260000d-01)
     $          +s3*(-3.930000d-01)
      C0g = 9.000000d+00+s*( 4.810000d-01)+s2*( 3.200000d+00)
     $          +s3*(-3.470000d-01)
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
cccccc   WHIT5 gluon (O100)
c
      A0g = 1.995000d+00+s*(-3.260000d+00)+s2*( 1.818000d+00)
     $          +s3*( 1.711000d+00)+s4*(-4.990000d+00)
      B0g =-4.660000d-01+s*(-6.100000d-01)+s2*( 1.691000d+00)
     $          +s3*(-6.680000d+00)+s4*( 1.019000d+01)
      C0g = 1.075000d+01+s*( 5.420000d+00)+s2*( 6.550000d+00)
     $          +s3*(-2.297000d+01)+s4*( 1.867000d+01)
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
         WHIT5G = prsccf/alinv/x*
     $            ( A0g * x**B0g * x1**C0g
     $             +(A1g+AA1g*x) * x**B1g * x1**C1g )
c
      return
      end
c
cccccccccccccccccccccccccccccc
c   QPM calculation
      subroutine WHIT5Q(x,r,cv,cs)
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
      DOUBLE PRECISION sum, y,z,beta,w,WHIT5G,L
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
         sum= sum + w * WHIT5G(y,mc**2/r)* L*XG(i)**(L-1.0d0)*XW(i)
c
   10 continue
c
         cs = 0.5d0/x * (x1-rx) * sum
c
      return
      end
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
c
cccccccccccccccccccccccccccccc
c   QPM calculation
      subroutine WHIT6Q(x,r,cv,cs)
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
      DOUBLE PRECISION sum, y,z,beta,w,WHIT6G,L
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
         sum= sum + w * WHIT6G(y,mc**2/r)* L*XG(i)**(L-1.0d0)*XW(i)
c
   10 continue
c
         cs = 0.5d0/x * (x1-rx) * sum
c
      return
      end
