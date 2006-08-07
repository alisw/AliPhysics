      subroutine DOGevolvep0(xin,qin,p2in,ip2in,pdf)
      include 'parmsetup.inc'
      real*8 xin,qin,q2in,p2in,pdf(-6:6),xval(45),qcdl4,qcdl5
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer ns
      
      save 

      call DOPHO1(xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu)

      pdf(-6)= 0.0d0
      pdf(6)= 0.0d0
      pdf(-5)= bot
      pdf(5 )= bot
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
      entry DOGevolvep1(xin,qin,p2in,ip2in,pdf)
      
      call DOPHO2(xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu)

      pdf(-6)= 0.0d0
      pdf(6)= 0.0d0
      pdf(-5)= bot
      pdf(5 )= bot
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
      entry DOGread(nset)
      read(1,*)nmem(nset),ndef(nset)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry DOGalfa(alfas,qalfa)
        call getnset(iset)
	call GetOrderAsM(iset,iord)
        call Getlam4M(iset,imem,qcdl4)
        call Getlam5M(iset,imem,qcdl5)
        call aspdflib(alfas,Qalfa,iord,qcdl5)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry DOGinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry DOGpdf(mem)
      imem = mem
      return
c
 1000 format(5e13.5)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       SUBROUTINE DOPHO1(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL)
C********************************************************************
C*                                                                  *
C*    Parametrization of parton distribution functions              *
C*    in the photon (LO analysis) - asymptotic solution of AP eq.!  *
C*                                                                  *
C* authors:  D.Duke and H.Owens (DO)                                *
C*           /Phys.Rev. D26 (1982) 1600/                            *
C*                                                                  *
C* Prepared by:                                                     *
C*             Krzysztof Charchula, DESY                            *
C*             bitnet: F1PCHA@DHHDESY3                              *
C*             decnet: 13313::CHARCHULA                             *
C*                                                                  *
C* Modified by:                                                     *
C*             H. Plothow-Besch/CERN-PPE                            *
C*                                                                  *
C********************************************************************
C
      implicit real*8 (a-h,o-z)
       double precision
     +        CQ(5),
     +        DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL
      PARAMETER (ALPEM=7.29927D-3, PI=3.141592D0)
      PARAMETER (ALAM=0.2D0)
      DATA CQ/0.33333D0,0.66666D0,0.33333D0,0.66666D0,0.33333D0/
C
       Q2 = DQ*DQ
       ALAM2=ALAM**2
       FQ=ALPEM/(2.*PI)*LOG(Q2/ALAM2)
C
C...gluons
       POMG=0.194*(1.-DX)**1.03/(DX**0.97)
       DGL=POMG*FQ
C
C...quarks
        POM1=(1.81-1.67*DX+2.16*DX**2)
        POM2=DX**0.7/(1.-0.4*LOG(1.-DX))
        POM3=38.D-4*(1.-DX)**1.82/(DX**1.18)
          DDB=(CQ(1)**2*POM1*POM2+POM3)*FQ
          DDV=DDB
          DUB=(CQ(2)**2*POM1*POM2+POM3)*FQ
          DUV=DUB
          DSB=(CQ(3)**2*POM1*POM2+POM3)*FQ
          DCB=(CQ(4)**2*POM1*POM2+POM3)*FQ
          DBB=(CQ(5)**2*POM1*POM2+POM3)*FQ
       RETURN
       END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       SUBROUTINE DOPHO2(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL)
C********************************************************************
C*                                                                  *
C*    Parametrization of parton distribution functions              *
C*    in the photon (LO analysis) - asymptotic solution of AP eq.!  *
C*                                                                  *
C* authors:  D.Duke and H.Owens (DO)                                *
C*           /Phys.Rev. D26 (1982) 1600/                            *
C*                                                                  *
C* Prepared by:                                                     *
C*             Krzysztof Charchula, DESY                            *
C*             bitnet: F1PCHA@DHHDESY3                              *
C*             decnet: 13313::CHARCHULA                             *
C*                                                                  *
C* Modified by:                                                     *
C*             H. Plothow-Besch/CERN-PPE                            *
C*                                                                  *
C********************************************************************
C
      implicit real*8 (a-h,o-z)
      double precision
     +        CQ(5),
     +        DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL
      PARAMETER (ALPEM=7.29927D-3,PI=3.141592D0)
      PARAMETER (ALAM=0.4D0)
      DATA CQ/0.33333D0,0.66666D0,0.33333D0,0.66666D0,0.33333D0/
C
       Q2 = DQ*DQ
       ALAM2=ALAM**2
       FQ=ALPEM/(2.*PI)*LOG(Q2/ALAM2)
C
C...gluons
       POMG=0.194*(1.-DX)**1.03/(DX**0.97)
       DGL=POMG*FQ
C
C...quarks
        POM1=(1.81-1.67*DX+2.16*DX**2)
        POM2=DX**0.7/(1.-0.4*LOG(1.-DX))
        POM3=38.D-4*(1.-DX)**1.82/(DX**1.18)
          DDB=(CQ(1)**2*POM1*POM2+POM3)*FQ
          DDV=DDB
          DUB=(CQ(2)**2*POM1*POM2+POM3)*FQ
          DUV=DUB
          DSB=(CQ(3)**2*POM1*POM2+POM3)*FQ
          DCB=(CQ(4)**2*POM1*POM2+POM3)*FQ
          DBB=(CQ(5)**2*POM1*POM2+POM3)*FQ
       RETURN
       END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
