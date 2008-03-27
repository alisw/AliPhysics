      subroutine OWPevolve(xin,qin,pdf)
      include 'parmsetup.inc'
      real*8 xin,qin,pdf(-6:6),xval(45),qcdl4,qcdl5    
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset)
      common/NAME/name,nmem,ndef,mmem
      integer mem,mmem
      integer nset
c      integer iset,iimem
c      common/SET/iset,iimem
      
      save 
      
      q2in = qin*qin
c      iset = imem
       
      if(imem.eq.0) then
         call strowp1(xin,Qin,upv,dnv,usea,str,chm,glu)
      elseif(imem.eq.1) then
         call strowp1(xin,Qin,upv,dnv,usea,str,chm,glu)
      elseif(imem.eq.2) then
         call strowp2(xin,Qin,upv,dnv,usea,str,chm,glu)
      else
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
      pdf(2 )= upv+usea
      pdf(-1)= usea
      pdf(1 )= dnv+usea
      pdf(0 )= glu

      return
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry OWPread(nset)
      read(1,*)nmem(nset),ndef(nset)
c      iset = nset
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry OWPalfa(alfas,qalfa)
        call getnset(iset)
	call GetOrderAsM(iset,iord)
c	print*,'from getorderasm',iord
        call Getlam4M(iset,imem,qcdl4)
c	print*,'from getorderasm',iord
        call Getlam5M(iset,imem,qcdl5)
c	print*,'from getorderasm',iord
        call aspdflib(alfas,Qalfa,iord,qcdl5)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry OWPinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry OWPpdf(mem)
      imem = mem
      return
c
 1000 format(5e13.5)
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C*********************************************************************

      SUBROUTINE STROWP1(X,SCALE,UPV,DNV,SEA,STR,CHM,GL)
C ::::::::::::::  OWENS SET 1 PION STRUCTURE FUNCTION  :::::::::::::::
      implicit real*8 (a-h,o-z)
      DOUBLE PRECISION DGAMMA_LHA
      double precision
     +       COW(3,5,4),TS(6),XQ(9)

C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1),IS=1,5),IP=1,3)/
     1  4.0000D-01,  7.0000D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -6.2120D-02,  6.4780D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     3 -7.1090D-03,  1.3350D-02,  0.0000D+00,  0.0000D+00,  0.0000D+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2),IS=1,5),IP=1,3)/
     1  8.8800D-01,  0.0000D+00,  3.1100D+00,  6.0000D+00,  0.0000D+00,
     2 -1.8020D+00, -1.5760D+00, -1.3170D-01,  2.8010D+00, -1.7280D+01,
     3  1.8120D+00,  1.2000D+00,  5.0680D-01, -1.2160D+01,  2.0490D+01/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3),IS=1,5),IP=1,3)/
     1  9.0000D-01,  0.0000D+00,  5.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -2.4280D-01, -2.1200D-01,  8.6730D-01,  1.2660D+00,  2.3820D+00,
     3  1.3860D-01,  3.6710D-03,  4.7470D-02, -2.2150D+00,  3.4820D-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4),IS=1,5),IP=1,3)/
     1  0.0000D+00, -2.2120D-02,  2.8940D+00,  0.0000D+00,  0.0000D+00,
     2  7.9280D-02, -3.7850D-01,  9.4330D+00,  5.2480D+00,  8.3880D+00,
     3 -6.1340D-02, -1.0880D-01, -1.0852D+01, -7.1870D+00, -1.1610D+01/

       DATA ZEROD/0.D0/, ONED/1.D0/, SIXD/6.D0/
       DATA ALAM/0.2D0/, Q02/4.D0/, QMAX2/2.D3/
C...Pion structure functions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.

C...Determine set, Lambda and s expansion variable.
        Q2 = SCALE*SCALE
        Q2IN = MIN( QMAX2,MAX( Q02,Q2))
        SD = LOG( LOG( Q2IN/ALAM**2)/ LOG( Q02/ALAM**2))

C...Calculate structure functions.
        DO 240 KFL=1,4
        DO 230 IS=1,5
  230   TS(IS)=COW(1,IS,KFL)+COW(2,IS,KFL)*SD+
     &  COW(3,IS,KFL)*SD*SD
        IF(KFL.EQ.1) THEN
cif defined(CERNLIB_SINGLE)
c          DENOM = GAMMA(TS(1))*GAMMA(TS(2)+ONED)/GAMMA(TS(1)+TS(2)+ONED)
cendif
cif defined(CERNLIB_DOUBLE)
          DENOM = DGAMMA_LHA(TS(1))*DGAMMA_LHA(TS(2)+ONED)/
     +                              DGAMMA_LHA(TS(1)+TS(2)+ONED)
cendif
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/DENOM
        ELSE
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2)
        ENDIF
  240   CONTINUE

C...Put into output arrays.
        UPV = XQ(1)
        DNV = XQ(1)
        SEA = XQ(3)/SIXD
        STR = XQ(3)/SIXD
        CHM = XQ(4)
        BOT = ZEROD
        TOP = ZEROD
        GL  = XQ(2)
C
        RETURN
        END

C*********************************************************************

      SUBROUTINE STROWP2(X,SCALE,UPV,DNV,SEA,STR,CHM,GL)
C ::::::::::::::  OWENS SET 2 PION STRUCTURE FUNCTION  :::::::::::::::
      implicit real*8 (a-h,o-z)
      DOUBLE PRECISION DGAMMA_LHA
      double precision
     +       COW(3,5,4),TS(6),XQ(9)

C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1),IS=1,5),IP=1,3)/
     1  4.0000D-01,  6.2800D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -5.9090D-02,  6.4360D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     3 -6.5240D-03,  1.4510D-02,  0.0000D+00,  0.0000D+00,  0.0000D+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2),IS=1,5),IP=1,3)/
     1  7.9400D-01,  0.0000D+00,  2.8900D+00,  6.0000D+00,  0.0000D+00,
     2 -9.1440D-01, -1.2370D+00,  5.9660D-01, -3.6710D+00, -8.1910D+00,
     3  5.9660D-01,  6.5820D-01, -2.5500D-01, -2.3040D+00,  7.7580D+00/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3),IS=1,5),IP=1,3)/
     1  9.0000D-01,  0.0000D+00,  5.0000D+00,  0.0000D+00,  0.0000D+00,
     2 -1.4170D-01, -1.6970D-01, -2.4740D+00, -2.5340D+00,  5.6210D-01,
     3 -1.7400D-01, -9.6230D-02,  1.5750D+00,  1.3780D+00, -2.7010D-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4),IS=1,5),IP=1,3)/
     1  0.0000D+00, -8.8200D-02,  1.9240D+00,  0.0000D+00,  0.0000D+00,
     2  6.2290D-02, -2.8920D-01,  2.4240D-01, -4.4630D+00, -8.3670D-01,
     3 -4.0990D-02, -1.0820D-01,  2.0360D+00,  5.2090D+00, -4.8400D-02/

       DATA ZEROD/0.D0/, ONED/1.D0/, SIXD/6.D0/
       DATA ALAM/0.4D0/, Q02/4.D0/, QMAX2/2.D3/
C...Pion structure functions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.

C...Determine set, Lambda and s expansion variable.
        Q2 = SCALE*SCALE
        Q2IN = MIN( QMAX2,MAX( Q02,Q2))
        SD = LOG( LOG( Q2IN/ALAM**2)/ LOG( Q02/ALAM**2))

C...Calculate structure functions.
        DO 10 KFL=1,4
        DO 20 IS=1,5
   20   TS(IS)=COW(1,IS,KFL)+COW(2,IS,KFL)*SD+
     &  COW(3,IS,KFL)*SD*SD
        IF(KFL.EQ.1) THEN
cif defined(CERNLIB_SINGLE)
c         DENOM = GAMMA(TS(1))*GAMMA(TS(2)+ONED)/GAMMA(TS(1)+TS(2)+ONED)
cendif
cif defined(CERNLIB_DOUBLE)
          DENOM = DGAMMA_LHA(TS(1))*DGAMMA_LHA(TS(2)+ONED)/
     +                              DGAMMA_LHA(TS(1)+TS(2)+ONED)
cendif
          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/DENOM
        ELSE
          XQ(KFL)=TS(1)*X**TS(2)*(1.-X)**TS(3)*(1.+TS(4)*X+TS(5)*X**2)
        ENDIF
   10   CONTINUE

C...output
        UPV = XQ(1)
        DNV = XQ(1)
        SEA = XQ(3)/SIXD
        STR = XQ(3)/SIXD
        CHM = XQ(4)
        GL  = XQ(2)
C
        RETURN
        END
c**************************************************************************
