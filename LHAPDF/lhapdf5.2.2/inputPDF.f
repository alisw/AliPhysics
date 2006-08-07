      
      subroutine weightPDF(x)
      implicit none
      integer nset,imem,nfmax
      real*8 x,Q
      nset=1
      call weightPDFM(nset,x)
      return

      entry GetNF(nfmax)
      nset=1
      call GetNFM(nset,nfmax)
      return

      entry GetThreshold(imem,Q)
      nset=1
      call GetThresholdM(nset,imem,Q)
      return
      
      end
      
      subroutine parmPDF(nset,x,pdf)
      implicit none
      include 'parmsetup.inc'
      integer nset
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mem
      common/NAME/name,nmem,ndef,mem
      real*4 xp(40),x4,fgdis4
      character*16 s1,s2
      integer i,j,imem,nop,parmN,Nfunc(nmxset),Fw(nmxset),nfmax,M
      real*8 x,N,b0,Poly,pdf(-6:6),Fparm(nopmax),F(nofmax)
      real*8 Ccoef(nmxset,-6:6,nofmax),Fpow(nmxset,nofmax)
      real*8 Q,Treshold(nmxset,-6:6)
c      data Treshold/39*0d0/
      integer Fmap(nmxset,nofmax,npfmax)
      integer Ftype(nmxset,nofmax),Fn(nmxset,nofmax),Ctype(nmxset,-6:6)
      integer lhasilent
      common/lhasilent/lhasilent
      logical first
      data first/.true./
      save Nfunc,Fn,Fw,Fpow,Fmap,Ccoef,Fparm,Ftype,Ctype,first,Treshold
            
      do i=1,Nfunc(nset)
         if (Ftype(nset,i).eq.1) then
            Poly=1.0

            do j=4,Fn(nset,i)
         Poly=Poly+Fparm(Fmap(nset,i,j))*x**(float(j-3)/Fpow(nset,i))
            enddo
            Poly=Fparm(Fmap(nset,i,1))*Poly
       F(i)=x**Fparm(Fmap(nset,i,2))*(1.0-x)**Fparm(Fmap(nset,i,3))*Poly
         endif
         if (Ftype(nset,i).eq.2) then
            if (x.lt.0.9999999) then
      Poly=Fparm(Fmap(nset,i,2))*log(x)+Fparm(Fmap(nset,i,3))*log(1.0-x)
     .      +Fparm(Fmap(nset,i,4))*x
     .      +Fparm(Fmap(nset,i,6))*log(1.0+x*exp(Fparm(Fmap(nset,i,5))))
         F(i)=Fparm(Fmap(nset,i,1))*exp(Poly)
            else
               F(i)=0d0
            endif
         endif
         if (Ftype(nset,i).eq.101) then
       Poly=exp(Fparm(Fmap(nset,i,1)))*x**(Fparm(Fmap(nset,i,2))-1)
     .       *(1d0-x)**Fparm(Fmap(nset,i,3))
       Poly=
     . Poly+(1d0+Fparm(Fmap(nset,i,4))*x)*(1d0-x)**Fparm(Fmap(nset,i,5))
            b0=10d0
            if (Poly.gt.b0) then
               F(i)=Poly
            elseif (Poly.lt.-b0) then
               F(i)=0d0
            else
               F(i)=Poly+log(1d0+exp(-b0*Poly)-exp(-b0))/b0
            endif
         endif
c - to add the mrst2004 gluon convolution
         if (Ftype(nset,i).eq.201) then
	   xp(2) = Fparm(Fmap(nset,i,1))
	   xp(3) = Fparm(Fmap(nset,i,2))
	   xp(23)= Fparm(Fmap(nset,i,3))
	   xp(16)= Fparm(Fmap(nset,i,4))
	   xp(5) = Fparm(Fmap(nset,i,5))
	   xp(40)= Fparm(Fmap(nset,i,6))
	   xp(24)= Fparm(Fmap(nset,i,7)) 
	   xp(20)= Fparm(Fmap(nset,i,8)) 
           xp(4) = Fparm(Fmap(nset,i,9)) 
	   xp(1) = Fparm(Fmap(nset,i,10))
	   x4 = x
           call gconv(x4,xp,fgdis4)
          F(i) = FGDIS4
	 endif
c--
      enddo
      do i=-6,6
         pdf(i)=0.0
         if (Ctype(nset,i).gt.0) then
            if (Ctype(nset,i).eq.1) then
               do j=1,Nfunc(nset)
                  pdf(i)=pdf(i)+Ccoef(nset,i,j)*F(j)
c		  print *,i,j,Ccoef(i,j),F(j)
               enddo
            endif
            if (Ctype(nset,i).eq.101) then
               if (i.eq.-2) then
       pdf(i)=F(int(Ccoef(nset,i,1)))/(F(int(Ccoef(nset,i,2)))+1d0)
               endif
               if (i.eq.-1) then
      pdf(i)=
     .F(int(Ccoef(nset,i,1)))/(1d0/F(int(Ccoef(nset,i,2)))+1d0)
               endif
               if (i.eq.1) then
                  pdf(i)=F(1)+pdf(-1)
               endif
               if (i.eq.2) then
                  pdf(i)=F(2)+pdf(-2)
               endif
            endif
         endif
      enddo
c      print *,pdf
      return
*
      entry weightPDFM(nset,x)
      if (Fw(nset).ge.0) then
         x=Fparm(Fw(nset))
      else
         call numberPDF(nop)
         x=1.0/float(nop)
      endif
      return
*
c      entry GetParmPDFM(nset,imem,x)
      entry GetParmPDF(nset,imem,x)
      x=Fparm(imem)
      return
*
      entry GetNfM(nset,nfmax)
      nfmax=0
      do i=1,6
         if (Treshold(nset,-i).ge.0d0) nfmax=nfmax+1
         if (Treshold(nset,i).ge.0d0) nfmax=nfmax+1
      enddo
      nfmax=nfmax/2
      return
*
      entry GetThresholdM(nset,imem,Q)
      Q=Treshold(nset,imem)
      return
*
c      entry InitEvolvePDFM(nset,imem)
      entry InitEvolvePDF(nset,imem)
c      print*, 'calling listPDF',nset,imem
      call listPDF(nset,imem,Fparm)
c      print *,Fparm
      return
*
c      entry initInputPDFM(nset)
      entry initInputPDF(nset)

      if(first) then
        do i=1,nmxset
	  do j=-6,6
             Treshold(i,j)=0.0d0
	  enddo
	enddo
        first=.false.
      endif
      read(1,*) s1,Fw(nset),Nfunc(nset)
c      print *,s1,Fw,Nfunc
      if(lhasilent.eq.0) then
        write(*,*) 'Parametrization: ',s1
        write(*,*)
      endif
      do i=1,Nfunc(nset)
         Ftype(nset,i)=-1
         read(1,*) s1,s2
         if (index(s2,'x-taylor').eq.1) then
            Ftype(nset,i)=1
            read(1,*) FPow(nset,i),Fn(nset,i)
            read(1,*) (Fmap(nset,i,j),j=1,Fn(nset,i))
         endif
         if (index(s2,'log-pade').eq.1) then
            Ftype(nset,i)=2
            FPow(nset,i)=0d0
            read(1,*) Fn(nset,i)
            read(1,*) (Fmap(nset,i,j),j=1,Fn(nset,i))
         endif
         if (index(s2,'cteq6-ratio').eq.1) then
            Ftype(nset,i)=101
            Fpow(nset,i)=0d0
            Fn(nset,i)=5
            read(1,*) (Fmap(nset,i,j),j=1,Fn(nset,i))
         endif
         if (index(s2,'convol').eq.1) then
            Ftype(nset,i)=201
            read(1,*) Fn(nset,i)
            read(1,*) (Fmap(nset,i,j),j=1,Fn(nset,i))
         endif
         if (Ftype(nset,i).lt.0) then
            write(*,*) 'File description error:'
            write(*,*) 'Unknown functional ',s2
            stop
         endif
      enddo
      read(1,*) s1
      do i=-6,6
         Ctype(nset,i)=-1
         read(1,*) s1,s2
c	 print *,s1,s2
         if (index(s2,'none').eq.1) then
            Ctype(nset,i)=0
            Treshold(nset,i)=-1d0
         endif
         if (index(s2,'treshold').eq.1) then
            Ctype(nset,i)=0
            read(1,*) Treshold(nset,i)
         endif
         if (index(s2,'composite').eq.1) then
            Ctype(nset,i)=1
            Treshold(nset,i)=0d0
            read(1,*) (Ccoef(nset,i,j),j=1,Nfunc(nset))
         endif
         if (index(s2,'cteq6-ratio').eq.1) then
            Ctype(nset,i)=101
            Treshold(nset,i)=0d0
            read(1,*) (Ccoef(nset,i,j),j=1,3)
         endif
         if (Ctype(nset,i).lt.0) then
            write(*,*) 'File description error:'
            write(*,*) 'Unknown composit type ',s2
            stop
         endif
      enddo
      if (Fw(nset).ge.0) then
         write(*,*) '***********************************************'
         write(*,*) '* Note that this is a weigthed PDF set.       *'
         write(*,*) '* See manual for proper use.                  *'
         write(*,*) '***********************************************'
      endif
      return
*
      end

      FUNCTION BETA_LHA(X1,X2)
      CALL GAMMA_LHA(X1,G1,IER)
      IF(IER.NE.0) write(16,*) 'GAMMA_LHA ERROR: IER= ',IER,X1,X2
      CALL GAMMA_LHA(X2,G2,IER)
      IF(IER.NE.0) write(16,*) 'GAMMA_LHA ERROR: IER= ',IER,X1,X2
      X3=X1+X2
      CALL GAMMA_LHA(X3,G3,IER)
      IF(IER.NE.0) write(16,*) 'GAMMA_LHA ERROR: IER= ',IER,X1,X2
      BETA_LHA=G1*G2/G3
      RETURN
      END

      SUBROUTINE GAMMA_LHA(XX,GX,IER)
      IF(XX-34.5)6,6,4
    4 IER=2
      GX=1.E38
      RETURN
    6 X=XX
      ERR=1.0E-6
      IER=0
      GX=1.0
      IF(X-2.0)50,50,15
   10 IF(X-2.0)110,110,15
   15 X=X-1.0
      GX=GX*X
      GO TO 10
   50 IF(X-1.0)60,120,110
C        SEE IF X IS NEAR NEGATIVE INTEGER OR ZERO
   60 IF(X-ERR)62,62,80
   62 K=X
      Y=FLOAT(K)-X
      IF(ABS(Y)-ERR)130,130,64
   64 IF(1.0-Y-ERR)130,130,70
C        X NOT NEAR A NEGATIVE INTEGER OR ZERO
   70 IF(X-1.0)80,80,110
   80 GX=GX/X
      X=X+1.0
      GO TO 70
  110 Y=X-1.0
      GY=1.0+Y*(-0.5771017+Y*(+0.9858540+Y*(-0.8764218+Y*(+0.8328212+
     1Y*(-0.5684729+Y*(+0.2548205+Y*(-0.05149930)))))))
      GX=GX*GY
  120 RETURN
  130 IER=1
      RETURN
      END

      FUNCTION	ALPHA(T,AL)
      COMMON/AINPUT/IORD,QSCT,QSDT
      COMMON/PARAM/PARA(40)
      DATA PI/3.14159/
      DATA TOL/.0005/
      ITH=0
      TT=T
      qsctt=qsct/4.
      qsdtt=qsdt/4.
c      AL=para(1)
      AL2=AL*AL
      FLAV=4.
      QS=AL2*EXP(T)

      if(qs.lt.0.5d0) then   !!  running stops below 0.5
          qs=0.5d0
          t=alog(qs/al2)
          tt=t
      endif

      IF(QS.gt.QSCTT) GO	TO 12  
      IF(QS.lt.QSDTT) GO	TO 312  
   11 CONTINUE
      B0=11-2.*FLAV/3. 
      IF(IORD)1,1,2
c     IF(IORD)2,2,2	!TAKE CARE !!
    1 CONTINUE
      ALPHA=4.*PI/B0/T
      RETURN
    2 CONTINUE
      X1=4.*PI/B0
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS=X1/T*(1.-X2*aLOG(T)/T)
    5 CONTINUE
      F=-T+X1/AS-X2*aLOG(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF(DEL-TOL)3,3,4
    3 CONTINUE
      ALPHA=AS2
      IF(ITH.EQ.0) RETURN
      GO TO (13,14,15) ITH
    4 CONTINUE
      AS=AS2
      GO TO 5
   12 ITH=1
      T=aLOG(QSCTT/AL2)
      GO TO 11
   13 ALFQC4=ALPHA
      FLAV=5.   
      ITH=2
      GO TO 11
   14 ALFQC5=ALPHA
      ITH=3
      T=TT
      GO TO 11
   15 ALFQS5=ALPHA
      ALFINV=1./ALFQS5+1./ALFQC4-1./ALFQC5
      ALPHA=1./ALFINV
      RETURN

  311 CONTINUE
      B0=11-2.*FLAV/3. 
      IF(IORD)31,31,32
c     IF(IORD)32,32,32	!TAKE CARE !!
   31 CONTINUE
      ALPHA=4.*PI/B0/T
      RETURN
   32 CONTINUE
      X1=4.*PI/B0
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS=X1/T*(1.-X2*aLOG(T)/T)
   35 CONTINUE
      F=-T+X1/AS-X2*aLOG(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF(DEL-TOL)33,33,34
   33 CONTINUE
      ALPHA=AS2
      IF(ITH.EQ.0) RETURN
      GO TO (313,314,315) ITH
   34 CONTINUE
      AS=AS2
      GO TO 35
  312 ITH=1
      T=aLOG(QSDTT/AL2)
      GO TO 311
  313 ALFQC4=ALPHA
      FLAV=3.   
      ITH=2
      GO TO 311
  314 ALFQC3=ALPHA
      ITH=3
      T=TT
      GO TO 311
  315 ALFQS3=ALPHA
      ALFINV=1./ALFQS3+1./ALFQC4-1./ALFQC3
      ALPHA=1./ALFINV
      RETURN
      END

      SUBROUTINE WATE96
C*******************************************************************
C*****							       *****
C***** THE X(I)	AND W(I) ARE THE DIRECT	OUTPUT FROM A PROGRAM  *****
C***** USING NAG ROUTINE D01BCF	TO CALCULATE THE	       *****
C***** GAUSS-LEGENDRE WEIGHTS FOR 96 POINT INTEGRATION.	       *****
C***** THEY AGREE TO TYPICALLY 14 DECIMAL PLACES WITH THE      *****
C***** TABLE IN	ABRAMOWITZ & STEGUN, PAGE 919.		       *****
C*****							       *****
C***** ---->   PETER HARRIMAN, APRIL 3RD 1990.		       *****
C*****							       *****
C*******************************************************************
      DIMENSION	X(48),W(48)
      COMMON/GAUS96/XI(96),WI(96),nterms,XX(97)
      NTERMS=96

      X( 1)=   0.01627674484960183561
      X( 2)=   0.04881298513604856015
      X( 3)=   0.08129749546442434360
      X( 4)=   0.11369585011066471632
      X( 5)=   0.14597371465489567682
      X( 6)=   0.17809688236761733390
      X( 7)=   0.21003131046056591064
      X( 8)=   0.24174315616383866556
      X( 9)=   0.27319881259104774468
      X(10)=   0.30436494435449495954
      X(11)=   0.33520852289262397655
      X(12)=   0.36569686147231213885
      X(13)=   0.39579764982890709712
      X(14)=   0.42547898840729897474
      X(15)=   0.45470942216774136446
      X(16)=   0.48345797392059470382
      X(17)=   0.51169417715466604391
      X(18)=   0.53938810832435567233
      X(19)=   0.56651041856139533470
      X(20)=   0.59303236477757022282
      X(21)=   0.61892584012546672523
      X(22)=   0.64416340378496526886
      X(23)=   0.66871831004391424358
      X(24)=   0.69256453664216964528
      X(25)=   0.71567681234896561582
      X(26)=   0.73803064374439816819
      X(27)=   0.75960234117664555964
      X(28)=   0.78036904386743123629
      X(29)=   0.80030874413913884180
      X(30)=   0.81940031073792957139
      X(31)=   0.83762351122818502758
      X(32)=   0.85495903343459936363
      X(33)=   0.87138850590929436968
      X(34)=   0.88689451740241818933
      X(35)=   0.90146063531585023110
      X(36)=   0.91507142312089592706
      X(37)=   0.92771245672230655266
      X(38)=   0.93937033975275308073
      X(39)=   0.95003271778443564022
      X(40)=   0.95968829144874048809
      X(41)=   0.96832682846326217918
      X(42)=   0.97593917458513455843
      X(43)=   0.98251726356301274934
      X(44)=   0.98805412632962202890
      X(45)=   0.99254390032376081654
      X(46)=   0.99598184298720747465
      X(47)=   0.99836437586317963722
      X(48)=   0.99968950388322870559
      W( 1)=   0.03255061449236316962
      W( 2)=   0.03251611871386883307
      W( 3)=   0.03244716371406427668
      W( 4)=   0.03234382256857594104
      W( 5)=   0.03220620479403026124
      W( 6)=   0.03203445623199267876
      W( 7)=   0.03182875889441101874
      W( 8)=   0.03158933077072719007
      W( 9)=   0.03131642559686137819
      W(10)=   0.03101033258631386231
      W(11)=   0.03067137612366917839
      W(12)=   0.03029991542082762553
      W(13)=   0.02989634413632842385
      W(14)=   0.02946108995816795100
      W(15)=   0.02899461415055528410
      W(16)=   0.02849741106508543861
      W(17)=   0.02797000761684838950
      W(18)=   0.02741296272602931385
      W(19)=   0.02682686672559184485
      W(20)=   0.02621234073567250055
      W(21)=   0.02557003600534944960
      W(22)=   0.02490063322248370695
      W(23)=   0.02420484179236479915
      W(24)=   0.02348339908592633665
      W(25)=   0.02273706965832950717
      W(26)=   0.02196664443874448477
      W(27)=   0.02117293989219144572
      W(28)=   0.02035679715433347898
      W(29)=   0.01951908114014518992
      W(30)=   0.01866067962741165898
      W(31)=   0.01778250231604547316
      W(32)=   0.01688547986424539715
      W(33)=   0.01597056290256253144
      W(34)=   0.01503872102699521608
      W(35)=   0.01409094177231515264
      W(36)=   0.01312822956696188190
      W(37)=   0.01215160467108866759
      W(38)=   0.01116210209983888144
      W(39)=   0.01016077053500880978
      W(40)=   0.00914867123078384552
      W(41)=   0.00812687692569928101
      W(42)=   0.00709647079115442616
      W(43)=   0.00605854550423662775
      W(44)=   0.00501420274292825661
      W(45)=   0.00396455433844564804
      W(46)=   0.00291073181793626202
      W(47)=   0.00185396078894924657
      W(48)=   0.00079679206555731759
      DO 1 I=1,48
      XI(I)=-X(49-I)
      WI(I)=W(49-I)
      XI(I+48)=X(I)
      WI(I+48)=W(I)
    1 CONTINUE
      DO 2 I=1,96
    2 XX(I)=0.5*(XI(I)+1.)
      XX(97)=1.0
      EXPON=1.0
      DO 3 I=1,96
      YI=2.*(0.5*(1.+XI(I)))**EXPON-1.
      WI(I)=WI(I)/(1.+XI(I))*(1.+YI)*EXPON
      XI(I)=YI
      XX(I)=0.5*(1.+YI)
    3 CONTINUE
      RETURN
      END
      
      subroutine gconv(x,xp,fgdis)
      COMMON/AINPUT/IORD,QSCT,QSDT
      common/GAUS96/XI(96),WI(96),NTERMS,XX(97)
      dimension xp(40)
      logical first
      data first/.true./
      if(first) then
         call wate96
	 first=.false.
      endif
      PI = 3.14159
      PI2 = PI*PI
      iord = 1
      qsdt=8.18    !!  This is the value of 4m_c^2
      qsct=74.0    !!  This is the value of 4m_b^2
      cf = 4./3.
      eta4 = xp(40)
      T=alog(1/xp(1)**2)
c      AL = 0.550/(4.* 3.14159)
      AL=ALPHA(T,xp(1))/(4.* pi)
      rx=sqrt(x)
      FF1=BETA_LHA(XP(2),XP(3)+1.)+XP(16)*BETA_LHA(XP(2)+1.,XP(3)+1.)+
     XXP(23)*BETA_LHA(XP(2)+0.5,XP(3)+1.)
      FF2=BETA_LHA(XP(2)+1.,XP(3)+1.)+
     XXP(16)*BETA_LHA(XP(2)+2.,XP(3)+1.)+
     XXP(23)*BETA_LHA(XP(2)+1.5,XP(3)+1.)
      FF3=BETA_LHA(XP(5),ETA4+1.)+XP(20)*BETA_LHA(XP(5)+1.,ETA4+1.)+
     XXP(24)*BETA_LHA(XP(5)+0.5,ETA4+1.)
      FF4=BETA_LHA(XP(5)+1.,ETA4+1.)+XP(20)*BETA_LHA(XP(5)+2.,ETA4+1.)+
     XXP(24)*BETA_LHA(XP(5)+1.5,ETA4+1.)
      COEFU=2.*XP(4)/FF1
      COEFD=   XP(4)/FF3
c	   print *,'coefu ',coefu
c	   print *,'coefd ',coefd
      UV=coefu*X**XP(2)*(1.-X)**XP(3)*(1.+XP(16)*X+XP(23)*SQRT(X))
      DV=coefd*X**XP(5)*(1.-X)**ETA4*(1.+XP(20)*X+XP(24)*SQRT(X))
      FGDIS=al*CF*(-9.-2.*PI2/3.+alog(1.-x)*(-3.+
     .2.*alog(1.-x)))*(UV+DV)
     
      DO 23 M=1,NTERMS
      Y=0.5*(1.-X)*XI(M)+0.5*(1.+X)
      XY=X/Y
      UVXY=coefu*XY**XP(2)*(1.-XY)**XP(3)*(1.+XP(16)*XY
     .+XP(23)*SQRT(XY))
      DVXY=coefd*XY**XP(5)*(1.-XY)**ETA4*(1.+XP(20)*XY
     .+XP(24)*SQRT(XY))
      AL1=ALOG(1.-Y)
      C22=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*ALOG(Y)-2.*(1.+Y)*ALOG(1.-Y))
      C23=CF*(-3.+4.*ALOG(1.-Y))/(1.-Y)

      FGDIS=FGDIS+.5*(1.-X)*WI(M)*al*
     .(C22*(uvxy+dvxy)+C23*(uvxy+dvxy-uv-dv))

   23 CONTINUE
      
      return
      end
