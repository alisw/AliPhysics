      subroutine GSGevolvep0(xin,qin,p2in,ip2in,pdf)
      include 'parmsetup.inc'
      real*8 xin,qin,q2in,p2in,pdf(-6:6),xval(45),qcdl4,qcdl5
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      real*8 SIG,QNS,GL
      real*8 holdit
      common/gsgdat/SIG(78,11,3),QNS(78,11,3),GL(78,11,3)
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer nset
      save 
      
      call getnset(iset)
      call getnmem(iset,iimem)
c -- this is LO 2 --> 3 0/1 --> 2  
      if(iimem.eq.2) iimem = 3
      if(iimem.eq.0) iimem = 2
      if(iimem.eq.1) iimem = 2

      call SFGSHL(iimem,xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu)

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
      entry GSGevolvep1(xin,qin,p2in,ip2in,pdf)
      
c--- this is HO --- iimem=1
      iimem = 1
       call SFGSHL(iimem,xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu)
      

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
      entry GSGread(nset)
      read(1,*)nmem(nset),ndef(nset)
      do j=1,3
      do k=1,78
      do m=1,11
         read(1,*)SIG(k,m,j),QNS(k,m,j),GL(k,m,j)
      enddo
      enddo
      enddo
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GSGalfa(alfas,qalfa)
        call getnset(iset)
	call getnmem(iset,imem)
	call GetOrderAsM(iset,iord)
        call Getlam4M(iset,imem,qcdl4)
        call Getlam5M(iset,imem,qcdl5)
        call aspdflib(alfas,Qalfa,iord,qcdl5)
c        call aspdflib(alfas,Qalfa,iord,qcdl5)

      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GSGinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GSGpdf(mem)
      call getnset(iset)
      call setnmem(iset,mem)
c      imem = mem
      return
c
 1000 format(5e13.5)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      SUBROUTINE SFGSHO(X,Q,U,D,US,DS,S,C,B,G)
      SUBROUTINE SFGSHL(iset,X,Q,U,D,US,DS,S,C,B,G)
C
*****************************************************************
* Subroutine returns the parton distributions in the photon in  *
* higher  order. u,d etc. gives the actual distributions and    *
* not x times the distributions; Q2 means Q 2. The distributions*
* are valid for 5.0e -4< x < 1.0 and 5.3 GeV 2 < Q 2 < 1.0e 8.  *
* if higher Q 2 or lower x is required, these may be obtained   *
* from the authors on request.                                  *
* Lionel Gordon July 1991 : Gordon@uk.ac.man.ph.v2              *
*****************************************************************
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NP=78,NQ=11,NARG=2)
      double precision
     +       DBFINT,
c     +       SIG(NP,NQ),QNS(NP,NQ),GL(NP,NQ),Y(NP),
     +       Y(NP),
     +       XT(NARG),A(NP+NQ),QT(NQ)
      common/gsgdat/SIG(78,11,3),QNS(78,11,3),GL(78,11,3)
      DIMENSION NA(NARG)
      EXTERNAL GSXCOR
c      SAVE SIG,QNS,GL,Y,ICALL
      SAVE Y,ICALL
      DATA QT /5.3D0,20.0D0,50.0D0,1.0D2,5.0D2,1.0D3,1.0D4,1.0D5,
     * 1.0D6,1.0D7,1.0D8/
      DATA ZEROD/0.D0/
      DATA ICALL/0/
******************************************************************
       U = ZEROD
       D = ZEROD
       S = ZEROD
       C = ZEROD
       B = ZEROD
       G = ZEROD
*if x is out of range
      if(iset.eq.1) then
        IF((X.LT.5.0D-4).OR.(X.GT.0.95D0)) GOTO 90
      else
        IF((X.LT.5.0D-4).OR.(X.GT.0.99D0)) GOTO 90
      endif
      
*******************************************************************
       IF (ICALL.NE.1) THEN
* get the x coordinates
          CALL GSXCOR(Y,NP)
          ICALL=1
        END IF
*
      DO 30 IX=1,NP
        A(IX)=Y(IX)
   30 CONTINUE
      DO 40 IQ=1,NQ
        A(NP+IQ)=QT(IQ)
   40 CONTINUE
*
      Q2 = Q*Q
      NA(1)=NP
      NA(2)=NQ
       XT(1)=X
       XT(2)=Q2
      XSIG=DBFINT(2,XT,NA,A,SIG(1,1,iset))
      XQNS=DBFINT(2,XT,NA,A,QNS(1,1,iset))
        G =DBFINT(2,XT,NA,A,GL(1,1,iset))
*
      IF (Q2.LT.50.0D0) THEN
C Use three flavour evolution.
       U=(XSIG+9.0D0*XQNS)/6.0D0
       D=(XSIG-4.5D0*XQNS)/6.0D0
       S=D
       C=ZEROD
       B=ZEROD
*
      ELSE IF((Q2.GT.50.0D0).AND.(Q2.LT.250.0D0)) THEN
C Use four flavour evolution
      U=(XSIG+6.0D0*XQNS)/8.0D0
      D=(XSIG-6.0D0*XQNS)/8.0D0
      S=D
      C=U
      B=ZEROD
      ELSE
C Use five flavour evolution
      U=(XSIG+7.5D0*XQNS)/10.0D0
      D=(XSIG-5.0D0*XQNS)/10.0D0
      S=D
      C=U
      B=D
      ENDIF
      U=X*U
      US=U
      D=X*D
      DS=D
      S=X*S
      C=X*C
      B=X*B
      G=X*G
 90   RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE GSXCOR(Y,NP)
C
      implicit real*8 (a-h,o-z)
      double precision
     +       Y(NP)
      N=1
      DO 10 IX=1,20,2
         Y(N)=    (IX)/2000.0D0
      N=N+1
   10 CONTINUE
      DO 20 IX=30,200,10
         Y(N)=    (IX)/2000.0D0
           N=N+1
   20 CONTINUE
      DO 30 IX=240,1600,40
         Y(N)=    (IX)/2000.0D0
      N=N+1
   30 CONTINUE
      DO 40 IX=1625,1980,25
         Y(N)=    (IX)/2000.0D0
      N=N+1
   40 CONTINUE
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
