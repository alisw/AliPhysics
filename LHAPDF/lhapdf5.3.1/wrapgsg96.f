      subroutine GSG96evolvep0(xin,qin,p2in,ip2in,pdf)
      include 'parmsetup.inc'
      real*8 xin,qin,q2in,p2in,pdf(-6:6),xval(45),qcdl4,qcdl5
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      real*8 SIG,QNS,GL
      common/gsgdat/SIG(78,11,3),QNS(78,11,3),GL(78,11,3)
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer nset
     
      save 
      
      iimem = 2
      call GS96HL(iimem,xin,qin,upv,dnv,str,chm,bot,glu)
      usea = upv
      dsea = dnv

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
      entry GSG96evolvep1(xin,qin,p2in,ip2in,pdf)

      iimem = 1
      call GS96HL(iimem,xin,qin,upv,dnv,str,chm,bot,glu)
      usea = upv
      dsea = dnv

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
      entry GSG96read(nset)
      read(1,*)nmem(nset),ndef(nset)
      
      do m = 1,2
      do k = 1,78
        read(1,*)(sig(k,j,m),j=1,4)
        read(1,*)(sig(k,j,m),j=5,8)
        read(1,*)(sig(k,j,m),j=9,11)
      enddo
      do k = 1,78
        read(1,*)(qns(k,j,m),j=1,4)
        read(1,*)(qns(k,j,m),j=5,8)
        read(1,*)(qns(k,j,m),j=9,11)
      enddo
      do k = 1,78
        read(1,*)(gl(k,j,m),j=1,4)
        read(1,*)(gl(k,j,m),j=5,8)
        read(1,*)(gl(k,j,m),j=9,11)
      enddo
      enddo
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GSG96alfa(alfas,qalfa)
       call getnset(iset)
       call getnmem(iset,imem)
       call GetOrderAsM(iset,iord)
       call Getlam4M(iset,imem,qcdl4)
       call Getlam5M(iset,imem,qcdl5)
       call aspdflib(alfas,Qalfa,iord,qcdl5)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GSG96init(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GSG96pdf(mem)
      call getnset(iset)
      call setnmem(iset,mem)
c      imem = mem
      return
c
 1000 format(5e13.5)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      SUBROUTINE GS96HO(X,Q,U,D,S,C,B,G)
      SUBROUTINE GS96HL(iset,X,Q,U,D,S,C,B,G)
      implicit real*8 (a-h,o-z)
      PARAMETER(NP=78,NQ=11,NARG=2)
*****************************************************************
* Subroutine returns the parton distributions in the photon in  *
* next-to-leading order. u,d etc. gives the actual distributions*
* not x times the distributions; Q2 means Q^2. The distributions*
* are valid for 5.0e^-4< x < 1.0 and 5.3 GeV^2 < Q^2 < 1.0e^8.  *
* if higher Q^2 or lower x is required, these may be obtained   *
* from the authors on request.                                  * 
* Lionel Gordon April 1996 : Gordon@hep.anl.gov                 *
* John Storrow             :johns@a3.ph.man.ac.uk               *
*****************************************************************
c      DIMENSION SIG(NP,NQ),QNS(NP,NQ),GL(NP,NQ),Y(NP)
      DIMENSION Y(NP)
      DIMENSION XT(NARG),NA(NARG),A(NP+NQ),QT(NQ)
      common/gsgdat/SIG(78,11,3),QNS(78,11,3),GL(78,11,3)
      EXTERNAL GS2XCOR
c      SAVE SIG,QNS,GL,Y,ICALL
      SAVE Y,ICALL
      DATA QT /3.0D0,20.0D0,50.0D0,100.0D0,500.0D0,1.0D3,1.0D4,1.0D5,
     + 1.0D6,1.0D7,1.0D8/
******************************************************************
        q2=q*q
*if x is out of range
       IF((X.LT.5.0D-4).OR.(X.GT.1.0D0)) GOTO 90
*******************************************************************
       IF (ICALL.NE.1) THEN
* get the x coordinates
       CALL GS2XCOR(Y,NP)
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
       C=0.0D0
       B=0.0D0
*      
      ELSE IF((Q2.GT.50.0).AND.(Q2.LT.250.0)) THEN
C Use four flavour evolution 
      U=(XSIG+6.0D0*XQNS)/8.0D0
      D=(XSIG-6.0D0*XQNS)/8.0D0
      S=D
      C=U
      B=0.0D0
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
      SUBROUTINE GS2XCOR(Y,NP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(NP)
      N=1
      DO 10 IX=1,20,2
         Y(N)=DBLE(IX)/2000.0D0
      N=N+1
   10 CONTINUE
      DO 20 IX=30,200,10
         Y(N)=DBLE(IX)/2000.0D0
      N=N+1
   20 CONTINUE
      DO 30 IX=240,1600,40
         Y(N)=REAL(IX)/2000.0D0
      N=N+1
   30 CONTINUE
      DO 40 IX=1625,1980,25
         Y(N)=REAL(IX)/2000.0D0
      N=N+1
   40 CONTINUE
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
