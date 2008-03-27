      subroutine SMRSPevolve(xin,qin,pdf)
      include 'parmsetup.inc'
      PARAMETER(NX=50)
      PARAMETER(NQ=19)
      real*8 xin,qin,pdf(-6:6),xval(45),qcdl4,qcdl5 
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      real*8 SIG,QNS,GL
      real*8 holdit
      real*8 f
      common /SMRSP/ F(7,NX,NQ,3)
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer nset
      save 
      
      call getnset(iset)
      call getnmem(iset,imem)
      
      iimem = imem
      if(iimem.eq.0) iimem = 2
      if(iimem.le.3) then
       call SMRSPxx(iimem,xin,qin,upv,dnv,usea,str,chm,bot,glu)
      endif  
      

      pdf(-6)= 0.0d0
      pdf(6)= 0.0d0
      pdf(-5)= bot
      pdf(5 )= bot
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry SMRSPread(nset)
      read(1,*)nmem(nset),ndef(nset)
      do j=1,3
      do k=1,NX
      do l=1,NQ
         read(1,*)(F(m,k,l,j),m=1,5),F(7,k,l,j),F(6,k,l,j)
      enddo
      enddo
      enddo
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry SMRSPalfa(alfas,qalfa)
        call getnset(iset)
	call getnmem(iset,imem)
	call GetOrderAsM(iset,iord)
        call Getlam4M(iset,imem,qcdl4)
        call Getlam5M(iset,imem,qcdl5)
        call aspdflib(alfas,Qalfa,iord,qcdl5)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry SMRSPinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry SMRSPpdf(mem)
c      imem = mem
      call getnset(iset)
      call setnmem(iset,mem)
      return
c
 1000 format(5e13.5)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      SUBROUTINE SMRSPxx(iset,X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C
C ::::::::::::  PION STRUCTURE FUNCTION :: 10% SEA :::::::::::::::::
C
      implicit real*8 (a-h,o-z)
      PARAMETER(NX=50)
      PARAMETER(NQ=19)
      PARAMETER(NTENTH=21)
      common /SMRSP/ F(7,NX,NQ,3)
      DIMENSION G(7),XX(NX),N0(7)
c      DIMENSION F(7,NX,NQ),G(7),XX(NX),N0(7)
      DATA XX/1.D-5,2.D-5,4.D-5,6.D-5,8.D-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.85D0,.9D0,.95D0,.975D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/0,0,3,5,0,5,0/
      DATA ZEROD/0.D0/,PONED/0.1D0/,ONED/1.D0/,ONEDO/1.1D0/,TWOD/2.D0/
      DATA INIT/0/
      XSAVE=X
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,NX-1
      DO 20 M=1,NQ
         DO 25 I=1,7
  25     F(I,N,M,iset)=F(I,N,M,iset)/(ONED-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,NTENTH-1
      XX(J)= LOG10(XX(J))+ONEDO
      DO 30 I=2,6
      DO 30 K=1,NQ
  30  F(I,J,K,iset)= LOG(F(I,J,K,iset))
     +              *F(I,NTENTH,K,iset)/ LOG(F(I,NTENTH,K,iset))
  50  FORMAT(7F10.5)
      DO 40 I=1,7
      DO 40 M=1,NQ
  40  F(I,NX,M,iset)=ZEROD
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.PONED) XXX= LOG10(X)+ONEDO
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM= LOG(QSQ/QSQMIN)/ LOG(TWOD)
      B=RM-AINT(RM)
      M=1+  INT(RM)
      DO 60 I=1,7
      G(I)= (ONED-A)*(ONED-B)*F(I,N,M,iset)+(ONED-A)*B*F(I,N,M+1,iset)
     .    + A*(ONED-B)*F(I,N+1,M,iset)  + A*B*F(I,N+1,M+1,iset)
      IF(N.GE.NTENTH) GOTO 65
      IF(I.EQ.7.OR.I.EQ.1) GOTO 65
          FAC=(ONED-B)*F(I,NTENTH,M,iset)+B*F(I,NTENTH,M+1,iset)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(ONED-X)**N0(I)
  60  CONTINUE
C UPBAR DISTRIBUTION = D DISTRIBUTION
      UPV=G(2)
      DNV=G(2)
C THIS SEA IS (UBAR+DBAR)/2
      SEA=G(4)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
      X=XSAVE
      RETURN
      END
