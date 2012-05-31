      subroutine MRST98evolve(x,Q,pdf)
      implicit real*8(a-h,o-z)
      include 'parmsetup.inc'
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      integer nset
      common/NAME/name,nmem,ndef,mmem
      parameter(nx=49,nq=37,ntenth=23,np=8,members=5)
      real*8 pdf(-6:6)
      real*8 f(0:members,np,nx,nq+1)
      real*8 qq(nq),xx(nx),xxin(nx),g(np),n0(np)
      data xxin/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      save
c
      xsave=x
      qsq = q*q
      q2save=qsq
c
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
c
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(imem,i,n,m)+(1d0-a)*b*f(imem,i,n,m+1)
     .	  +       a*(1d0-b)*f(imem,i,n+1,m)+a*b*f(imem,i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(imem,i,ntenth,m)+b*f(imem,i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
c
      pdf(0)  = glu
      pdf(1)  = dnv+dsea
      pdf(-1) = dsea
      pdf(2)  = upv+usea
      pdf(-2) = usea
      pdf(3)  = str
      pdf(-3) = str
      pdf(4)  = chm
      pdf(-4) = chm
      pdf(5)  = bot
      pdf(-5) = bot
      pdf(6)  = 0.0d0
      pdf(-6) = 0.0d0
      
      x=xsave
      qsq=q2save
      return
*
      entry MRST98read(nset)
      read(1,*)nmem(nset),ndef(nset)
c - first resotre the xx array
      do j=1,nx
        xx(j)=xxin(j)
      enddo
c - next read in the data points
      do nm = 0,nmem(nset)
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(nm,1,n,m),f(nm,2,n,m),f(nm,3,n,m),f(nm,4,n,m),
     .		  f(nm,5,n,m),f(nm,7,n,m),f(nm,6,n,m),f(nm,8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(nm,i,n,m)=f(nm,i,n,m)/(1d0-xx(n))**n0(i)
  20    continue
c        write(*,*)'PDF set ',nm,' first element ',f(nm,1,1,1)
        do 31 j=1,ntenth-1
c        xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
        do 31 i=1,8
        if(i.eq.5.or.i.eq.7) goto 31
        do 30 k=1,nq
  30    f(nm,i,j,k)=dlog10(f(nm,i,j,k)/f(nm,i,ntenth,k))
     .              +f(nm,i,ntenth,k)
  31    continue
  50    format(8f10.5)
        do 40 i=1,np
        do 40 m=1,nq
  40  f(nm,i,nx,m)=0d0
       enddo
      do 32 j=1,ntenth-1
        xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
  32  continue
      return
*
      entry MRST98alfa(alfas,Qalfa)
        call alphamrs(5,alfas,Qalfa)
      return
*
      entry MRST98init(Eorder,Q2fit)
      return
*
      entry MRST98pdf(mem)
c      if(mem.eq.0) mem=ndef
      imem = mem
c      print *,imem

      return
c
      end

