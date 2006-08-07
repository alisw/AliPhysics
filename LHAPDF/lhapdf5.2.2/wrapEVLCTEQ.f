      subroutine EVLCTEQevolve(x,Q,f)
      implicit none
      include 'parmsetup.inc'
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      real*8 f(-6:6)
      real*8 x,Q,xx,qq,up,dn,ubar,dbar,gluon,sbar,cbar,bbar,tbar,Blam
      real*8 xmin,qmax,anx,anq,Q2fit,qini,Ahdn,ordi,anf,mass,pi,alfas
      real*8 Q0,alfas0
      real*8 CtLhpardis,CtLhALPI,parmflavor,CtLhastolam
      integer nflavor,Eorder,j,nx,nq
      integer nset
      data pi / 3.141592653589793d0 /
      external parmflavor
      save nflavor,xmin,qmax,anx,anq

*
	xx = x
	qq = Q
	up   = xx*CtLhpardis(1,xx,qq)
	dn   = xx*CtLhpardis(2,xx,qq) 
	ubar = xx*CtLhpardis(-1,xx,qq)
	dbar = xx*CtLhpardis(-2,xx,qq)
	gluon= xx*CtLhpardis( 0,xx,qq)

	if(nflavor .ge. 3) then
	   sbar = xx*CtLhpardis(-3,xx,qq)
	else
	   sbar = 0.d0
	endif

	if(nflavor .ge. 4) then
	   cbar = xx*CtLhpardis(-4,xx,qq)
	else
	   cbar = 0.d0
	endif

	if(nflavor .ge. 5) then
	   bbar = xx*CtLhpardis(-5,xx,qq)
	else
	   bbar = 0.d0
	endif

	if(nflavor .eq. 6) then
	   tbar = xx*CtLhpardis(-6,xx,qq)
	else
	   tbar = 0.d0
	endif

	f(6) = tbar 
	f(5) = bbar 
	f(4) = cbar 
	f(3) = sbar 
	f(2) = up 
	f(1) = dn 
	f(0) = gluon 
	f(-1) = dbar 
	f(-2) = ubar 
	f(-3) = sbar 
	f(-4) = cbar 
	f(-5) = bbar 
	f(-6) = tbar 
*
	return
*
      entry EVLCTEQread(nset)
      call CtLhbldat1
      call CtLhbldat2
      read(1,*) xmin,qmax,nx,nq
      anx=nx
      anq=nq-1
      return
*
      entry EVLCTEQalfa(alfas,Q)
      alfas = pi*CtLhALPI(Q)
      return
*
      entry EVLCTEQinit(nset,Eorder,Q2fit)
*
      Ahdn=1d0
      call CtLhParPdf(1,'IHDN',Ahdn,j)
      call CtLhParPdf(1,'QMAX',qmax,j)
      call CtLhParPdf(1,'XMIN',xmin,j)
      call GetOrderAs(j)
      ordi=j+1.0
      call CtLhParQcd(1,'ORDR',ordi,j)
      call GetThresholdM(nset,4,mass)
      call CtLhParQcd(1,'M4',mass,j)
      call GetThresholdM(nset,5,mass)
      call CtLhParQcd(1,'M5',mass,j)
      mass=180d0
      call CtLhParQcd(1,'M6',mass,j)
      call GetNfM(nset,j)
      nflavor=j
      Anf=nflavor
      call CtLhParQcd(1,'NFL',Anf,j)
      qini=sqrt(Q2fit)
      call CtLhParPdf(1,'QINI',qini,j)
      call CtLhEvlpar(1,'NFMX', Anf, j)
      ordi=Eorder+1.0
      call CtLhParPdf(1,'IKNL',ordi,j)
      return
*
      entry EVLCTEQpdf(nset)
      call GetOrderAsM(nset,j)
      call GetAlfas(nset,alfas0,Q0)
      Blam=CtLhastolam(alfas0,Q0,j+1,nflavor)
      call CtLhSetLam (nflavor,Blam,j+1)
      call CtLhParPdf(1,'NX',anx,j)
      call CtLhParPdf(1,'NT',anq,j)
      call CtLhEvolve (parmflavor,j)
      if (j .ne. 0) then
         write(*,*) 'EVLCTEQ Evolve Error code :',j
         stop
      endif
      return
*
      end
*
      double precision function parmflavor(i,x)
      implicit none
      real*8 x,f(-6:6)
      integer i,i0
      integer iset
*
      call getnset(iset)
      call parmPDF(iset,x,f)
      i0=i
      if (i.eq.-2) i0=-1
      if (i.eq.-1) i0=-2
      if (i.eq.1) i0=2
      if (i.eq.2) i0=1
      parmflavor=f(i0)/x
      return
*
      end
      
      function CtLhastolam(as,q,nloop,nf) 
      implicit double precision (a-h,o-z)
      data pi / 3.141592653589793d0 /
      xlp = nloop-1d0
      b  = (33.0-2.0*nf)/pi/12.0 
      bp = (153.0 - 19.0*nf) / pi / 2.0 / (33.0 - 2.0*nf) * xlp 
      t  = 1.0/b/as 
c----------------------------------------------------------- 
c Solve the equation 
c 
    1 xlt = log(t) 
      ot = t 
c----------------------------------------------------------- 
c Solve the equation 
c Value and Derivative of alfa with respect to t 
c 
      as0  = 1/b/t - bp*xlt/(b*t)**2 
      as1  = - 1/b/t**2 -bp/b**2*(1-2*xlt)/t**3 
      t  = (as-as0)/as1 + t 
      if(abs(ot-t)/ot.gt..00001)goto 1 
      CtLhastolam = q/exp(t/2.0) 
      return 
      end 

