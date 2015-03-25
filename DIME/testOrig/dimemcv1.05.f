      program DimeMC

      implicit double precision(a-y)       
      implicit complex*16(z)
      double precision pt1(2),pt2(2),ptx(2),q(4,20)
      double precision pboo(4),pcm(4),plb(4)
      double precision tvec(4),uvec(4),svec(4),vv1(4),vv2(4)
      integer d,h,i,j,k,l,m,n,o,p,r,iset,af,nhist,ntotal,eflag,ichi,
     &icut,ncut,id(20),id0,id1,id2,pdf,num,ptgen,mode,iord,ii,ip,ll,nev
     &,hh,lmax
      character prefix*50,fsp*10,order*10,pflag*10,fsi*10,formf*10
     &,ppbar*10,output*10,mregge*10,cuts*10,unw*10
      integer iin,nch
      double precision cc0(5),bm(5),bb0(5),pp0(5),bex(5),gaa(5)
      common/pars/cc0,bm,bb0,pp0,bex,asp,sigo,gaa,ep,norm
      common/ipars/nch
      common/mom/q
      common/int/ip
      common/mompt/pt1,pt2,ptx
      common/vars/s,rts,mmes,yx
      common/cuts/etaelmax,etaelmin,ptelmin,ptphmin,ecut,rmax,rmin,mcut
      common/flags/pflag
      common/gluons/g1,g2,kt1,kt2
      common/levicivita/epsilon
      common/ang/cost,costa,costb
      common/alph/alfas
      common/hvars/sh,th,uh
      common/wvars/sig0,bb,t1,t2
      common/alphas/alphap,alpha0,alphapr,alpha0r,alphapm,alpha0m
      common/sec/cpom,cf,crho,aff,ar
      integer iii
      common/ivar/iii
      common/ff/formf
      common/ffpars/bexp,ao,bo,aoo,boo
      common/regge/mregge
ccccc hepevt output
      integer nmxhep,kk
      parameter (nmxhep=4000)
      integer nevhep,nhep,isthep,idhep,jmohep,jdahep
      double precision phep,vhep
      common /hepevt/ nevhep,nhep,isthep(nmxhep),idhep(nmxhep),
     &jmohep(2,nmxhep),jdahep(2,nmxhep),phep(5,nmxhep),vhep(4,nmxhep)
ccccc Les Houches  Event Common Block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)
ccccc rambo variables
      integer npart,nrun
      double precision pin(4),am(100),pout(4,100),wt,ein
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                         c
c     Dime MC for the central exclusive production of     c
c     meson pairs by double Pomeron exchange:             c
c                                                         c
c     p(1) p(2) --> p(3) + M_1(5) M_2(6) + p(4)           c
c                                                         c
c     Momenta for each event in array q(i,j), where j is  c
c     the particle label and i is the 4-momentum          c
c     component, with:                                    c
c                                                         c
c     1,2 = transverse components                         c
c     3   = beam component                                c 
c     4   = energy                                        c
c                                                         c
c     PDG number for ith particle in arrary ID(i)         c
c                                                         c
c     Also gives event record in HEPEVT or LHE format     c
c     (others are available upon request)                 c
c                                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                         c
c     Particles generated:                                c
c                                                         c
c     pflag='pipm'  :  M_1=pi^+ M_2=pi^-                  c
c     pflag='pi0'   :  M_1=M_2=pi^0                       c
c     pflag='kpkm'  :  M_1=K^+  M_2=K^-                   c
c     pflag='ks'    :  M_1=M_2=K_0                        c
c     pflag='rho'   :  M_1=M_2=rho_0                      c
c                                                         c
c     with decay: rho(5) --> pi^+(7)pi^-(8)               c
c                 rho(6) --> pi^+(9)pi^-(10)              c
c     according to phase space, with finite rho           c
c     width included                                      c
c                                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                         c
c     User defined cuts can readily be implemented        c
c     in subroutine 'icut'                                c
c                                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                         c
c     This is version 1.05           : 17 Sep 2014        c 
c                                                         c
c     Comments etc to: lucian.harland-lang@durham.ac.uk   c
c                                                         c
c     If you use the code please reference (and see for   c
c     details of model used):                             c
c                                                         c
c     L.A. Harland-Lang, V.A. Khoze, M.G. Ryskin          c
c     and W.J. Stirling arXiv:1105.1626                   c
c                                                         c
c     L.A. Harland-Lang, V.A. Khoze, M.G. Ryskin          c
c     and W.J. Stirling arXiv:1204.4803                   c
c                                                         c
c     L.A. Harland-Lang, V.A. Khoze and  M.G. Ryskin      c
c     arXiv:1312.4553                                     c
c                                                         c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(35,file='exerec1.dat') ! event record
      
      rts=7.0d3         ! centre of mass energy

cccc  Some basic cuts, imposed in subtroutine 'icut'. Other user defined cuts can readily be implemented in subroutine
cccc  note: in rhorho case these cuts are imposed on rho's, and not their decay productions. Cuts on the decay products
cccc  can be imposed by editing 'icut'

      rmax=2.0d0        ! max meson rapidity
      rmin=-2.0d0       ! min meson rapidity
      ecut=0.2d0        ! min meson p_t

c -- physics parameters 

      pflag='pipm'    ! Process generated - see preamble for options
      fsi='true'      ! phenomenological model for exclusive suppression in Pom Pom --> meson pair. To turn on/off --> 'true/false'
      formf='orexp'   ! meson - pomeron form factor.
      ppbar='false'   ! set true if ppbar collisions
      output='lhe'    ! Event record style, HEPEVT or LHE
      cuts='true'     ! Impose cuts or not
      unw='true'      ! Set = 'true' for unweighted events
      iin=1           ! Model for soft survival factor, as described in arXiv:1306.2149. Default = 1

cccccccc

      ntotal=100000          ! no. of runs for weighted events
      nev=100                 ! no. of unweighted events generated to event record

ccccccccc

      if(formf.eq.'exp')then           ! Set parameters for Pomeron-Meson form factor in function 'fpi(x)'
         bexp=1d0/2.2d0
      elseif(formf.eq.'orexp')then
         bo=1d0/1.1d0
         ao=dsqrt(0.5d0)
      elseif(formf.eq.'power')then
         aoo=1.7d0
      endif

ccccccccccccccccccccccccccccccccccccccc

ccccccccc   Start of main body of code

cccccccccccccccccccccccccccccccccccccc

      call initpars(iin)   ! Initialise soft survival parameters
      call calcop          ! proton opacity 
      call calcscreen      ! screening amplitude

      if(pflag.eq.'pipm')then
         mmes=0.13957018d0      ! pi+/- mass, PDG 2011 value
         sig0=13.63d0
         alphapm=0.7d0
         alpha0m=-0.7d0*mmes**2
         cf=31.79d0/0.389d0
         crho=4.23d0/0.389d0

         nup=6
         istup(5)=1
         idup(5)=211
         mothup(1,5)=1
         mothup(2,5)=2
         icolup(1,5)=0
         icolup(2,5)=0
         vtimup(5)=0
         spinup(5)=9.

         istup(6)=1
         idup(6)=-211
         mothup(1,6)=1
         mothup(2,6)=2
         icolup(1,6)=0
         icolup(2,6)=0
         vtimup(6)=0
         spinup(6)=9.
         
      elseif(pflag.eq.'pi0')then
         mmes=0.1349766d0       ! pi0 mass, PDG 2011 value
         sig0=13.63d0
         alphapm=0.7d0
         alpha0m=-0.7d0*mmes**2
         cf=31.79d0/0.389d0
         crho=0d0

         nup=6
         istup(5)=1
         idup(5)=111
         mothup(1,5)=1
         mothup(2,5)=2
         icolup(1,5)=0
         icolup(2,5)=0
         vtimup(5)=0
         spinup(5)=9.

         istup(6)=1
         idup(6)=-111
         mothup(1,6)=1
         mothup(2,6)=2
         icolup(1,6)=0
         icolup(2,6)=0
         vtimup(6)=0
         spinup(6)=9.

      elseif(pflag.eq.'kpkm')then
         mmes=0.493677d0         ! K+/- mass, PDG 2011 value
         sig0=11.82d0
         cf=17.255d0/0.389d0
         crho=9.105d0/0.389d0

         nup=6
         istup(5)=1
         idup(5)=321
         mothup(1,5)=1
         mothup(2,5)=2
         icolup(1,5)=0
         icolup(2,5)=0
         vtimup(5)=0
         spinup(5)=9.

         istup(6)=1
         idup(6)=-321
         mothup(1,6)=1
         mothup(2,6)=2
         icolup(1,6)=0
         icolup(2,6)=0
         vtimup(6)=0
         spinup(6)=9.

      elseif(pflag.eq.'ks')then
         mmes=0.497614d0         ! K_0 mass, PDG 2011 value
         sig0=11.82d0
         cf=17.255d0/0.389d0
         crho=0d0

         nup=6
         istup(5)=1
         idup(5)=310
         mothup(1,5)=1
         mothup(2,5)=2
         icolup(1,5)=0
         icolup(2,5)=0
         vtimup(5)=0
         spinup(5)=9.

         istup(6)=1
         idup(6)=310
         mothup(1,6)=1
         mothup(2,6)=2
         icolup(1,6)=0
         icolup(2,6)=0
         vtimup(6)=0
         spinup(6)=9.

      elseif(pflag.eq.'rho')then
         mmes0=0.77549d0          ! rho mass, PDG 2013 value  
         mwidth=0.1491d0         ! rho width PDG 2013 value
         sig0=10d0
         cf=0d0
         crho=0d0

         nup=10
         istup(5)=2
         idup(5)=113
         mothup(1,5)=1
         mothup(2,5)=2
         icolup(1,5)=0
         icolup(2,5)=0
         vtimup(5)=0
         spinup(5)=9.

         istup(6)=2
         idup(6)=113
         mothup(1,6)=1
         mothup(2,6)=2
         icolup(1,6)=0
         icolup(2,6)=0
         vtimup(6)=0
         spinup(6)=9

         do k=7,10
            istup(k)=1
            mothup(2,k)=0
            icolup(1,k)=0
            icolup(2,k)=0
            vtimup(k)=0
            spinup(k)=9.
         enddo

         idup(7)=211
         idup(8)=-211
         idup(9)=211
         idup(10)=-211
         mothup(1,7)=5
         mothup(1,8)=5
         mothup(1,9)=6
         mothup(1,10)=6
         
      endif

c -- other parameters

      ebeam=rts/2d0
      s=rts*rts
      zi=(0d0,1d0)
      rt2=dsqrt(2d0)
      pi=dacos(-1d0)
      bp=rts/dsqrt(2.0d0)
      mp=0.93827d0
      beta=dsqrt(1d0-4d0*mp**2/s)
      s0=1d0

c     Pomeron + t-slope

      bb=4d0
      bjac=6d0
      bjac1=2d0

      alphap=0.25d0    ! D-L '92 fit
      alpha0=1.0808d0
      alphapr=0.93d0
      alpha0r=0.55d0

      mf127=1.275d0
      mf1525=1.525d0

      cpom=sig0/0.389d0
      aff=-0.860895d0
      ar=-1.16158d0

      mmes1=mmes
      mmes2=mmes

ccccc Initialise rambo (rho0 decay)

      if(pflag.eq.'rho')then

         mmes=mmes0
         npart=2
         do j=1,4
            am(j)=0.13957018d0
         enddo
      endif

ccccccc

c     initialise counters etc

c      do hh=1,20

      nhist=1
      sum=0d0
      sum1=0d0
      ncut=0

      weightm=0d0

c     initialise histograms

      do i=1,nhist
         call histo3(i)
      enddo

      num=0

cccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Incoming protons to event record array

cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ID(1)=2212
      q(1,1)=0d0
      q(2,1)=0d0
      q(3,1)=ebeam*beta
      q(4,1)=ebeam

      istup(1)=-1
      idup(1)=2212
      mothup(1,1)=0
      mothup(2,1)=0
      icolup(1,1)=0
      icolup(2,1)=0
      do j=1,4
         pup(j,1)=q(j,1)
      enddo
      pup(5,1)=dsqrt(q(4,1)**2-q(3,1)**2-q(2,1)**2-q(1,1)**2)
      vtimup(1)=0
      spinup(1)=9.

      q(1,2)=0d0
      q(2,2)=0d0
      q(3,2)=-ebeam*beta
      q(4,2)=ebeam

      istup(2)=-1
      if(ppbar.eq.'true')then
         idup(2)=-2212
      else
         idup(2)=2212
      endif
      mothup(1,2)=0
      mothup(2,2)=0
      icolup(1,2)=0
      icolup(2,2)=0
      do j=1,4
         pup(j,2)=q(j,2)
      enddo
      pup(5,2)=dsqrt(q(4,2)**2-q(3,2)**2-q(2,2)**2-q(1,2)**2)
      vtimup(2)=0
      spinup(2)=9.

cccc  outgoing initial info

      istup(3)=1
      if(ppbar.eq.'true')then
         idup(2)=-2212
      else
         idup(2)=2212
      endif
      idup(3)=2212
      mothup(1,3)=1
      mothup(2,3)=2
      icolup(1,3)=0
      icolup(2,3)=0
      vtimup(3)=0
      spinup(3)=9

      istup(4)=1
      idup(4)=2212
      mothup(1,4)=1
      mothup(2,4)=2
      icolup(1,4)=0
      icolup(2,4)=0
      vtimup(4)=0
      spinup(4)=9

ccc   HEPEVT

      if(output.eq.'hepevt')then

      nhep=nup
      
      do k=1,5
         phep(k,1)=pup(k,1)
         phep(k,2)=pup(k,2)
      enddo

      do k=1,nup
         isthep(k)=istup(k)
         idhep(k)=idup(k)
         jmohep(1,k)=mothup(1,k)
         jmohep(2,k)=jmohep(1,k)
         jdahep(1,k)=0
         jdahep(2,k)=0
         vhep(1,k)=0d0
         vhep(2,k)=0d0
         vhep(3,k)=0d0
         vhep(4,k)=0d0
      enddo

      if(pflag.eq.'rho')then
      jdahep(1,5)=7
      jdahep(2,5)=8
      jdahep(1,6)=9
      jdahep(2,6)=10
      endif
      
      jmohep(2,5)=2
      jmohep(2,6)=2

      endif

      if(unw.eq.'true')then
         lmax=2
      else
         lmax=1
      endif

      do ll=1,lmax

      if(ll.eq.2)then
c         ntotal=nev*10
         ntotal=10000000
      endif

      ip=ntotal+1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c     Start of event loop                               c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(ll.eq.1)then
         write(6,*)'Generating weighted events...'
      else
        write(6,*)'Generating unweighted events...'
      endif

      do i=1,ntotal

      weight=0d0

      call r2455(ran0)
      call r2455(ran1)
      call r2455(ran2)
      call r2455(ran3)
      call r2455(ran4)
      call r2455(ran5)
      call r2455(ran6)

ccccccccccccccccccccccccccccccccccccccccccccccccccc

c     incoming protons

      ID(1)=2212
      q(1,1)=0d0
      q(2,1)=0d0
      q(3,1)=ebeam
      q(4,1)=ebeam

      ID(2)=2212
      q(1,2)=0d0
      q(2,2)=0d0
      q(3,2)=-ebeam
      q(4,2)=ebeam    

      phi1=2d0*pi*ran0
      phi2=2d0*pi*ran1

      pt1sq=-dlog(ran2)/(bjac)
      pt2sq=-dlog(ran3)/(bjac)

      weight=(dexp(bjac*pt1sq)*dexp(bjac*pt2sq))/bjac**2 

      pt1(1)=dsqrt(pt1sq)*dsin(phi1)
      pt1(2)=dsqrt(pt1sq)*dcos(phi1)
      pt2(1)=dsqrt(pt2sq)*dsin(phi2)
      pt2(2)=dsqrt(pt2sq)*dcos(phi2)

      pt1x=pt1(1)
      pt1y=pt1(2)
      pt2x=pt2(1)
      pt2y=pt2(2)

cccccccccccccccccccccccccccccccc

ccc non-zero rho width

cccccccccccccccccccccccccccccccc

       if(pflag.eq.'rho')then  ! new part here

 677     call r2455(rm1)
         call r2455(rm2)

         msmax=mmes0+4d0*mwidth
         msmin=2d0*0.13957018d0

cccccccccccccccc

         almin=datan(-(mmes0**2-msmin**2)/mwidth/mmes0)
         almax=datan(-(mmes0**2-msmax**2)/mwidth/mmes0)

         al1=almin+(almax-almin)*rm1
         al2=almin+(almax-almin)*rm2

         mmes1=dsqrt(dtan(al1)*mmes0*mwidth+mmes0**2)
         mmes2=dsqrt(dtan(al2)*mmes0*mwidth+mmes0**2)

         weight=weight*(almax-almin)**2
         weight=weight*mwidth**2*mmes0**2
         weight=weight*(1d0+dtan(al1)**2)*(1d0+dtan(al2)**2)
         weight=weight/4d0/mmes1/mmes2

ccccccccccccccccc

         mwidth1=mwidth*((1d0-4d0*0.13957018d0**2/mmes1**2)/
     &        (1d0-4d0*0.13957018d0**2/mmes0**2))**1.5d0
         mwidth2=mwidth*((1d0-4d0*0.13957018d0**2/mmes2**2)/
     &        (1d0-4d0*0.13957018d0**2/mmes0**2))**1.5d0

         if(mmes1.lt.2d0*0.13957018d0) goto 677
         if(mmes2.lt.2d0*0.13957018d0) goto 677

         weight=weight*2d0*mmes0*mmes1*mwidth1/pi
         weight=weight*2d0*mmes0*mmes2*mwidth2/pi
         weight=weight*mmes1**2*mmes2**2/mmes0**4
         weight=weight/((mmes0**2-mmes1**2)**2+mmes1**4*
     &        mwidth1**2/mmes0**2)
         weight=weight/((mmes0**2-mmes2**2)**2+mmes2**4*
     &        mwidth2**2/mmes0**2)
         
      endif


cccccccccccccccccccccccccccccccccccccccccccccccccccc

c     pi rapidities

         call r2455(ranx1)
         call r2455(ranx2)
         call r2455(ranx3)
         call r2455(ranx4)

         phix1=2d0*pi*ranx1

c         r2max=1d0
c         r2=r2max*ranx1
c         ptxsq1=-dlog(r2)/bjac1  ! generate exp p_t^2 (can be more efficient)

         ptxsqmin=0d0     
         ptxsqmax=10d0         ! generate linear p_t^2

         ptxsq1=ptxsqmin+(ptxsqmax-ptxsqmin)*ranx2
 
         q(1,5)=dsqrt(ptxsq1)*dcos(phix1)
         q(2,5)=dsqrt(ptxsq1)*dsin(phix1)
         q(1,6)=-pt1(1)-pt2(1)-q(1,5)
         q(2,6)=-pt1(2)-pt2(2)-q(2,5)
         ptxsq2=q(1,6)**2+q(2,6)**2
         
         rmx1=dsqrt(mmes1**2+ptxsq1) 
         rmx2=dsqrt(mmes2**2+ptxsq2) 

         ymax1=dlog(rts/rmx1)

         if(cuts.eq.'true')then
            if(rmax.lt.ymax1)then
               ymax1=rmax
            endif
         endif

         ymin1=-ymax1

         yx1=ymin1+(ymax1-ymin1)*ranx3
    
         ymax2=(dlog((rts-rmx1*dexp(yx1))/rmx2))
         ymin2=-(dlog((rts-rmx1*dexp(-yx1))/rmx2))

         if(cuts.eq.'true')then
            if(ymax2.gt.rmax)then
               ymax2=rmax
            endif
            if(ymin2.lt.-rmax)then
               ymin2=-rmax
            endif
         endif

         if(ymax2.lt.ymin2)then
            weight=0d0
            goto 700
         endif

         yx2=ymin2+(ymax2-ymin2)*ranx4
         x1=(rmx2*dexp(yx2)+rmx1*dexp(yx1))/rts
         x2=(rmx2*dexp(-yx2)+rmx1*dexp(-yx1))/rts

         weight=weight*(ptxsqmax-ptxsqmin)
c         weight=weight*dexp(bjac1*ptxsq1)*r2max/bjac1
         weight=weight*(ymax1-ymin1)*(ymax2-ymin2)


         q(3,5)=rmx1*(dexp(yx1)-dexp(-yx1))/2d0
         q(4,5)=rmx1*(dexp(yx1)+dexp(-yx1))/2d0

         q(3,6)=rmx2*(dexp(yx2)-dexp(-yx2))/2d0
         q(4,6)=rmx2*(dexp(yx2)+dexp(-yx2))/2d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     impose massive on-shell condition by solving
c                   p1+ + cc1/p2- = aa1
c                   p2- + cc2/p1+ = aa2 

      aa1=bp*(1d0-x1)
      aa2=bp*(1d0-x2)
      cc1=0.5d0*(pt2sq+mp**2)
      cc2=0.5d0*(pt1sq+mp**2)

      root1sq=(cc1-cc2-aa1*aa2)**2-4d0*cc2*aa1*aa2
      root2sq=(cc2-cc1-aa1*aa2)**2-4d0*cc1*aa1*aa2
      if(root1sq.le.0d0.or.root2sq.le.0d0)then
         weight=0d0
         goto 700
      endif
      p1p=(cc2-cc1+aa1*aa2+dsqrt(root1sq))/(2d0*aa2)
      p2m=(cc1-cc2+aa1*aa2+dsqrt(root2sq))/(2d0*aa1)
      p1m=(pt1sq+mp**2)/(2d0*p1p)
      p2p=(pt2sq+mp**2)/(2d0*p2m)

      if(p1p.lt.0d0.or.p1m.lt.0d0.or.p2p.lt.0d0.or.p2m.lt.0d0)then
         weight=0d0
         goto 700
      endif

      t1=-rts*p1m*rt2
      t2=-rts*p2p*rt2

      q(1,3)=pt1(1)
      q(2,3)=pt1(2)
      q(3,3)=(p1p-p1m)/rt2
      q(4,3)=(p1p+p1m)/rt2

      q(1,4)=pt2(1)
      q(2,4)=pt2(2)
      q(3,4)=(p2p-p2m)/rt2
      q(4,4)=(p2p+p2m)/rt2

ccccccccccccccccccccccccccccccccccccccccccccccccc

      do k=1,4
         svec(k)=q(k,5)+q(k,6)
      enddo

         call dot(svec,svec,sh)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     rho0 decays

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
         if(pflag.eq.'rho')then
            
            kk=6

            do k=5,6
               
               if(k.eq.5)then
                  mmesp=mmes1
               else
                  mmesp=mmes2
               endif

               call rambo(npart,mmesp,am,pout,wt) ! call RAMBO  
              
               do j=1,4
                  pboo(j)=q(j,k)
               enddo
               
               do h=1,2
                  do j=1,4
                     pcm(j)=pout(j,h)
                  enddo
                  call boost(mmesp,pboo,pcm,plb)
                  kk=kk+1
                  do j=1,4
                     q(j,kk)=plb(j)
                  enddo
               enddo
               
            enddo

         endif


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      place cuts

         if(cuts.eq.'true')then

            call cut(icut)
            
            if(icut.eq.0)then
               weight=0d0
               weight0=0d0
               weight1=0d0
               weight2=0d0
               ncut=ncut+1
               goto 700
            endif

         endif
      
c         print*,i

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Write 4-momenta to event record array

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if(output.eq.'lhe')then

         do k=3,4
            do j=1,4
               pup(j,k)=q(j,k)
            enddo
            pup(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
         enddo
         
         if(pflag.eq.'rho')then            

            do k=5,10
               do j=1,4
                  pup(j,k)=q(j,k)
               enddo
               pup(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
            enddo

         else
            
            do k=5,6
               do j=1,4
                  pup(j,k)=q(j,k)
               enddo
               pup(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
            enddo
            
         endif

         else

         do k=3,4
            do j=1,4
               phep(j,k)=q(j,k)
            enddo
            phep(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
         enddo
         
         if(pflag.eq.'rho')then            

            do k=5,10
               do j=1,4
                  phep(j,k)=q(j,k)
               enddo
               phep(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
            enddo

         else
            
            do k=5,6
               do j=1,4
                  phep(j,k)=q(j,k)
               enddo
               phep(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
            enddo
            
         endif

         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Weight evaluation

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c     pp-->pp(M_1 M_2) matrix element

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      call schimc(pt1x,pt1y,pt2x,pt2y,zmat)

      weight=weight*cdabs(zmat)**2
      weight=weight/norm**2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      weight=weight/(s**2*(16d0)**3*pi**5)*389d3

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccc   Probability of no additional particles produced in production subprocess

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(fsi.eq.'true')then
         
         if(pflag.eq.'pipm'.or.pflag.eq.'pi0')then
            
            if(dsqrt(sh).gt.mf127)then
               m0=mf127 
               c=0.7d0
               cont=1d0-(m0-2d0*mmes)**2/m0**2
               exn=c*dlog(((dsqrt(sh)-2d0*mmes)**2)/m0**2+cont)
               snp=dexp(-exn)
            else
               snp=1d0
            endif
            
         elseif(pflag.eq.'kpkm'.or.pflag.eq.'ks')then
            
            if(dsqrt(sh).gt.mf1525)then
               m0=mf1525
               c=0.7d0
               cont=1d0-(m0-2d0*mmes)**2/m0**2
               exn=c*dlog(((dsqrt(sh)-2d0*mmes)**2)/m0**2+cont)
               snp=dexp(-exn)
            else
               snp=1d0
            endif

         else

            snp=1d0
            
         endif
         
      else
         
         snp=1d0
         
      endif
      
      if(snp.gt.1d0)then
         snp=1d0
      endif
      
      weight=weight*snp
 
ccccccccccccccc

      if(pflag.eq.'pi0'.or.pflag.eq.'ks'.or.pflag.eq.'rho')then
         weight=weight/2d0    ! symmetry factor
      endif
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

 666  weight=weight

      if(weightm.lt.weight)then
         weightm=weight
      endif


      if(ll.eq.1)then
         
         xwgtup=weight

         wthist=weight/dfloat(ntotal)
         
         call binit(wthist)
         
      else
         
         call r2455(ranhist)
         
         if(ranhist.lt.weight/weightm)then
            
            xwgtup=1d0

            num=num+1

            
c            call binit(sumt/nev)

            write(35,302)num,nup
            
            if(output.eq.'lhe')then

            do j=1,nup
               write(35,300)j,istup(j),idup(j),mothup(1,j),mothup(2,j),
     &              icolup(1,j),icolup(2,j),pup(1,j),pup(2,j),
     &              pup(3,j),pup(4,j),pup(5,j),vtimup(j),spinup(j)
            enddo
            
            else

               do j=1,nup
                  write(35,301)j,idhep(j),isthep(j),jmohep(1,j),
     &                 jmohep(2,j),jdahep(1,j),jdahep(2,j),
     &                 phep(1,j),phep(2,j),phep(3,j),phep(4,j),phep(5,j)
     &                 ,vhep(1,j),vhep(2,j),vhep(3,j),vhep(4,j)
               enddo

            endif


         endif

      endif

      if(ll.eq.2)then

      if(num.gt.nev-1)then  ! exit loop once required no. of unweighted events generated

         ntotal=i
         goto 888

      endif
      endif

 700  sum=sum+weight
      sum1=sum1+weight**2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c                                                                      
c     End of event loop                                              c
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      enddo
   
 888  sum=sum/dfloat(ntotal)
      sum1=sum1/dfloat(ntotal)
      var=dsqrt((sum1-sum**2)/dfloat(ntotal))
      EFF=1d0*(ntotal-ncut)/ntotal  

      if(ll.eq.1)then
         print*,'...done!'
         write(6,*)
         sumt=sum
      else
         print*,'...done!'
         write(6,*)
      endif


 300  format(i4,1x,i4,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,E13.6,1x,
     &E13.6,1x,E13.6,1x,E13.6,1x,E13.6,1x,E13.6,1x,E13.6)
 301  format(i4,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4,1x,E13.6,1x,E13.6
     &,1x,E13.6,1x,E13.6,1x,E13.6,1x,E13.6,1x,E13.6,1x,E13.6,1x
     &,E13.6)
 302  format(1x,i8,1x,i5,1x,E13.6)

ccccc     create histograms

      if(ll.eq.1)then

      do i=1,nhist
           call histo2(i,0)
      enddo

      endif
      
      if(ll.eq.1)then
      if(pflag.eq.'pipm')then
         write(6,*) 
         write(6,*)'pi+pi- production'
      elseif(pflag.eq.'pi0')then
         write(6,*) 
         write(6,*)'pi0pi0 production'
      elseif(pflag.eq.'kpkm')then
         write(6,*)
         write(6,*)'K+K- production'
      elseif(pflag.eq.'ks')then
         write(6,*)
         write(6,*)'K0K0 production'
      elseif(pflag.eq.'rho')then
         write(6,*)
         write(6,*)'rho0rho0 production'
      endif
      endif

      if(ll.eq.1)then
         if(pflag.eq.'rho')then
         write(6,221)rts,ntotal,sum,var
      else
         write(6,222)rts,mmes,ntotal,sum,var
      endif
      else
         if(output.eq.'lhe')then
            print*,'LHE output'
         else
            print*,'HEPEVT output'
         endif
         write(6,223)nev,sum,var
      endif

      enddo

      close(35)

 221  format(/,
     .       3x,'collider energy  (GeV)    : ',f10.3,/,
     .       3x,'number of events          : ',i9,/,
     .       3x,'sigma (nb)                : ',f15.5,'  +-  ',f15.5,
     .   /)

 222  format(/,
     .       3x,'collider energy  (GeV)    : ',f10.3,/,
     .       3x,'meson mass (GeV)          : ',f10.5,/,
     .       3x,'number of events          : ',i9,/,
     .       3x,'sigma (nb)                : ',f15.5,'  +-  ',f15.5,
     .   /)

 223  format(3x,'number of events          : ',i9,/,
     .       3x,'sigma (nb)                : ',f15.5,'  +-  ',f15.5,
     .   /)

      call cpu_time(t2)
      print*,'time elapsed = ', t2, ' s'

      stop
      end

ccccc Pomeron -- (off-shell) meson form factor

      function fpi(x)             
      implicit double precision(a-y)
      character formf*10
      common/vars/s,rts,mmes,yx
      common/ff/formf
      common/ffpars/bexp,ao,bo,aoo,boo

      if(formf.eq.'exp')then
         fpi=exp((x-mmes**2)*bexp)
      elseif(formf.eq.'orexp')then
         fpi=exp(-(dsqrt(-x+mmes**2+ao**2)-ao)*bo)
      elseif(formf.eq.'power')then
         fpi=1d0/(1d0-(x-mmes**2)/aoo)
      endif

      return
      end

ccccccc  Pom Pom --> meson pair amplitude

      subroutine wev(matt,v1,v2)
      implicit double precision(a-y)
      implicit complex*16(z)
      double precision q(4,20),svec(4),tvec(4),uvec(4),v1(4),v2(4)
      integer k
      character*10 mregge,pflag
      complex*16 matt
      common/mom/q
      common/vars/s,rts,mmes,yx
      common/wvars/sig0,bb,t1,t2
      common/alphas/alphap,alpha0,alphapr,alpha0r,alphapm,alpha0m
      common/sec/cpom,cf,crho,aff,ar
      common/hvars/sh,th,uh
      common/regge/mregge
      common/flags/pflag

      zi=(0d0,1d0)
      pi=dacos(-1d0)

      do k=1,4
         q(k,10)=v1(k)
         q(k,11)=v2(k)
      enddo

      do k=1,4
         tvec(k)=q(k,10)-q(k,5)
         uvec(k)=q(k,10)-q(k,6)
      enddo

         call dot(tvec,tvec,th)
         call dot(uvec,uvec,uh)

         sht1=2d0*(q(4,5)*q(4,3)-q(3,5)*q(3,3)-q(2,5)*q(2,3)-
     &q(1,5)*q(1,3))+mmes**2

         sht2=2d0*(q(4,6)*q(4,4)-q(3,6)*q(3,4)-q(2,6)*q(2,4)-
     &q(1,6)*q(1,4))+mmes**2

         shu1=2d0*(q(4,6)*q(4,3)-q(3,6)*q(3,3)-q(2,6)*q(2,3)-
     &q(1,6)*q(1,3))+mmes**2

         shu2=2d0*(q(4,5)*q(4,4)-q(3,5)*q(3,4)-q(2,5)*q(2,4)-
     &q(1,5)*q(1,4))+mmes**2

       bt1=(2d0*alphap*dlog(sht1))
       bt2=(2d0*alphap*dlog(sht2))
       bu1=(2d0*alphap*dlog(shu1))
       bu2=(2d0*alphap*dlog(shu2))

       bt1r=(2d0*alphapr*dlog(sht1))
       bt2r=(2d0*alphapr*dlog(sht2))
       bu1r=(2d0*alphapr*dlog(shu1))
       bu2r=(2d0*alphapr*dlog(shu2))

       sig0t1=sig0*sht1**0.0808d0/0.389d0
       sig0t2=sig0*sht2**0.0808d0/0.389d0
       sig0u1=sig0*shu1**0.0808d0/0.389d0
       sig0u2=sig0*shu2**0.0808d0/0.389d0

       t11=v1(1)**2+v1(2)**2
       t22=v2(1)**2+v2(2)**2

          zmt=(zi*dexp(-bt1*t11/2d0)*cpom*sht1**alpha0
     &         +((aff+zi)*cf*sht1**alpha0r+(ar-zi)*crho*sht1**alpha0r)
     &         *dexp(-bt1r*t11/2d0))
          zmt=zmt*(zi*dexp(-bt2*t22/2d0)*cpom*sht2**alpha0
     &         +((aff+zi)*cf*sht2**alpha0r-(ar-zi)*crho*sht2**alpha0r)
     &         *dexp(-bt2r*t11/2d0))
          zmt=zmt*fpi(th)**2/(mmes**2-th)
          
          zmu=(zi*dexp(-bu1*t11/2d0)*cpom*shu1**alpha0
     &         +((aff+zi)*cf*shu1**alpha0r-(ar-zi)*crho*shu1**alpha0r)
     &         *dexp(-bu1r*t11/2d0))
          zmu=zmu*(zi*dexp(-bu2*t22/2d0)*cpom*shu2**alpha0
     &         +((aff+zi)*cf*shu2**alpha0r+(ar-zi)*crho*shu2**alpha0r)
     &         *dexp(-bu2r*t22/2d0))
          zmu=zmu*fpi(uh)**2/(mmes**2-uh)
          
          matt=zmu+zmt
      
       return
       end

c     binning subroutine

      subroutine binit(wt)
      implicit double precision(a-y)
      double precision q(4,20),pt1(2),pt2(2),ptx(2)
      common/mom/q
      common/mompt/pt1,pt2,ptx     
      common/vars/s,rts,mmes,yx
      common/vars1/ptgam,etagam,ptel2,ptel1,etael1,etael2
      common/vars2/ptpi1,etapi1,ptpi2,etapi2 
      common/ang/cost,costa,costb
      common/hvars/sh,th,uh

      ypisum=0.5d0*dlog((q(4,5)+q(4,6)+q(3,5)+q(3,6))
     &     /(q(4,5)+q(4,6)-q(3,5)-q(3,6)))

      ypi1=0.5d0*dlog((q(4,5)+q(3,5))
     &     /(q(4,5)-q(3,5)))
      ypi2=0.5d0*dlog((q(4,6)+q(3,6))
     &     /(q(4,6)-q(3,6)))

      call histo1(1,30,0d0,4.5d0,dsqrt(sh),wt)

      return
      end

ccccc Cutting subroutine
ccccc
ccccc User can define their own cuts on any particle momenta, stored in array q(i,j)
ccccc

      subroutine cut(icut)
      implicit double precision(a-y)
      double precision q(4,20)
      integer icut
      common/mom/q
      common/vars/s,rts,mmes,yx
      common/vars1/ptgam,etagam,ptel2,ptel1,etael1,etael2
      common/vars2/ptpi1,etapi1,ptpi2,etapi2      
      common/cuts/etaelmax,etaelmin,ptelmin,ptphmin,ecut,rmax,rmin,mcut
      common/hvars/sh,th,uh
      integer iii
      common/ivar/iii

      icut=0

c -- meson 1 variables
      ptpi1=dsqrt(q(1,5)**2+q(2,5)**2)
      pmodpi1=dsqrt(q(1,5)**2+q(2,5)**2+q(3,5)**2)
      etapi1=.5d0*dlog((pmodpi1+q(3,5))
     .              /(pmodpi1-q(3,5)))

      ypi1=.5d0*dlog((q(4,5)+q(3,5))
     .              /(q(4,5)-q(3,5)))

        etpi1=q(4,5)*dsqrt(pmodpi1**2-q(3,5)**2)/pmodpi1

c -- meson 2 variables
      ptpi2=dsqrt(q(1,6)**2+q(2,6)**2)
      pmodpi2=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
      etapi2=.5d0*dlog((pmodpi2+q(3,6))
     .              /(pmodpi2-q(3,6)))

      ypi2=.5d0*dlog((q(4,6)+q(3,6))
     .              /(q(4,6)-q(3,6)))

        etpi2=q(4,6)*dsqrt(pmodpi2**2-q(3,6)**2)/pmodpi2

c        print*,ptpi1,ptpi2
        
ccccc    example cuts....

c        if(dsqrt(sh).lt.0.8d0)return

        xf1=2d0*q(3,3)/rts
        xf2=2d0*q(3,4)/rts

c        if(dabs(xf1).lt.0.9d0)return
c        if(dabs(xf2).lt.0.9d0)return
        
ccccccccc  Basic cuts described at start of code

        if(ptpi1.lt.ecut)return
        if(ptpi2.lt.ecut)return
        if(ypi1.gt.rmax)return
        if(ypi2.gt.rmax)return
        if(ypi1.lt.rmin)return
        if(ypi2.lt.rmin)return

        icut=1

      return
      end
 
c     prints histograms

      subroutine histo1(ih,ib,x0,x1,x,w)
      implicit real*8(a-h,o-z)
      character*1 regel(30),blank,star
      dimension h(20,100),hx(20),io(20),iu(20),ii(20)
      dimension y0(20),y1(20),ic(20)
      data regel / 30*' ' /,blank /' ' /,star /'*'/
      save
      open(10,file="output.dat")
      y0(ih)=x0
      y1(ih)=x1
      ic(ih)=ib
      if(x.lt.x0) goto 11
      if(x.gt.x1) goto 12
      ix=idint((x-x0)/(x1-x0)*dble(ib))+1
      h(ih,ix)=h(ih,ix)+w
      if(h(ih,ix).gt.hx(ih)) hx(ih)=h(ih,ix)
      ii(ih)=ii(ih)+1
      return
   11 iu(ih)=iu(ih)+1
      return
   12 io(ih)=io(ih)+1
      return
      entry histo2(ih,il)
      ib1=ic(ih)
      x01=y0(ih)
      x11=y1(ih)
      bsize=(x11-x01)/dble(ib1)
      hx(ih)=hx(ih)*(1.d0+1.d-06)
      if(il.eq.0) write(6,21) ih,ii(ih),iu(ih),io(ih)
      if(il.eq.1) write(6,22) ih,ii(ih),iu(ih),io(ih)
   21 format(' no.',i3,' lin : inside,under,over ',3i6)
   22 format(' no.',i3,' log : inside,under,over ',3i6)
      if(ii(ih).eq.0) goto 28
      write(6,23)
   23 format(35(1h ),3(10h----+----i))
      do 27 iv=1,ib1
      z=(dble(iv)-0.5d0)/dble(ib1)*(x11-x01)+x01
      if(il.eq.1) goto 24
      iz=idint(h(ih,iv)/hx(ih)*30.)+1
      goto 25
   24 iz=-1
      if(h(ih,iv).gt.0.d0)
     .iz=idint(dlog(h(ih,iv))/dlog(hx(ih))*30.)+1
   25 if(iz.gt.0.and.iz.le.30) regel(iz)=star
      write(6,26) z,h(ih,iv)/bsize,(regel(i),i=1,30)
c      write(10,*)z-0.125d0/2d0,h(ih,iv)/bsize        ! Print histogram to file
      write(10,*)z,h(ih,iv)/bsize        ! Print histogram to file
c      write(10,*)z,h(ih,iv)/bsize        ! Print histogram to file
   26 format(1h ,2g15.6,4h   i,30a1,1hi)
   36 format(1h ,2g15.6)
      if(iz.gt.0.and.iz.le.30) regel(iz)=blank
   27 continue
      write(6,23)
      return
   28 write(6,29)
   29 format('  no entries inside histogram')
      return
      entry histo3(ih)
      do 31 i=1,100
   31 h(ih,i)=0.
      hx(ih)=0.
      io(ih)=0
      iu(ih)=0
      ii(ih)=0
      close(10)
      return 
      end


cccc  Initializes soft model parameters

      subroutine initpars(in)
      implicit double precision(a-y)
      integer in,i1,i2
      integer nch
      double precision cc0(5),bm(5),bb0(5),pp0(5),bex(5),gaa(5)
      common/pars/cc0,bm,bb0,pp0,bex,asp,sigo,gaa,ep,norm
      common/ipars/nch
      common/vars/s,rts,mmes,yx

      pi=dacos(-1d0)

      if(in.eq.1)then

         ep=0.13d0
         asp=0.08d0
         ep1=0d0
         sigo=23d0/0.39d0
         gaa(3)=0.4d0
         gd2=0.3d0
         nch=2
         ntf=0
         
         cc0(1)=0.45d0
         bm(1)=3d0
         bb0(1)=0.1d0
         pp0(1)=0.92d0
         bex(1)=8.5d0
         
         cc0(2)=0.45d0
         bm(2)=1.5d0
         bb0(2)=0.5d0
         pp0(2)=0.1d0
         bex(2)=4.5d0
         
         cc0(3)=1d0
         bm(3)=0d0
         bb0(3)=0.8d0
         pp0(3)=0.5d0
         bex(3)=0.5d0

      elseif(in.eq.2)then
         
         ep=0.115d0
         asp=0.11d0
         ep1=0d0
         sigo=33d0/0.39d0
         gaa(3)=0.6d0
         gd2=0.16d0
         nch=2d0
         ntf=0d0

         cc0(1)=0.63d0
         bm(1)=3d0
         bb0(1)=0.1d0
         pp0(1)=0.5d0
         bex(1)=8d0
         
         cc0(2)=0.47d0
         bm(2)=1.5d0
         bb0(2)=0.5d0
         pp0(2)=0.1d0
         bex(2)=6d0
         
         cc0(3)=1d0
         bm(3)=0d0
         bb0(3)=0.8d0
         pp0(3)=0.5d0
         bex(3)=0.5d0
        
      elseif(in.eq.3)then

         ep=0.093d0
         asp=0.075d0
         ep1=0d0
         sigo=60d0/0.39d0
         gaa3=4.8d0
         gd2=1.03d0
         nch=2
         nga=1
         cc0(1)=0.55d0
         bm(1)=3d0
         bb0(1)=0.27d0
         pp0(1)=0.48d0
         bex(1)=5.3d0
         
         cc0(2)=0.48d0
         bm(2)=1.5d0
         bb0(2)=0.1d0
         pp0(2)=1d0
         bex(2)=3.8d0
         
         cc0(3)=0.24d0
         
      elseif(in.eq.4)then

         ep=0.11d0              !!! Capital delta
         asp=0.06d0             !!! alpha'
         ep1=0d0                !!! Zero in all models, matters (?)
         sigo=50d0/0.39d0       !!! sigma_0 (GeV^-2)
         gaa3=6d0               !!! k2/k(1.8 TeV)
         gd2=1.3d0              !!! k1/k(1.8 TeV)
         nch=2
         nga=1                  !!! A flag (?)
         cc0(1)=0.6d0           !!! d1
         bm(1)=3d0              !!! Doesn't matter
         bb0(1)=0.45d0          !!! c1-0.08 (added back later)
         pp0(1)=0.5d0           !!! 2*|a_1|^2
         bex(1)=7.2d0           !!! b1
         
         cc0(2)=0.48d0          !!! d2
         bm(2)=1.5d0            !!! Doesn't matter
         bb0(2)=0.16d0          !!! c2-0.08 (added back later)
         pp0(2)=1d0             !!! |a_2|^2 is set later
         bex(2)=4.2d0           !!! b2
         
         cc0(3)=0.12d0          !!! Beta : k^2_min ~ s^Beta
      
      endif

      if(nch.eq.3) pp0(3)=3d0-pp0(2)-pp0(1)
      if(nch.eq.2) pp0(2)=2d0-pp0(1) !!! Set |a_2|^2
      if(nch.eq.1) pp0(1)=1d0

      if(in.eq.3.or.in.eq.4)then

      gamm=(1800d0/rts)**cc0(3)
      ga1=1d0/(1d0+gamm*gd2)
      ga2=1d0/(1d0+gamm*gaa3)
      gaa(1)=2d0*ga1/(ga1+ga2)   !!! gamma_1 (+) Eq (15)
      gaa(2)=2d0*ga2/(ga1+ga2)   !!! gamma_2 (-) Eq (15)

      elseif(in.eq.1.or.in.eq.2)then

         if(nch.eq.2)then
            gaa(1)=1d0+dsqrt(gd2)
            gaa(2)=1d0-dsqrt(gd2)
            gaa(3)=0.
         elseif(nch.eq.1)then
            gaa(1)=1d0
            gaa(2)=0d0
            gaa(3)=0d0
         endif

      endif

      sum=0d0

      do i1=1,nch
         do i2=1,nch 
            sum=sum+gaa(i1)*gaa(i2)*pp0(i1)*pp0(i2)/dble(nch)**2
         enddo
      enddo
      
      norm=sum

      return
      end

ccccc Calculates proton opacity

      subroutine calcop
      implicit double precision(a-y)
      integer nb,i1,i2,nch,ib
      common/ipars/nch
      double precision op(5,5,10000,2),oph(5,5,10000,2)
      common/opac/op,oph

      nb=900
      hb=100d0/dble(nb)

      print*,'Calculating opacity'


      do ib=1,nb+1
         bt=dble(ib-1)*hb

         do i1=1,nch
            do i2=1,nch

               call opacity(i1,i2,bt,fr,fr1)


               op(i1,i2,ib,1)=bt
               op(i1,i2,ib,2)=fr
               oph(i1,i2,ib,1)=bt
               oph(i1,i2,ib,2)=fr1

               write(40,*)bt,fr,fr1

            enddo
         enddo
      enddo

      return
      end

cccc  Calculates screening amplitude (in k_t space)

      subroutine calcscreen
      implicit double precision(a-y)
      integer ns,i1,i2,nch,ib
      common/ipars/nch
      double precision sca(5,5,0:40000,2),sca1(5,5,0:40000,2)
      common/spac/sca,sca1

      ns=900
      ksqma=8.2d0
      inck=ksqma/dble(ns)
      ksqmin=0.001d0
      lginck=(dlog(ksqma/ksqmin))/dble(ns)

      print*,'Calculating screening amplitude'

      do ib=0,ns+1

         ksq=dble(ib-1)*inck
         lgksq=ksq

         if(ib.eq.0)then
            ksq=0d0
            lgksq=0d0
         else
            lgksq=dble(ib-1)*lginck+dlog(ksqmin)
            ksq=dexp(lgksq)
         endif

         do i1=1,nch
            do i2=1,nch

               call screening(i1,i2,ksq,sc,sc1)
      
               sca(i1,i2,ib,1)=lgksq
               sca(i1,i2,ib,2)=sc
               sca1(i1,i2,ib,1)=lgksq
               sca1(i1,i2,ib,2)=sc1

            enddo
         enddo

      enddo

      return
      end

ccccc Interpolation....

      subroutine opacityint(i,j,bt,out,out1)
      implicit double precision(a-y)
      double precision op(5,5,10000,2),oph(5,5,10000,2)
      integer i,j,it
      common/opac/op,oph

      incbt=op(1,1,2,1)-op(1,1,1,1)
      it=nint(bt/incbt)
      if(dble(it).gt.(bt/incbt))then
         it=it-1
      endif

      m=(op(i,j,it+2,2)-op(i,j,it+1,2))
     &/(op(i,j,it+2,1)-op(i,j,it+1,1))
      del=bt-op(1,1,it+1,1)
      mh=(oph(i,j,it+2,2)-oph(i,j,it+1,2))
     &/(oph(i,j,it+2,1)-oph(i,j,it+1,1))
      delh=bt-oph(1,1,it+1,1)

      out=m*del+op(i,j,it+1,2)
      out1=mh*delh+oph(i,j,it+1,2)

      return
      end

      subroutine screeningint(i,j,ktsq,out,out1)
      implicit double precision(a-y)
      double precision sca(5,5,0:40000,2),sca1(5,5,0:40000,2)
      integer i,j,it,ns
      common/spac/sca,sca1

      if(ktsq.lt.0.001d0)then

         it=0

         m=(sca(i,j,it+1,2)-sca(i,j,it,2))
     &        /(dexp(sca(i,j,it+1,1))-sca(i,j,it,1))
         del=ktsq-sca(1,1,it,1)
         m1=sca1(i,j,it+1,2)-sca1(i,j,it,2)
     &        /(dexp(sca1(i,j,it+1,1))-sca1(i,j,it,1))
         del1=ktsq-sca1(1,1,it,1)

         out=m*del+sca(i,j,it,2)
         out1=m1*del1+sca1(i,j,it,2)

      elseif(ktsq.lt.8d0)then

         ktmin=sca(1,1,1,1)
         inckt=sca(1,1,2,1)-sca(1,1,1,1)
         it=nint((dlog(ktsq)-ktmin)/inckt)
         if(dble(it).gt.((dlog(ktsq)-ktmin)/inckt))then
            it=it-1
         endif

         m=(sca(i,j,it+2,2)-sca(i,j,it+1,2))
     &        /(sca(i,j,it+2,1)-sca(i,j,it+1,1))
         del=dlog(ktsq)-sca(1,1,it+1,1)
         m1=sca1(i,j,it+2,2)-sca1(i,j,it+1,2)
     &        /(sca1(i,j,it+2,1)-sca1(i,j,it+1,1))
         del1=dlog(ktsq)-sca1(1,1,it+1,1)

         out=m*del+sca(i,j,it+1,2)
         out1=m1*del1+sca1(i,j,it+1,2)

      else

         out=0d0
         out1=0d0

      endif

      return
      end

cccc  Integrates round Pomeron loop (to calculate screened amplitude)

      subroutine schimc(p1x,p1y,p2x,p2y,out)
      implicit double precision(a-y)
      integer nx,jx,jy,nch,i1,i2,it,k
      double precision cc0(5),bm(5),bb0(5),pp0(5),bex(5),gaa(5)
      double precision sca(5,5,0:40000,2),sca1(5,5,0:40000,2)
      double precision a1(4),a2(4),q(4,20)
      complex*16 out,x0,x00
      common/pars/cc0,bm,bb0,pp0,bex,asp,sigo,gaa,ep,norm
      common/ipars/nch
      common/spac/sca,sca1
      common/mom/q
      common/wvars/sig0,bb,t1,t2
      common/alphas/alphap,alpha0,alphapr,alpha0r,alphapm,alpha0m

      do k=1,4
         a1(k)=q(k,1)-q(k,3)
         a2(k)=q(k,2)-q(k,4)
      enddo

      nx=10
      hx=2d0/dble(nx)

      out=0d0

      do jx=-nx,nx
         do jy=-nx,nx

            tpx=dble(jx)*hx
            tpy=dble(jy)*hx
            tp2=tpx**2+tpy**2
            wt=hx*hx

            p1xp=p1x-tpx
            p1yp=p1y-tpy
            t12=p1xp**2+p1yp**2
            p2xp=tpx+p2x
            p2yp=tpy+p2y
            t22=p2xp**2+p2yp**2

            a1(1)=-p1xp
            a1(2)=-p1yp
            a2(1)=-p2xp
            a2(2)=-p2yp

      do i1=1,nch
         do i2=1,nch

           call screeningint(i1,i2,tp2,sc,sc1)     
           call wev(x0,a1,a2)

            x0=x0*pp0(i1)*pp0(i2)/dble(nch)**2
            x0=x0*dexp(-((t12+0.08d0+bb0(i1))*bex(i1))**cc0(i1)+
     &     (bex(i1)*(bb0(i1)+0.08d0))**cc0(i1))
            x0=x0*dexp(-((t22+0.08d0+bb0(i2))*bex(i2))**cc0(i2)+
     &     (bex(i2)*(bb0(i2)+0.08d0))**cc0(i2)) 

            out=out+x0*wt*sc
            
         enddo
      enddo


      enddo
      enddo

      a1(1)=-p1x
      a1(2)=-p1y
      a2(1)=-p2x
      a2(2)=-p2y

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2
      
      call formfac(t11,t22,x00p)
      
      call wev(x00,a1,a2)

      

      out=out+x00*x00p

      return
      end

ccccc Pomeron -- diffrative eignstate i form factor

      subroutine formfac(t1,t2,out)
      implicit double precision(a-y)
      double precision cc0(5),bm(5),bb0(5),pp0(5),bex(5),gaa(5)
      common/pars/cc0,bm,bb0,pp0,bex,asp,sigo,gaa,ep,norm
      integer nch,i1,i2
      common/ipars/nch

      out=0d0

      delta1=dexp(-t1/2d0)
      delta2=dexp(-t2/2d0)

      do i1=1,nch
         do i2=1,nch
            
            wt=gaa(i1)*gaa(i2)*pp0(i1)*pp0(i2)/dble(nch)**2
            wt=wt*dexp(-((t1+0.08d0+bb0(i1))*bex(i1))**cc0(i1)+
     &     (bex(i1)*(bb0(i1)+0.08d0))**cc0(i1))
            wt=wt*dexp(-((t2+0.08d0+bb0(i2))*bex(i2))**cc0(i2)+
     &     (bex(i2)*(bb0(i2)+0.08d0))**cc0(i2))

            out=out+wt

         enddo
      enddo
      
      return
      end

      subroutine screening(i,j,ktsq,out,out1)
      implicit double precision(a-y)
      integer i,j,ib,nb,it,nt,nch
      double precision cc0(5),bm(5),bb0(5),pp0(5),bex(5),gaa(5)
      common/pars/cc0,bm,bb0,pp0,bex,asp,sigo,gaa,ep,norm
      common/ipars/nch
      common/vars/s,rts,mmes,yx

      pi=dacos(-1d0)

      nb=5000
      hb=99d0/dble(nb)  

      out=0d0
      out1=0d0

      do ib=1,nb
         bt=dble(ib)*hb   
         wt=-bt/2d0/pi*hb
 
         call opacityint(i,j,bt,fr,fr1)
      
         sige=sigo*dexp(dlog(rts)*2d0*ep)
         fr=fr*gaa(i)*gaa(j)*sige
         fr1=fr1*gaa(i)*gaa(j)*sige
         
         out=out+wt*(1d0-dexp(-fr/2d0))*besj0(bt*dsqrt(ktsq))
     &        *gaa(i)*gaa(j)
         out1=out1+wt*(1d0-dexp(-fr1/2d0))*besj0(bt*dsqrt(ktsq))
     &        *gaa(i)*gaa(j)

      enddo

      return
      end


      subroutine opacity(i,j,bt,out,out1)
      implicit double precision(a-y)
      integer i,j,ib,nb,it,nt,nch
      double precision cc0(5),bm(5),bb0(5),pp0(5),bex(5),gaa(5)
      common/pars/cc0,bm,bb0,pp0,bex,asp,sigo,gaa,ep,norm
      common/ipars/nch
      common/vars/s,rts,mmes,yx

      pi=dacos(-1d0)

      ampi=0.02d0
      amro=1d0
      a4=4d0*ampi
      alo=dlog(amro/ampi)
      bpol=2.4d0

      nt=6000
      htt=6d0/dble(nt)**2

         out=0d0
         out1=0d0

         do it=0,nt
            t=dble(it)**2*htt
            if(it.eq.0) t=1d-8
            wt=htt*2d0*dble(it)/4d0/pi
            if(it.eq.0) wt=wt/2d0
            bes0=besj0(bt*dsqrt(t))

            ffi=dexp(-((t+0.08d0+bb0(i))*bex(i))**cc0(i)+
     &           (bex(i)*(bb0(i)+0.08d0))**cc0(i))
            ffj=dexp(-((t+0.08d0+bb0(j))*bex(j))**cc0(j)+
     &           (bex(j)*(bb0(j)+0.08d0))**cc0(j))

            asp1=asp
            form1=dlog(ffi*ffj)-2d0*t*asp1*dlog(rts)
cccccccccccccccc
            ah=dsqrt(1d0+a4/t)
            h1pi=2d0*a4+t*(alo-ah*ah*ah*dlog((1d0+ah)/(ah-1d0))) !!! Pion loop insertion
            h1pi=h1pi*sigo/(72d0*pi**3*(1d0+t/bpol)**2)
ccccccccccccccc
            ww=bes0*dexp(form1-2d0*h1pi*dlog(rts))
            aspt=t*asp+h1pi

            out=out+ww*wt
            out1=out1+bes0*wt*ffi*ffj

         enddo
         

      return
      end

      function rgauss(r1,r2)
      implicit double precision (a-y)

      pi=dacos(-1d0)
      rgauss=dsqrt(-2d0*dlog(r1))*dsin(2d0*pi*r2)

      return
      end

c     boosting subroutine

      subroutine boost(p,pboo,pcm,plb)
      real*8  pboo(4),pcm(4),plb(4),p,fact
         plb(4)=(pboo(4)*pcm(4)+pboo(3)*pcm(3)
     &             +pboo(2)*pcm(2)+pboo(1)*pcm(1))/p
         fact=(plb(4)+pcm(4))/(p+pboo(4))
         do 10 j=1,3
  10     plb(j)=pcm(j)+fact*pboo(j)
      return
      end

c     calculates lorentzian dot product for real 4-vectors

      subroutine dot(v1,v2,dt)
      double precision v1(4),v2(4),dt

      dt=-(v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
     &-v1(4)*v2(4))

      return
      end
c

*
* subtractive mitchell-moore generator
* ronald kleiss - october 2, 1987
*
* the algorithm is N(i)=[ N(i-24) - N(i-55) ]mod M,
* implemented in a cirucular array with identifcation
* of NR(i+55) and nr(i), such that effectively:
*        N(1)   <---   N(32) - N(1)
*        N(2)   <---   N(33) - N(2)  ....
*   .... N(24)  <---   N(55) - N(24)
*        N(25)  <---   N(1)  - N(25) ....
*   .... N(54)  <---   N(30) - N(54)
*        N(55)  <---   N(31) - N(55)
*
* in this version  M =2**30  and  RM=1/M=2.D0**(-30.D0)
*
* the array NR has been initialized by putting NR(i)=i
* and subsequently running the algorithm 100,000 times.
*

      subroutine R2455(Ran)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION N(55)
      DATA N/
     . 980629335, 889272121, 422278310,1042669295, 531256381,
     . 335028099,  47160432, 788808135, 660624592, 793263632,
     . 998900570, 470796980, 327436767, 287473989, 119515078,
     . 575143087, 922274831,  21914605, 923291707, 753782759,
     . 254480986, 816423843, 931542684, 993691006, 343157264,
     . 272972469, 733687879, 468941742, 444207473, 896089285,
     . 629371118, 892845902, 163581912, 861580190,  85601059,
     . 899226806, 438711780, 921057966, 794646776, 417139730,
     . 343610085, 737162282,1024718389,  65196680, 954338580,
     . 642649958, 240238978, 722544540, 281483031,1024570269,
     . 602730138, 915220349, 651571385, 405259519, 145115737/
      DATA M/1073741824/
      DATA RM/0.9313225746154785D-09/
      DATA K/55/,L/31/
      IF(K.EQ.55) THEN
         K=1
      ELSE
         K=K+1
      ENDIF
      IF(L.EQ.55) THEN
         L=1
      ELSE
         L=L+1
      ENDIF
      J=N(L)-N(K)
      IF(J.LT.0) J=J+M
      N(K)=J
      RAN=J*RM
      END

      

      SUBROUTINE RAMBO(N,ET,XM,P,WT)
C------------------------------------------------------
C
C                       RAMBO
C
C    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
C
C    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
C    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
C    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
C
C    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
C    ET = TOTAL CENTRE-OF-MASS ENERGY
C    XM = PARTICLE MASSES ( DIM=100 )
C    P  = PARTICLE MOMENTA ( DIM=(4,100) )
C    WT = WEIGHT OF THE EVENT
C
C------------------------------------------------------
      implicit none
!     boost arguments
      integer n
      double precision et,xm(100),p(4,100),wt
!     function
      double precision rn
!     internal variables
      double precision q(4,100),z(100),r(4),
     &     b(3),p2(100),xm2(100),e(100),v(100)
      integer iwarn(5)
      double precision acc
      integer itmax,ibegin
      data acc/1.0d-14/,itmax/6/,ibegin/0/,iwarn/5*0/

      integer i,k,iter
      double precision twopi,po2log,xmt,c,s,f,rmas,g,a,x,bq,
     &     accu,f0,g0,wt2,wt3,wtm,x2,xmax
      integer nm
      save


C
C INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
      IF(IBEGIN.NE.0) GOTO 103
      IBEGIN=1
      TWOPI=8.*DATAN(1.D0)
      PO2LOG=DLOG(TWOPI/4.)
      Z(2)=PO2LOG
      DO 101 K=3,100
  101 Z(K)=Z(K-1)+PO2LOG-2.*DLOG(DFLOAT(K-2))
      DO 102 K=3,100
  102 Z(K)=(Z(K)-DLOG(DFLOAT(K-1)))
C
C CHECK ON THE NUMBER OF PARTICLES
  103 IF(N.GT.1.AND.N.LT.101) GOTO 104
      PRINT 1001,N
      STOP
C
C CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
  104 XMT=0.
      NM=0
      DO 105 I=1,N
      IF(XM(I).NE.0.D0) NM=NM+1
  105 XMT=XMT+DABS(XM(I))
      IF(XMT.LE.ET) GOTO 201
      PRINT 1002,XMT,ET
      STOP
C
C THE PARAMETER VALUES ARE NOW ACCEPTED
C
C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
  201 DO 202 I=1,N
      C=2.*RN(1)-1.
      S=DSQRT(1.-C*C)
      F=TWOPI*RN(2)
      Q(4,I)=-DLOG(RN(3)*RN(4))
      Q(3,I)=Q(4,I)*C
      Q(2,I)=Q(4,I)*S*DCOS(F)
 202  Q(1,I)=Q(4,I)*S*DSIN(F)
 
C
C CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
      DO 203 I=1,4
  203 R(I)=0.
      DO 204 I=1,N
      DO 204 K=1,4
  204 R(K)=R(K)+Q(K,I)
      RMAS=DSQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
      DO 205 K=1,3
  205 B(K)=-R(K)/RMAS
      G=R(4)/RMAS
      A=1./(1.+G)
      X=ET/RMAS
C
C TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
      DO 207 I=1,N
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
      DO 206 K=1,3
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
  207 P(4,I)=X*(G*Q(4,I)+BQ)
C
C CALCULATE WEIGHT AND POSSIBLE WARNINGS
      WT=PO2LOG
      IF(N.NE.2) WT=(2.*N-4.)*DLOG(ET)+Z(N)
      IF(WT.GE.-180.D0) GOTO 208
      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
  208 IF(WT.LE. 174.D0) GOTO 209
      IF(IWARN(2).LE.5) PRINT 1005,WT
      IWARN(2)=IWARN(2)+1
C
C RETURN FOR WEIGHTED MASSLESS MOMENTA
  209 IF(NM.NE.0) GOTO 210
      WT=DEXP(WT)
      RETURN
C
C MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
  210 XMAX=DSQRT(1.-(XMT/ET)**2)
      DO 301 I=1,N
      XM2(I)=XM(I)**2
  301 P2(I)=P(4,I)**2
      ITER=0
      X=XMAX
      ACCU=ET*ACC
  302 F0=-ET
      G0=0.
      X2=X*X
      DO 303 I=1,N
      E(I)=DSQRT(XM2(I)+X2*P2(I))
      F0=F0+E(I)
  303 G0=G0+P2(I)/E(I)
      IF(DABS(F0).LE.ACCU) GOTO 305
      ITER=ITER+1
      IF(ITER.LE.ITMAX) GOTO 304
C      PRINT 1006,ITMAX,ACCU,DABS(F0)
      WRITE(99,1006)ITMAX,ACCU,DABS(F0)
      IF (DABS(F0).GT.1.0E-6) THEN
         WRITE(*,1007)DABS(F0)
      ENDIF
      GOTO 305
  304 X=X-F0/(X*G0)
      GOTO 302
  305 DO 307 I=1,N
      V(I)=X*P(4,I)
      DO 306 K=1,3
  306 P(K,I)=X*P(K,I)
  307 P(4,I)=E(I)
C
C CALCULATE THE MASS-EFFECT WEIGHT FACTOR
      WT2=1.
      WT3=0.
      DO 308 I=1,N
      WT2=WT2*V(I)/E(I)
  308 WT3=WT3+V(I)**2/E(I)
      WTM=(2.*N-3.)*DLOG(X)+DLOG(WT2/WT3*ET)
C
C RETURN FOR  WEIGHTED MASSIVE MOMENTA
      WT=WT+WTM
      IF(WT.GE.-180.D0) GOTO 309
      IF(IWARN(3).LE.5) PRINT 1004,WT
      IWARN(3)=IWARN(3)+1
  309 IF(WT.LE. 174.D0) GOTO 310
      IF(IWARN(4).LE.5) PRINT 1005,WT
      IWARN(4)=IWARN(4)+1
  310 WT=DEXP(WT)
      RETURN
C
 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',
     . ' DESIRED ACCURACY =',D10.3,' (STOPS AT',D10.3,')')
 1007 FORMAT(' RAMBO WARNS: END OF ITERATION ACCURACY TOO LOW =',D10.3)
      END
C

      function rn(idum)
*
* SUBTRACTIVE MITCHELL-MOORE GENERATOR
* RONALD KLEISS - OCTOBER 2, 1987
*
* THE ALGORITHM IS N(I)=[ N(I-24) - N(I-55) ]MOD M,
* IMPLEMENTED IN A CIRUCULAR ARRAY WITH IDENTIFCATION
* OF NR(I+55) AND NR(I), SUCH THAT EFFECTIVELY:
*        N(1)   <---   N(32) - N(1)
*        N(2)   <---   N(33) - N(2)  ....
*   .... N(24)  <---   N(55) - N(24)
*        N(25)  <---   N(1)  - N(25) ....
*   .... N(54)  <---   N(30) - N(54)
*        N(55)  <---   N(31) - N(55)
*
* IN THIS VERSION  M =2**30  AND  RM=1/M=2.D0**(-30D0)
*
* THE ARRAY NR HAS BEEN INITIALIZED BY PUTTING NR(I)=I
* AND SUBSEQUENTLY RUNNING THE ALGORITHM 100,000 TIMES.
*

      implicit none
      double precision rn
      integer idum
      integer n(55)
      data n/
     . 980629335, 889272121, 422278310,1042669295, 531256381,
     . 335028099,  47160432, 788808135, 660624592, 793263632,
     . 998900570, 470796980, 327436767, 287473989, 119515078,
     . 575143087, 922274831,  21914605, 923291707, 753782759,
     . 254480986, 816423843, 931542684, 993691006, 343157264,
     . 272972469, 733687879, 468941742, 444207473, 896089285,
     . 629371118, 892845902, 163581912, 861580190,  85601059,
     . 899226806, 438711780, 921057966, 794646776, 417139730,
     . 343610085, 737162282,1024718389,  65196680, 954338580,
     . 642649958, 240238978, 722544540, 281483031,1024570269,
     . 602730138, 915220349, 651571385, 405259519, 145115737/
      double precision eps
      double precision rm
      integer j,k,l,m

      data eps/1D-9/
      data m/1073741824/
      data rm/0.9313225746154785D-09/
      data k/55/,l/31/
   

  1   CONTINUE
      IF(K.EQ.55) THEN
         K=1
      ELSE
         K=K+1
      ENDIF
      IF(L.EQ.55) THEN
         L=1
      ELSE
         L=L+1
      ENDIF
      J=N(L)-N(K)
      IF(J.LT.0) J=J+M
      N(K)=J
      RN=J*RM
      IF(RN.LT.EPS) GOTO 1
      IF(RN.GT.1D0-EPS) GOTO 1
      RETURN
      END
