c----------------------------------------------------------------------------------
c   urqmdinit ... transfer EPOS -> UrQMD
c   urqmdepos ... main uUrQMD routine
c                    essentially copied from urqmdepos (urqmd.f)
c   urqmdexit ... transfer UrQMD -> EPOS
c   file14outx .. analysis, essentially copied from file14out (output.f)
c----------------------------------------------------------------------------------
c modifs of urqmd routines:
c   upmerge.f: 
c      comment the line "write(6,*) '(Info) pdg2ityp: ..."
c----------------------------------------------------------------------------------

c----------------------------------------------------------------------------------
      subroutine urqmdinit
c----------------------------------------------------------------------------------

      implicit none

      
      
      !urqmd common block

c     debug and validity range

      integer nmax, nspl
      real*8 hit_sphere
      parameter (nmax = 40000) ! maximum number of particles
      parameter (nspl = 500)  ! dimension of spline arrays
      parameter (hit_sphere = 8.d0)  ! hard collision cutoff: 251 mbarn

      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
c 7 integer

      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm



      real*8 
     +     gw, sgw, delr, fdel, dt,
     +     da, db,
     +     Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, gamYuk, drPau, dpPau,
     +     dtimestep
c 19 real*8
     
      real*8 cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      common /cuts/ cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww

      real*8 spx(nspl), spPauy(nspl), outPau(nspl), 
     +                spCby(nspl),  outCb(nspl),
     +                spYuky(nspl), outYuk(nspl),
     +                spSkyy(nspl), outSky(nspl),
     +                spdwwy(nspl), outdww(nspl)

      common /spdata/ spx, spPauy, outPau, spCby,  outCb,
     +                     spYuky, outYuk, spSkyy, outSky,
     +                     spdwwy, outdww

      real*8 
     +     r0(nmax), rx(nmax), ry(nmax), rz(nmax),
     +     p0(nmax), px(nmax), py(nmax), pz(nmax),
     +     airx(nmax), airy(nmax), airz(nmax),
     +     aipx(nmax), aipy(nmax), aipz(nmax),
     +     aorx(nmax,4), aory(nmax,4), aorz(nmax,4),
     +     aopx(nmax,4), aopy(nmax,4), aopz(nmax,4),
     +     fmass(nmax), rww(nmax), 
     +     dectime(nmax), tform(nmax), xtotfac(nmax)
      
      
      integer spin(nmax),ncoll(nmax),charge(nmax),strid(nmax),
     +        ityp(nmax),lstcoll(nmax),iso3(nmax),origin(nmax),uid(nmax)
      common/isys/spin,ncoll,charge,ityp,lstcoll,iso3,origin,strid,
     +            uid
     
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime
      common /frag/ tform, xtotfac


      common /aios/ airx, airy, airz, aipx, aipy, aipz,
     +              aorx, aory, aorz, aopx, aopy, aopz

      common /pots/ Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, 
     +              gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
     +              dt,da, db,dtimestep


c spectator arrays
	integer smax
	parameter(smax=500)  ! maximum number of spectators
	real*8 r0s(smax), rxs(smax), rys(smax), rzs(smax),
     +	       p0s(smax), pxs(smax), pys(smax), pzs(smax),
     +	       sfmass(smax)
	

        integer sspin(smax), scharge(smax), sityp(smax), siso3(smax),
     +          suid(smax)

	integer nspec

	common /scoor/ r0s, rxs, rys, rzs, p0s, pxs ,pys, pzs, sfmass

	common /sisys/ sspin, scharge, sityp, siso3, suid

	common /ssys/ nspec


        real*8 p0td(2,nmax),pxtd(2,nmax),pytd(2,nmax),pztd(2,nmax),
     +         fmasstd(2,nmax)
        integer ityptd(2,nmax),iso3td(2,nmax)
        integer itypt(2),uidt(2),origint(2),iso3t(2)

        common /rtdelay/p0td,pxtd,pytd,pztd,fmasstd
        common /itdelay/ityptd,iso3td
        common /svinfo/itypt,uidt,origint,iso3t
        real*8 ffermpx(nmax), ffermpy(nmax), ffermpz(nmax)
        real*8 peq1, peq2
        common /ffermi/ ffermpx, ffermpy, ffermpz
        common /peq/ peq1,peq2
        
	
	integer numcto,numctp,maxstables
        parameter(numcto=400) ! maximum number of options
        parameter(numctp=400) ! maximum number of parameters
        parameter(maxstables=20) ! maximum number of stable particles
	integer   CTOption(numcto)
	real*8    CTParam(numctp)
        common /options/CTOption,CTParam
	
	real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     +     frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)

      common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz 
	      integer maxbar,maxbra,minbar
      integer offmeson,maxmeson,pimeson,maxbrm,minnuc,mindel
      integer maxbrs1,maxbrs2
      integer numnuc,numdel,nucleon,maxnuc,maxdel
      integer minmes,maxmes


      parameter (minnuc=1) ! lowest baryon particle ID 
      parameter (minmes=100) ! lowest meson particle ID
      parameter (maxmes=132) ! hightest meson particle ID

c number of resonances of a kind
      parameter (numnuc=16) ! number of nucleon resonances
      parameter (numdel=10) ! number of delta resonances
c indices of minimal and maximal itype of a kind (redundant but nice)
      parameter (maxnuc=minnuc+numnuc-1) ! highest nucleon ID
      parameter (mindel=minnuc+maxnuc)   ! lowest delta ID
      parameter (maxdel=mindel+numdel-1) ! highest delta ID

c minres & maxres define the range of nonstable & nonstrange baryons
      integer minres,maxres
      parameter (minres=minnuc+1) ! lowest baryon resonance ID
      parameter (maxres=maxdel)   ! highest (nonstrange) baryon 
                                  ! resonance ID

c strangenes.ne.0 baryon resonances
      integer minlam,minsig,mincas,minome
      integer numlam,numsig,numcas,numome
      integer maxlam,maxsig,maxcas,maxome
      parameter (numlam=13) ! number of lambda states
      parameter (numsig=9)  ! number of sigma states
      parameter (numcas=6)  ! number of cascade states
      parameter (numome=1)  ! number of omega states
      parameter (minlam=mindel+numdel)   ! ID of lowest lambda state
      parameter (maxlam=minlam+numlam-1) ! ID of highest lambda state
      parameter (minsig=minlam+numlam)   ! ID of lowest sigma state
      parameter (maxsig=minsig+numsig-1) ! ID of highest sigma state
      parameter (mincas=minsig+numsig)   ! ID of lowest cascade state
      parameter (maxcas=mincas+numcas-1) ! ID of highest cascade state
      parameter (minome=mincas+numcas)   ! ID of lowest omega state
      parameter (maxome=minome+numome-1) ! ID of highest omega state

c minbar & maxbar define the range of all baryons
      parameter (minbar=minnuc) ! ID of lowest baryon state
      parameter (maxbar=maxome) ! ID of highest baryon state

      parameter (offmeson=minmes) ! offset between zero and lowest 
                                  ! meson state
      parameter (maxmeson=maxmes) ! ID of highest meson state
c... these variables are in principal obsolete and should be exchanged 
c were referenced 

c... avoid hard coded itypes
      integer itrho,itome,iteta,itkaon,itphi,itetapr
      parameter (itkaon=106)   ! ID of kaon
      parameter (itrho=104)    ! ID of rho meson 
      parameter (itome=103)    ! ID of omega meson
      parameter (iteta=102)    ! ID of eta
      parameter (itphi=109)    ! ID of phi
      parameter (itetapr=107)  ! ID of eta'
      parameter (pimeson=101)  ! ID of $\pi$
      parameter (nucleon=minnuc) ! ID of nucleon

      integer itmin,itmax
      parameter (itmin=minnuc)  ! lowest defined ID
      parameter (itmax=maxmes)  ! highest defined ID
c
      parameter (maxbra=11)  ! decay channels for $s=0$ baryon resonances
      parameter (maxbrm=25) ! decay channels for meson resonances
      parameter (maxbrs1=10)! decay channels for $s=1$ baryon resonances
      parameter (maxbrs2=3) ! decay channels for $s=2$ baryon resonances

c 
       integer mlt2it(maxmes-minmes) ! meson IDs sorted by multipletts


      real*8 massoff,mresmin,mresmax
      parameter (massoff=1d-4)      ! offset for mass generation
      parameter (mresmin=1.0765d0)  ! minimum baryon resonance mass
      parameter (mresmax=5d0)       ! maximum baryon resonance mass

      character*45 versiontag
      common /versioning/ versiontag

      real*8 massres(minbar:maxbar),widres(minbar:maxbar)
      real*8 branmes(0:maxbrm,minmes+1:maxmes)
      real*8 branres(0:maxbra,minnuc+1:maxdel)
      real*8 branbs1(0:maxbrs1,minlam+1:maxsig)
      real*8 branbs2(0:maxbrs2,mincas+1:maxcas)
      integer Jres(minbar:maxbar)
      integer Jmes(minmes:maxmes)
      integer pares(minbar:maxbar),pames(minmes:maxmes)
      integer Isores(minbar:maxbar), Isomes(minmes:maxmes)
      integer brtype(4,0:maxbra),bmtype(4,0:maxbrm)
      integer bs1type(4,0:maxbrs1),bs2type(4,0:maxbrs2)
      real*8 massmes(minmes:maxmes)
      real*8 mmesmn(minmes:maxmes)
      real*8 widmes(minmes:maxmes)
      integer strres(minbar:maxbar),strmes(minmes:maxmes)

      integer lbr(0:maxbra,minnuc+1:maxdel)
      integer lbs1(0:maxbrs1,minlam+1:maxsig)
      integer lbs2(0:maxbrs2,mincas+1:maxcas)
      integer lbm(0:maxbrm,minmes+1:maxmes)

      common /resonances/ massres,widres,massmes,widmes,mmesmn,
     ,                    branres,branmes,branbs1,branbs2,
     ,                    bs1type,bs2type,lbs1,lbs2,lbm,
     ,                    jres,jmes,lbr,brtype,pares,pames,
     ,                    bmtype,
     ,                    Isores,Isomes,strres,strmes,mlt2it
     
     
     
     
     
 
      !#############################################################################     
      !epos common blocks for particle list
      
      integer mmry,mxptl
      parameter (mmry=1)   !memory saving factor
      parameter (mxptl=200000/mmry) !max nr of particles in epos ptl list

      integer iorptl(mxptl),idptl(mxptl),istptl(mxptl),
     *  ifrptl(2,mxptl),jorptl(mxptl),ibptl(4,mxptl),ityptl(mxptl)
      real pptl(5,mxptl),tivptl(2,mxptl),xorptl(4,mxptl)

      common/cptl/nptl,pptl,iorptl,idptl
     *,istptl,tivptl,ifrptl,jorptl
     *,xorptl,ibptl,ityptl
   
      integer nptl
   
      integer nevt,kolevt,koievt,npjevt,ntgevt,npnevt,nppevt,
     *        ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt,
     *        nglevt,minfra,maxfra
      real phievt,bimevt,pmxevt,egyevt,xbjevt,qsqevt,zppevt,zptevt
   
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra

      integer laproj,maproj,latarg,matarg,nptlini,nptlbd
      real core,fctrmx
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      common/c4ptl/nptlbd
      
      integer istore,istmax,irescl,ntrymx,nclean,iopdg
      real gaumx
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg
      !#############################################################################     
      
      integer i,nn,iok,idtmp,ityptmp,iso3tmp,itmp
      


      integer idtrafo
      external idtrafo
      integer fchg
      external fchg
      real*8 dectim
      external dectim

      real*8 mintime,eb
      integer j,k,icount,npold
      integer strcount
      common /inewpartxx/ strcount
       
      call uinit(0)
      call osc_header
      call osc99_header

      npart = 0
      npold = 0
      nbar=0
      nmes=0
      uid_cnt=0
c reset counters
c     all collisions/decays
      ctag  = 0
c     all decays
      dectag = 0
c     number of prod. hard resonances
      NHardRes=0
c     number of prod. soft resonances
      NSoftRes=0
c     number of prod. resonances via decay
      NDecRes=0
c     number of blocked collisions
      NBlColl=0
c     number of elastic collisions
      NElColl=0
c     number of strings
      strcount=1
c
      eb=0D0
c icount is the number of EXTRAordinary pro/tar combinations (i.e. pion ...)
      icount = 0
c reset particle vectors
      do 20 j=1,nmax
	spin(j)  = 0
	ncoll(j) = 0
	lstcoll(j)=0
	r0(j) = 0.0
	rx(j)	 = 0.0
	ry(j)	 = 0.0
	rz(j)	 = 0.0
	p0(j)	 = 0.0
	px(j)	 = 0.0
	py(j)	 = 0.0
	pz(j)	 = 0.0
	frr0(j) = 0.0
	frrx(j)    = 0.0
	frry(j)    = 0.0
	frrz(j)    = 0.0
	frp0(j)    = 0.0
	frpx(j)    = 0.0
	frpy(j)    = 0.0
	frpz(j)    = 0.0
	fmass(j) = 0.0
	charge(j)= 0
	iso3(j)  = 0
	ityp(j)  = 0
	dectime(j)= 0.0
	origin(j)=0
	tform(j)=0.0
	xtotfac(j)=1.0
	strid(j)=0
	uid(j)=0
	 ffermpx(j) = 0.0
	 ffermpy(j) = 0.0
	 ffermpz(j) = 0.0

	 do 21 k=1,2
	    p0td(k,j)=0.d0
	    pxtd(k,j)=0.d0
	    pytd(k,j)=0.d0
	    pztd(k,j)=0.d0
	    fmasstd(k,j)=0.d0
	    ityptd(k,j)=0
	   iso3td(k,j)=0
 21	 continue
 20   continue


c epos event info to urqmd event info      
      bimp = bimevt
      
c convert epos list to urqmd list here. or in a seperate routine?
      npart=0  !number of output particles
      mintime = 1d2 !the minimum formation time
      nptlini=maproj+matarg+1
      
c fill in the baryons first      
      nbar = 0
      do nn=nptlini,nptlbd       
        iok=1
        if(istptl(nn).gt.istmax)iok=0
        if(iok.eq.1)then
           idtmp=idtrafo('nxs','pdg',idptl(nn))
           call pdg2id(ityptmp,iso3tmp,idtmp)
           if(abs(ityptmp).le.maxbar) then
               nbar=nbar+1
               r0(nbar)=xorptl(4,nn)
               rx(nbar)=xorptl(1,nn)
               ry(nbar)=xorptl(2,nn)
               rz(nbar)=xorptl(3,nn)
               p0(nbar)=pptl(4,nn)
               px(nbar)=pptl(1,nn)
               py(nbar)=pptl(2,nn)
               pz(nbar)=pptl(3,nn)
               fmass(nbar)=pptl(5,nn)
               ityp(nbar)=ityptmp
               iso3(nbar)=iso3tmp
               charge(nbar)=fchg(iso3(nbar),ityp(nbar))
               lstcoll(nbar)=0
               ncoll(nbar)=0
               origin(nbar)=0
               tform(nbar)=r0(nbar)
               dectime(nbar)=dectim(nbar,1)+tform(nbar)
               xtotfac(nbar)=0d0
               if(r0(nbar).lt.mintime) mintime = r0(nbar)
           endif
        endif
      enddo

c then fill in the mesons
      nmes = 0
      do nn=nptlini,nptlbd
        iok=1
        if(istptl(nn).gt.istmax)iok=0
        if(iok.eq.1)then
           idtmp=idtrafo('nxs','pdg',idptl(nn))
           call pdg2id(ityptmp,iso3tmp,idtmp)
           if(abs(ityptmp).ge.minmes) then
               nmes=nmes+1
               itmp=nbar+nmes
               r0(itmp)=xorptl(4,nn)
               rx(itmp)=xorptl(1,nn)
               ry(itmp)=xorptl(2,nn)
               rz(itmp)=xorptl(3,nn)
               p0(itmp)=pptl(4,nn)
               px(itmp)=pptl(1,nn)
               py(itmp)=pptl(2,nn)
               pz(itmp)=pptl(3,nn)
               fmass(itmp)=pptl(5,nn)
               ityp(itmp)=ityptmp
               iso3(itmp)=iso3tmp
               charge(itmp)=fchg(iso3(itmp),ityp(itmp))
               lstcoll(itmp)=0
               ncoll(itmp)=0
               origin(itmp)=0
               tform(itmp)=r0(itmp)
               dectime(itmp)=dectim(itmp,1)+tform(itmp)
               xtotfac(itmp)=0d0
               if(r0(itmp).lt.mintime) mintime = r0(itmp)
           endif
        endif
      enddo

      npart = nbar + nmes


c back to the same starting time
      do i = 1, npart
         !save freeze-out configuration, in case of no further
         !rescatterings
         frr0(i) = r0(i)
         frrx(i) = rx(i)
         frry(i) = ry(i)
         frrz(i) = rz(i)
         frp0(i) = p0(i)
         frpx(i) = px(i)
         frpy(i) = py(i)
         frpz(i) = pz(i)
         rx(i)=rx(i)-px(i)/p0(i)*(r0(i)-mintime)
         ry(i)=ry(i)-py(i)/p0(i)*(r0(i)-mintime)
         rz(i)=rz(i)-pz(i)/p0(i)*(r0(i)-mintime)
         r0(i)=mintime
      enddo
      
      acttime=mintime

c      write(*,*)'DEBUG INFO (epos.f): ',mintime,npart,istmax,nbar,nmes

      return
      end


c----------------------------------------------------------------------------------
      subroutine urqmdexit
c----------------------------------------------------------------------------------
c  transfer UrQMD -> EPOS
c----------------------------------------------------------------------------------

      implicit none


!urqmd common block

c     debug and validity range

      integer nmax, nspl
      real*8 hit_sphere
      parameter (nmax = 40000) ! maximum number of particles
      parameter (nspl = 500)  ! dimension of spline arrays
      parameter (hit_sphere = 8.d0)  ! hard collision cutoff: 251 mbarn

      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
c 7 integer

      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm



      real*8 
     +     gw, sgw, delr, fdel, dt,
     +     da, db,
     +     Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, gamYuk, drPau, dpPau,
     +     dtimestep
c 19 real*8
     
      real*8 cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      common /cuts/ cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww

      real*8 spx(nspl), spPauy(nspl), outPau(nspl), 
     +                spCby(nspl),  outCb(nspl),
     +                spYuky(nspl), outYuk(nspl),
     +                spSkyy(nspl), outSky(nspl),
     +                spdwwy(nspl), outdww(nspl)

      common /spdata/ spx, spPauy, outPau, spCby,  outCb,
     +                     spYuky, outYuk, spSkyy, outSky,
     +                     spdwwy, outdww

      real*8 
     +     r0(nmax), rx(nmax), ry(nmax), rz(nmax),
     +     p0(nmax), px(nmax), py(nmax), pz(nmax),
     +     airx(nmax), airy(nmax), airz(nmax),
     +     aipx(nmax), aipy(nmax), aipz(nmax),
     +     aorx(nmax,4), aory(nmax,4), aorz(nmax,4),
     +     aopx(nmax,4), aopy(nmax,4), aopz(nmax,4),
     +     fmass(nmax), rww(nmax), 
     +     dectime(nmax), tform(nmax), xtotfac(nmax)
      
      
      integer spin(nmax),ncoll(nmax),charge(nmax),strid(nmax),
     +        ityp(nmax),lstcoll(nmax),iso3(nmax),origin(nmax),uid(nmax)
      common/isys/spin,ncoll,charge,ityp,lstcoll,iso3,origin,strid,
     +            uid
     
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime
      common /frag/ tform, xtotfac


      common /aios/ airx, airy, airz, aipx, aipy, aipz,
     +              aorx, aory, aorz, aopx, aopy, aopz

      common /pots/ Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, 
     +              gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
     +              dt,da, db,dtimestep


c spectator arrays
	integer smax
	parameter(smax=500)  ! maximum number of spectators
	real*8 r0s(smax), rxs(smax), rys(smax), rzs(smax),
     +	       p0s(smax), pxs(smax), pys(smax), pzs(smax),
     +	       sfmass(smax)
	

        integer sspin(smax), scharge(smax), sityp(smax), siso3(smax),
     +          suid(smax)

	integer nspec

	common /scoor/ r0s, rxs, rys, rzs, p0s, pxs ,pys, pzs, sfmass

	common /sisys/ sspin, scharge, sityp, siso3, suid

	common /ssys/ nspec


        real*8 p0td(2,nmax),pxtd(2,nmax),pytd(2,nmax),pztd(2,nmax),
     +         fmasstd(2,nmax)
        integer ityptd(2,nmax),iso3td(2,nmax)
        integer itypt(2),uidt(2),origint(2),iso3t(2)

        common /rtdelay/p0td,pxtd,pytd,pztd,fmasstd
        common /itdelay/ityptd,iso3td
        common /svinfo/itypt,uidt,origint,iso3t
        real*8 ffermpx(nmax), ffermpy(nmax), ffermpz(nmax)
        real*8 peq1, peq2
        common /ffermi/ ffermpx, ffermpy, ffermpz
        common /peq/ peq1,peq2
        
	
	integer numcto,numctp,maxstables
        parameter(numcto=400) ! maximum number of options
        parameter(numctp=400) ! maximum number of parameters
        parameter(maxstables=20) ! maximum number of stable particles
	integer   CTOption(numcto)
	real*8    CTParam(numctp)
        common /options/CTOption,CTParam
	
	real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     +     frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)

      common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz 
	      integer maxbar,maxbra,minbar
      integer offmeson,maxmeson,pimeson,maxbrm,minnuc,mindel
      integer maxbrs1,maxbrs2
      integer numnuc,numdel,nucleon,maxnuc,maxdel
      integer minmes,maxmes


      parameter (minnuc=1) ! lowest baryon particle ID 
      parameter (minmes=100) ! lowest meson particle ID
      parameter (maxmes=132) ! hightest meson particle ID

c number of resonances of a kind
      parameter (numnuc=16) ! number of nucleon resonances
      parameter (numdel=10) ! number of delta resonances
c indices of minimal and maximal itype of a kind (redundant but nice)
      parameter (maxnuc=minnuc+numnuc-1) ! highest nucleon ID
      parameter (mindel=minnuc+maxnuc)   ! lowest delta ID
      parameter (maxdel=mindel+numdel-1) ! highest delta ID

c minres & maxres define the range of nonstable & nonstrange baryons
      integer minres,maxres
      parameter (minres=minnuc+1) ! lowest baryon resonance ID
      parameter (maxres=maxdel)   ! highest (nonstrange) baryon 
                                  ! resonance ID

c strangenes.ne.0 baryon resonances
      integer minlam,minsig,mincas,minome
      integer numlam,numsig,numcas,numome
      integer maxlam,maxsig,maxcas,maxome
      parameter (numlam=13) ! number of lambda states
      parameter (numsig=9)  ! number of sigma states
      parameter (numcas=6)  ! number of cascade states
      parameter (numome=1)  ! number of omega states
      parameter (minlam=mindel+numdel)   ! ID of lowest lambda state
      parameter (maxlam=minlam+numlam-1) ! ID of highest lambda state
      parameter (minsig=minlam+numlam)   ! ID of lowest sigma state
      parameter (maxsig=minsig+numsig-1) ! ID of highest sigma state
      parameter (mincas=minsig+numsig)   ! ID of lowest cascade state
      parameter (maxcas=mincas+numcas-1) ! ID of highest cascade state
      parameter (minome=mincas+numcas)   ! ID of lowest omega state
      parameter (maxome=minome+numome-1) ! ID of highest omega state

c minbar & maxbar define the range of all baryons
      parameter (minbar=minnuc) ! ID of lowest baryon state
      parameter (maxbar=maxome) ! ID of highest baryon state

      parameter (offmeson=minmes) ! offset between zero and lowest 
                                  ! meson state
      parameter (maxmeson=maxmes) ! ID of highest meson state
c... these variables are in principal obsolete and should be exchanged 
c were referenced 

c... avoid hard coded itypes
      integer itrho,itome,iteta,itkaon,itphi,itetapr
      parameter (itkaon=106)   ! ID of kaon
      parameter (itrho=104)    ! ID of rho meson 
      parameter (itome=103)    ! ID of omega meson
      parameter (iteta=102)    ! ID of eta
      parameter (itphi=109)    ! ID of phi
      parameter (itetapr=107)  ! ID of eta'
      parameter (pimeson=101)  ! ID of $\pi$
      parameter (nucleon=minnuc) ! ID of nucleon

      integer itmin,itmax
      parameter (itmin=minnuc)  ! lowest defined ID
      parameter (itmax=maxmes)  ! highest defined ID
c
      parameter (maxbra=11)  ! decay channels for $s=0$ baryon resonances
      parameter (maxbrm=25) ! decay channels for meson resonances
      parameter (maxbrs1=10)! decay channels for $s=1$ baryon resonances
      parameter (maxbrs2=3) ! decay channels for $s=2$ baryon resonances

c 
       integer mlt2it(maxmes-minmes) ! meson IDs sorted by multipletts


      real*8 massoff,mresmin,mresmax
      parameter (massoff=1d-4)      ! offset for mass generation
      parameter (mresmin=1.0765d0)  ! minimum baryon resonance mass
      parameter (mresmax=5d0)       ! maximum baryon resonance mass

      character*45 versiontag
      common /versioning/ versiontag

      real*8 massres(minbar:maxbar),widres(minbar:maxbar)
      real*8 branmes(0:maxbrm,minmes+1:maxmes)
      real*8 branres(0:maxbra,minnuc+1:maxdel)
      real*8 branbs1(0:maxbrs1,minlam+1:maxsig)
      real*8 branbs2(0:maxbrs2,mincas+1:maxcas)
      integer Jres(minbar:maxbar)
      integer Jmes(minmes:maxmes)
      integer pares(minbar:maxbar),pames(minmes:maxmes)
      integer Isores(minbar:maxbar), Isomes(minmes:maxmes)
      integer brtype(4,0:maxbra),bmtype(4,0:maxbrm)
      integer bs1type(4,0:maxbrs1),bs2type(4,0:maxbrs2)
      real*8 massmes(minmes:maxmes)
      real*8 mmesmn(minmes:maxmes)
      real*8 widmes(minmes:maxmes)
      integer strres(minbar:maxbar),strmes(minmes:maxmes)

      integer lbr(0:maxbra,minnuc+1:maxdel)
      integer lbs1(0:maxbrs1,minlam+1:maxsig)
      integer lbs2(0:maxbrs2,mincas+1:maxcas)
      integer lbm(0:maxbrm,minmes+1:maxmes)

      common /resonances/ massres,widres,massmes,widmes,mmesmn,
     ,                    branres,branmes,branbs1,branbs2,
     ,                    bs1type,bs2type,lbs1,lbs2,lbm,
     ,                    jres,jmes,lbr,brtype,pares,pames,
     ,                    bmtype,
     ,                    Isores,Isomes,strres,strmes,mlt2it

      
      !epos common blocks for particle list
      integer mmry,mxptl
      parameter (mmry=1)   !memory saving factor
      parameter (mxptl=200000/mmry) !max nr of particles in epos ptl list

      integer iorptl(mxptl),idptl(mxptl),istptl(mxptl),
     *  ifrptl(2,mxptl),jorptl(mxptl),ibptl(4,mxptl),ityptl(mxptl)
      real pptl(5,mxptl),tivptl(2,mxptl),xorptl(4,mxptl)

      common/cptl/nptl,pptl,iorptl,idptl
     *,istptl,tivptl,ifrptl,jorptl
     *,xorptl,ibptl,ityptl
     
      
      integer nptl
      
      
      integer nevt,kolevt,koievt,npjevt,ntgevt,npnevt,nppevt,
     *        ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt,
     *        nglevt,minfra,maxfra
      real phievt,bimevt,pmxevt,egyevt,xbjevt,qsqevt,zppevt,zptevt
      
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra
      
      integer istore,istmax,irescl,ntrymx,nclean,iopdg
      real gaumx
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg
      
      integer nn,idpdgg,idepos
      


      integer idtrafo
      external idtrafo
      integer fchg
      external fchg
      integer pdgid
      external pdgid
      real*8 dectim
      external dectim

      nbar=0
      do nn=1,npart  
csp        write(6,*) 'ityp', ityp(nn), 'iso3', iso3(nn)     
	   idpdgg=pdgid(ityp(nn),iso3(nn))
           idepos=idtrafo('pdg','nxs',idpdgg)
              
	       nbar=nbar+1
               
               xorptl(4,nbar)= r0(nn)
               xorptl(1,nbar)=rx(nn)
               xorptl(2,nbar)=ry(nn)
               xorptl(3,nbar)=rz(nn)
               pptl(4,nbar)= p0(nn)
               pptl(1,nbar)=px(nn)
               pptl(2,nbar)=py(nn)
               pptl(3,nbar)=pz(nn)
               pptl(5,nbar)=fmass(nn)
               idptl(nbar)=idepos
               
                                  
	nptl=nbar
      enddo
      
      end

c----------------------------------------------------------------------------------
      subroutine  urqmdepos
c----------------------------------------------------------------------------------
c   main urqmd modul, essentially copied from urqmdepos (urqmd.f)
c----------------------------------------------------------------------------------

      implicit none
      include '../urqmd23/coms.f'
      include '../urqmd23/comres.f'
      include '../urqmd23/options.f'
      include '../urqmd23/colltab.f'
      include '../urqmd23/inputs.f'
      include '../urqmd23/newpart.f'
      include '../urqmd23/boxinc.f'
      integer i,j,k,steps,ii,ocharge,ncharge, mc, mp, noc, it1,it2
      real*8 sqrts,otime,xdummy,st
      logical isstable
      integer stidx
      real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      integer cti1sav,cti2sav

      mc=0
      mp=0
      noc=0

      time = 0.0  !time is the system time at the BEGINNING of every timestep

      !initialize random number generator
      !call auto-seed generator only for first event and if no seed was fixed
      if(.not.firstseed.and.(.not.fixedseed)) then
         ranseed=-(1*abs(ranseed))
         call sseed(ranseed)
      else
         firstseed=.false.
      endif

      !old time if an old fort.14 is used 
      if(CTOption(40).eq.1)time=acttime
      if(CTOption(40).eq.3)time=acttime

      !write headers to file
      call output(13)
      !call output(14)
      call output(15)
      call output(16)
      !if(event.eq.1)call output(17)
      call osc99_event(-1)

      !for CTOption(4)=1 : output of initialization configuration
      if(CTOption(4).eq.1)call file14out(0)
      !participant/spectator model:
      if(CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

      otime = outsteps*dtimestep  !compute time of output

      steps = 0  !reset time step counter

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! loop over all timesteps
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      do 20  steps=1,nsteps  

         if (eos.ne.0) then
           do j=1,npart
               r0_t(j) = r0(j)
               rx_t(j) = rx(j)
               ry_t(j) = ry(j)
               rz_t(j) = rz(j)
           enddo
         end if

         !we are at the beginning of the timestep, set current time (acttime)
         acttime = time

         if(CTOption(16).ne.0) goto 103  !option for MD without collision term

         call colload  ! Load collision table with next collisions in current timestep

         ! check for collisions in time-step, nct = # of collisions in table
         if (nct.gt.0) then

 101        continue              !entry-point for collision loop in case  
            k = 0                 !      of full colload after every coll
 100        continue              !normal entry-point for collision loop 
            call getnext(k)       !get next collision
            if (k.eq.0) goto 102  !exit collision loop if no collisions are left

            !propagate all particles to next collision time
            !store actual time in acttime, propagation time st=cttime(k)-acttime
	    st=cttime(k)-acttime
            call cascstep(acttime,st)
            acttime = cttime(k)   !new actual time (for upcoming collision)

            !perform collision 

            if(cti2(k).gt.0.and.
     .           abs(sqrts(cti1(k),cti2(k))-ctsqrts(k)).gt.1d-3)then
               write(6,*)' ***(E) wrong collision update (col) ***'
               write(6,*)cti1(k),cti2(k),
     .              ctsqrts(k),sqrts(cti1(k),cti2(k))
            else if(cti2(k).eq.0.and.
     .              abs(fmass(cti1(k))-ctsqrts(k)).gt.1d-3) then
               write(6,*)' *** main(W) wrong collision update (decay)'
               write(6,*)ctag,cti1(k),ityp(cti1(k)),dectime(cti1(k)),
     .              fmass(cti1(k)),ctsqrts(k)
            endif

            ocharge=charge(cti1(k))
            if(cti2(k).gt.0) ocharge=ocharge+charge(cti2(k))

            !store quantities in local variables for charge conservation check
            it1= ityp(cti1(k))
            if(cti2(k).gt.0)it2= ityp(cti2(k))

            !increment "dirty" collision counter
            if(cti2(k).gt.0)then !scatter
               mc=mc+1
            endif
            !perform scattering/decay
            cti1sav = cti1(k)               
            cti2sav = cti2(k)    
            call scatter(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k),
     &                   ctcolfluc(k))

            !update collision table 

            !normal update mode
            if(CTOption(17).eq.0) then
               if(nexit.eq.0) then
                 !new collision partners for pauli-blocked states (nexit=0)
                 if (cti1(k).ne.cti1sav.or.cti2(k).ne.cti2sav) then
                   cti1(k) = cti1sav 
                   cti2(k) = cti2sav 
                 endif

                 call collupd(cti1(k),1)
                 if(cti2(k).gt.0) call collupd(cti2(k),1)
               else
                 ncharge=0
                 !new collision partners for scattered/produced particles (nexit><0)
                 do i=1,nexit
                   !ncharge is used for charge conservation check
                   ncharge=ncharge+charge(inew(i))
                   call collupd(inew(i),1)
                 enddo
               endif
               !update collisions for partners of annihilated particles
               do ii=1,nsav
                  call collupd(ctsav(ii),1)
               enddo
               nsav=0
            else ! (CTOption(17).ne.0)
              !full collision load
              call colload
            endif

            if (CTOption(17).eq.0) goto 100
            goto 101

            !this is the point to jump to after all collisions in the timestep
            !have been taken care of
 102        continue

         endif ! (nct.gt.0)


         !After all collisions in the timestep are done, propagate to end of 
         !the timestep.

         !point to jump to in case of MD without collision term
 103     continue

         time = time+dtimestep  !increment timestep

         !After all collisions in the timestep are done, propagate to end of 
         !the timestep.
         call cascstep(acttime,time-acttime)

         !in case of potential interaction, do MD propagation step
         if (eos.ne.0) then
            ! set initial conditions for MD propagation-step
            do j=1,npart
               r0(j) = r0_t(j)
               rx(j) = rx_t(j)
               ry(j) = ry_t(j)
               rz(j) = rz_t(j)
            enddo
            !now molecular dynamics trajectories
            call proprk(time,dtimestep)
         endif ! (eos.ne.0)

         !perform output if desired
         if(mod(steps,outsteps).eq.0.and.steps.lt.nsteps)then 
            if(CTOption(28).eq.2)call spectrans(otime)
            call file14outx(steps)
         endif ! output handling

 20   continue ! time step loop

      acttime=time
      
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! optional decay of all unstable 
  !  particles before final output
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      !DANGER: pauli-blocked decays are not performed !!!
      if(CTOption(18).eq.0.and.CTOption(51).eq.0) then
         !no do-loop is used because npart changes in loop-structure
         i=0
         nct=0
         actcol=0
         CTOption(10)=1  !disable Pauli-Blocker for final decays
 40      continue  !decay loop structure starts here
         i=i+1
         if(dectime(i).lt.1.d30) then !if particle unstable
 41         continue
            isstable = .false.
            do stidx=1,nstable
               if (ityp(i).eq.stabvec(stidx)) then
                  !write (6,*) 'no decay of particle ',ityp(i)
                  isstable = .true.
               endif
            enddo
            if (.not.isstable) then
               call scatter(i,0,0.d0,fmass(i),xdummy) !perform decay
               !backtracing if decay-product is unstable itself
               if(dectime(i).lt.1.d30) goto 41
            endif
         endif
         !check next particle
         if(i.lt.npart) goto 40
      endif ! final decay

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !     final output
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(CTOption(28).eq.2)call spectrans(otime)

      call file13out(nsteps)
      !call file14out(nsteps)
      call file16out
      call osc_event
      call osc99_event(1)
      call osc99_eoe
      
      mp=mp+npart
      if(ctag.eq.0)then
         write(*,*)'(W) No collision in event ',event
         noc=noc+1
      endif

      end

c-------------------------------------------------------------------------------------
      subroutine file14outx(timestep)
c-------------------------------------------------------------------------------------
      implicit none
      include '../urqmd23/comres.f'
      include '../urqmd23/coms.f'
      include '../urqmd23/options.f'
      include '../urqmd23/inputs.f'
      include '../urqmd23/newpart.f'
      include '../urqmd23/freezeout.f'
      include '../urqmd23/boxinc.f'
      integer i,ttime,timestep,itotcoll,iinelcoll
      real*8 sigmatot,t
      common /outco2/sigmatot
      include '../urqmd23/outcom.f'
      save 

      if(bf14)return
      ttime=int(timestep*dtimestep+0.01)
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl
      print*,'(file14outx)',ttime,npart
     @ ,itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes
       !----------------------------------------------------------------------------
       !  r0(i), rx(i), ry(i), rz(i)   ........................................ x4
       !  p0(i),px(i)+ffermpx(i),py(i)+ffermpy(i),pz(i)+ffermpz(i),fmass(i) ... p5
       !  ityp(i)  ..... particle id (listed in the urqmd user guide)
       !  iso3(i) ...... 2 times the isospin of a particle
       !  charge(i) .... charge of the particle
       !  lstcoll(i) ... index of last collision partner
       !  ncoll(i) ..... number of collisions
       !  origin(i) ....
       !  dectime(i) ...
       !  tform(i) ..... formation time
       !  xtotfac(i) ... cross section (zero if the particle is not yet formed)
       !  uid(i) ......
       !---------------------------------------------------------------------------
       do i=1,npart
         t=r0(i)
       enddo

      end

