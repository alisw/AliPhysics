
	subroutine hiwlnew(ng,mg,hidmin,hidmax,prb,mmax)
	integer mg(*)
	include 'comalwl.for'
	include 'comgam.for'
	data wlmin /0.02/

        hi=0.
	himax=0.
	hidmax=0.
	hidmin=1000000.
	if(ng.eq.0) return
        nb=0.
	do n=1,ng
	m=mg(n)
	if(wl(m).gt.wlmin) then
        nb=nb+1
c        write (*,*)  'swlt(m),swlt0(m)',swlt(m),swlt0(m)
	dspe=swlt(m)+swlt0(m)
	de=(wl(m)-wlt(m))**2
	dde=wl(m)-wlt(m)
c fix by Federico Carminati CERN-EP, just to avoid a crash
c it should be checked by the authors
	dspe=max(dspe,1e-11)
	hig=de/dspe
	hidg=dde/sqrt(dspe)
c	if(wlt(m).eq.0) write (*,*) ' n,m,wl(m),hi ',n,m,wl(m),sqrt(hig)
	hi=hi+hig

c	if(hig.gt.himax) then
c		himax=hig
c		mmax=m
c	endif

	if(hidg.gt.hidmax) then
		hidmax=hidg
		mmax=m
	endif
	if(hidg.lt.hidmin) then
		hidmin=hidg
	endif

	endif
	enddo
        nbp=nb-3*kg
c        nbp=nb

        if(nbp.gt.0) then
                aln=alog(float(nbp))
c               himax0=-0.088+0.548*aln-0.01396*aln**2

                hidmax0=0.318+0.5942*aln-0.0258*aln**2
                hidmax=hidmax-hidmax0
                prb=prob(hi,nbp)
        else
                himax0=0.
                hidmax0=0.
                prb=1.
        endif
c        write (*,*) ' kg,nb, nbp,prb ',kg,nb,nbp,prb

c	hi=sqrt(h)
c	himax=sqrt(himax)
	return
	end
	subroutine fitgaml(ng,mg,eg,hi,himax,mmax,level)
	dimension mg(*),eg(*)
	if(level.eq.-1) then
		call fitgam0(ng,mg,eg,hi,himax,mmax)
	elseif(level.eq.0) then
		call fitgam1(ng,mg,eg,hi,himax,mmax)
	elseif(level.eq.1) then
c		call fitgam2(ng,mg,eg,hi,himax,mmax)
		write (*,*) ' next time '
        else
                write (*,*) ' oshibsja parenek? '
	endif
	return
	end

	subroutine fitgam0(ng,mg,eg,hi,himax,mmax)
        common /error/ier
	include 'comalwl.for'
	include 'comgam.for'
	include 'comwlgen.for'
	dimension mg(*),eg(*)

c simplgam
	call filet
	call filwt
	call diven
	call ecgam
	call clrwt

	call filet
	call filwt
	call diven
	call ecgam
	call clrwt

	call filet
	call filwt
	call hiwlnew(ng,mg,himin,himax,prb,mmax)
	call clrwt

	return
	end

	subroutine fitgam1(ng,mg,eg,prb,himax,mmax)
        common /error/ier
	include 'comalwl.for'
	include 'comgam.for'
	include 'comwlgen.for'
	dimension mg(*),eg(*)

c	call fitgam(ng,mg,eg)
	call again_gmfit
	call filet
	call filwt
	call hiwlnew(ng,mg,himin,himax,prb,mmax)
	call clrwt
	return
	end

	subroutine gmplus(ng,mg,eg,gmporog,egamporog,rnenear)
	common /error/ier

	include 'comgl.for'
	include 'comgam.for'
	include 'comalwl.for'
	include 'comwlgen.for'
	common /comwldif/ swldiff(40000),eadd(2000)

	dimension mg(*),eg(*)
	data dmax /27./
c	data gmporog/0.1/
	integer ijdd(2)

c
	call filwt
	if(ng.gt.2000) then
		write (*,*) ' incr.in gmplus '
		return
	endif

	do n=1,ng
	m=mg(n)
	sig=sqrt(swlt(m)+swlt0(m))
		efar=0.
		do 1 kk=1,kg
		dx=abs(x(1,mw(kk))-x(1,m))-(d(1,mw(kk))+d(1,m))/2.
		dy=abs(x(2,mw(kk))-x(2,m))-(d(2,mw(kk))+d(2,m))/2.
		if(dx.le.1..and.dy.le.1.) goto 1
		if(dx.gt.dmax.or.dy.gt.dmax) goto 1
		xd=xw(kk)-x(1,m)
		yd=yw(kk)-x(2,m)
		adefar=e(kk)*aclnew(xd,yd,d(1,m),m)
		efar=efar+adefar
1		continue
	sefar=0.15*sqrt(efar)
	swldiff(m)=wl(m)-wlt(m)-efar-sqrt(swlt(m)+swlt0(m))-sefar
	if(swldiff(m).le.0.) swldiff(m)=0.
	enddo
	call clrwt
c
c  LET'S TRY TO FINDE NEW GAMMA

	do n=1,ng
	em=eg(n)
	m=mg(n)
	adsum=swldiff(m)
	do k=1,kg
	if(m.eq.mw(k)) goto 2
	rrk=sqrt((x(1,m)-xw(k))**2+(x(2,m)-yw(k))**2)
	if(rrk.lt.rnenear) goto 2
	enddo

	if(em.lt.gmporog) goto 2
	esumgg1=em
	esumgg2=em
	do is=2,ns(m)
	ms=ks(is,m)
	esumgg1=esumgg1+swldiff(ms)
	esumgg2=esumgg2+wl(ms)
	enddo
	if(esumgg1.lt.egamporog) goto 2

		do  3 n1=1,ng
		if(n.ne.n1) then
		m1=mg(n1)
		em1=swldiff(m1)
		dx=abs(x(1,m1)-x(1,m))-(d(1,m)+d(1,m1))/2.
		dy=abs(x(2,m1)-x(2,m))-(d(2,m)+d(2,m1))/2.
		if(dx.le.1..and.dy.le.1.) then
			if(em.gt.em1) then
				adsum=adsum+em1
			endif
			goto 3
		endif
		if(dx.gt.dmax.or.dy.gt.dmax) goto 3
		xd=dx+d(1,m)/2.
		yd=dy+d(2,m)/2.
		if(xd.le.0.) xd=0.
		if(yd.le.0.) yd=0.
		emax=1.5*em*aclnew(xd,yd,d(1,m),m)
		if(emax.gt.em1) then
			adsum=adsum+em1
		endif
		endif
3		continue
2	eadd(n)=adsum
	enddo

	nmmm=1
	summax=eadd(nmmm)
	mmax=mg(nmmm)
	emax=eg(nmmm)

	do n=2,ng
	if(eadd(n).gt.summax) then
		nmmm=n
		summax=eadd(nmmm)
		mmax=mg(nmmm)
		emax=eg(nmmm)
	endif
	if(eadd(n).eq.summax) then
	if(eg(n).gt.emax) then
		nmmm=n
		summax=eadd(nmmm)
		mmax=mg(nmmm)
		emax=eg(nmmm)
	endif
	endif
	enddo

c	write (*,*) ' wl '
c	call prwall(wl)
c	write (*,*) ' wlt '
c	call prwall(wlt)
c	write (*,*) ' swldiff '
c	call prwall(swldiff)

	if(emax.gt.gmporog.and.summax.gt.0.) then
c		kg=kg+1
c		e(kg)=emax+summax-swldiff(mmax)
c		mw(kg)=mmax
		ee=emax+summax-swldiff(mmax)
		xx=x(1,mmax)
		yy=x(2,mmax)
		see=ee/4.
		sxx=d(1,mmax)
		syy=d(2,mmax)
		ijdd(1)=1
		ijdd(2)=1
	call adgmnew(mmax,ee,xx,yy,see,sxx,syy,ijdd)

c		write (*,*) ' NEW GAMMA e,m',kg,e(kg),mw(kg)
c		write(*,*) ' swldiff kg,mg',kg,mw(1)
c		call prwlm(swldiff,mw(1),5,100.)
	else
c		write (*,*) ' NO NEW GAMMA ',kg
c		write(*,*) ' swldiff kg,mg',kg,mg(1)
c		call prwlm(swldiff,mw(1),5,100.)
	endif

	call vzero(swldiff,nt)

	return
	end
		
	subroutine filwltgl(ngg,mgg)
	include 'comalwl.for'
	include 'comgam.for'
	include 'comwlgen.for'
        common /comdevpar/sigmaph,sigmapd,sigphsq,sigpdsq

        integer mgg(*)

        if(ngg.eq.0) return

	do k=1,kg
	do n=1,ngg	
	mg=mgg(n)

	if(mg.gt.0) then
	etadd=e(k)*ampcelnew(e(k),raxay(1,k),x(1,mg),d(1,mg),mg)
	wlt(mg)=wlt(mg) +etadd
c	swlt(ms)=swlt(ms)+(emimx(2,n,k)-emimx(1,n,k))**2/4.
	swlt(mg)=swlt(mg)+
     +  (sigwlgam0(etadd,e(k)))**2+sigmphsq*etadd
	endif

	enddo
	enddo

	return
	end

	subroutine clrwltgl(ngg,mgg)
	include 'comalwl.for'
	include 'comgam.for'
	include 'comwlgen.for'

        integer mgg(*)

        if(ngg.eq.0) return

	do n=1,ngg	
	mg=mgg(n)
	if(mg.gt.0) then
	wlt(mg)=0.
	swlt(mg)=0.
	endif
	enddo

	return
	end

	subroutine again_gmfit
	common /error/ier
	include 'comwlgen.for'
	include 'comalwl.for'
	include 'comgam.for'
	common /par/npar,e0,x0,y0,se0,sx0,sy0,Hi,rrr0(5)
	real ws(100),pinc(4)


c	call deriv_test

	call filet
	call filwt
	call filwsnew
	do k=1,kg
	m=mw(k)
	e0=e(k)

c	x0=x(1,m)+id(k)*d(1,k)/4.
c	y0=x(2,m)+jd(k)*d(2,k)/4.
c	write (*,*) ' k,xm,ym,id,jd ',x(1,m),x(2,m),id(k),jd(k)

	x0=xw(k)
	y0=yw(k)

c only for anglx and angly
	call ucopy(raxay(1,k),rrr0(1),5)	
	call ucopy(wsw(1,k),ws,ns(m))

c	write (*,*) ' es k',k,(es(ii,k),ii=1,ns(m))
c	write (*,*) ' ws k',k,(sqrt(1./wsw(ii,k)),ii=1,ns(m))

c	write (*,*) ' e0,x0,y0 before ',e0,x0,y0

	se0=0.
	sx0=0.
	sy0=0.

	call gmfit(m,es(1,k),ws)

c	write (*,*) ' e0,x0,y0 after ',e0,x0,y0
c	write (*,*) ' se0,sx0,sy0    ',se0,sx0,sy0


        e(k)=e0
	xw(k)=x0
	yw(k)=y0
	raxay(1,k)=x0
	raxay(2,k)=y0

	sigexy(1,k)=se0
	sigexy(2,k)=sx0
	sigexy(3,k)=sy0


	enddo
	call clrwt

	return
	end	

	subroutine gmfit(m,eg,weg)
cold	external acl,dxacl	
cnew
	external ampcelnew,dampcelnew
	include 'comalwl.for'
	common /error/ier
	common /par/npar,e0,x0,y0,se0,sx0,sy0,Hi,rrr0(5)
	common /comkey/keykey
	real eg(25),weg(25)
cnew	
c	real rrr0(5)

	real f(4),ff(4,4),de(5,5)
	equivalence (aij,f(1)),(dxij,f(2)),(dyij,f(3)),(deij,f(4))
	equivalence (ade,ff(1,4)),(dxde,ff(2,4)),(dyde,ff(3,4))
	equivalence (aq,ff(1,1)),(dxq,ff(2,2)),(dyq,ff(3,3))
	equivalence (adx,ff(1,2)),(ady,ff(1,3)),(dxdy,ff(2,3))
	iefc=0
	ixfc=0
	iyfc=0	
	nccmx=100
c	write (*,*) ' se0,sx0,sy0 ',se0,sx0,sy0

	if(se0.gt.0.) then
		iefc=1
		ec=e0
		wec=1./se0**2
	endif
	if(sx0.gt.0.) then
		ixfc=1
		xc=x0
		wxc=1./sx0**2
	endif
	if(sy0.gt.0.) then
		iyfc=1
		yc=y0
		wyc=1./sy0**2
	endif
	hiold=10000000.
	do 10 ncc=1,nccmx

	call vzero(f,4)
	call vzero(ff,16)
	hi=0.	
	ssww=0.
        DO N=1,NS(M)
        MS=KS(N,M)
        XS=X(1,MS)
        YS=x(2,MS)
cold        aij=ACL(X0-XS,Y0-YS,D(1,ms))
cnew
	rrr0(1)=x0
	rrr0(2)=y0
c	rrr0(3)=z0
c	write (*,*) ' ax0,ay0 ',rrr0(4),rrr0(5)
c	rrr0(4)=ax0
c	rrr0(5)=ay0
	aij=ampcelnew(e0,rrr0,x(1,ms),d(1,ms),ms)
cnewend
	deij=eg(n)-e0*aij
cold	dxij=dxacl(x0-xs,y0-ys,d(1,ms))
cold	dyij=dxacl(y0-ys,x0-xs,d(1,ms))

cnew
	dxij=dampcelnew(1,e0,rrr0,x(1,ms),d(1,ms),ms)
	dyij=dampcelnew(2,e0,rrr0,x(1,ms),d(1,ms),ms)
cnewend
c	if(keykey.eq.1) then
c
c	write (*,*) ' deij,dxij,dyij ',deij,dxij,dyij
c
c	endif
	do 2 k=1,4
	do 2 l=1,4
	ff(k,l)=ff(k,l)+f(k)*f(l)*weg(n)
2	continue
	ssww=ssww+weg(n)
        ENDDO
	if(iefc.ne.0) then
		ade=ade+wec*(ec-e0)
		aq=aq+wec
	endif
	if(ixfc.ne.0) then
		dxde=dxde+wxc*(xc-x0)/e0
		dxq=dxq+wxc/e0**2
	endif
	if(iyfc.ne.0) then
		dyde=dyde+wyc*(yc-y0)/e0
		dyq=dyq+wyc/e0**2
	endif
	d11=dxq*dyq-dxdy**2
	d12=ady*dxdy-adx*dyq
	d13=adx*dxdy-ady*dxq
	d22=aq*dyq-ady**2
	d23=adx*ady-aq*dxdy
	d33=aq*dxq-adx**2
	det=aq*d11+adx*d12+ady*d13
	if(e0.le.0.) then
	write (*,*) ' e0=0 !!!!!!!!!!!!!!!!!!'
		e0=0.
		return
	endif
	if(det.eq.0.) then 
	write (*,*) ' det=0 !!!!!!!!!!!!!!!!!!'
		e0=0.
		ier=1
		return
	endif
	de0=(ade*d11+dxde*d12+dyde*d13)/det
	dx0=(ade*d12+dxde*d22+dyde*d23)/(e0*det)
	dy0=(ade*d13+dxde*d23+dyde*d33)/(e0*det)


	if(abs(de0).gt.3.) de0=3.*de0/abs(de0)
	if(abs(dx0).gt.d(1,ms)/3.) dx0=d(1,ms)/3.*dx0/abs(dx0)
	if(abs(dy0).gt.d(2,ms)/3.) dy0=d(2,ms)/3.*dy0/abs(dy0)
c	write (*,*) ' de0,dx0,dy0 ',de0,dx0,dy0
	e0=e0+de0
	x0=x0+dx0
	y0=y0+dy0
c e=0!!!!!
	if(e0.le.0.03) then
c	write (*,*) ' return e0 !!!',e0
c	write (*,*) ' main ',m,ns(m),x(1,m),x(2,m)
c	write (*,*) ' eg ',(eg(ne),ne=1,ns(m))
c	write (*,*) ' weg ',(sqrt(1./weg(ne)),ne=1,ns(m))
	e0fit=e0
	e0=0.
        DO Ne=1,NS(M)
	e0=e0+eg(ne)
        ENDDO
	if(e0.gt.1000.) write (*,*) ' return e0 from ',e0fit,' to',e0
c	write (*,*) 'nc,e0,x0,y0',ncc,e0,x0,y0	
	endif
	Hi=ff(4,4)
	if(iefc.ne.0) Hi=Hi+wec*(ec-e0)**2	
	if(ixfc.ne.0) Hi=Hi+wxc*(xc-x0)**2	
	if(iyfc.ne.0) Hi=Hi+wyc*(yc-y0)**2	
	hi=sqrt(hi)
	delthi=(hi-hiold)/hi
c	write (*,*) '       Hi delthi hiold',hi,delthi,hiold
	hiold=hi
	if(abs(delthi).lt.0.01) goto 11
10	continue
11	continue
c	write (*,*) ' Hi delthi ',hi,delthi
	if(d11/det.gt.0.) se0=sqrt(d11/det)
	if(d22/det.gt.0.) sx0=sqrt(d22/det/e0)
	if(d33/det.gt.0.) sy0=sqrt(d33/det/e0)
	if(se0.eq.0) write (*,*) ' se0=0.'
	if(sx0.eq.0) write (*,*) ' sx0=0.'
	if(sy0.eq.0) write (*,*) ' sy0=0.'
c	write (*,*) ' se0,sx0,sy0 ',se0,sx0,sy0
	return
	end

	subroutine filwsnew
	include 'comalwl.for'
	include 'comgam.for'
	do k=1,kg
	m=mw(k)
c	write (*,*) ' new gamma print '
        DO N=1,NS(M)
        MS=KS(N,M)
c???	wsw(n,k)=swlt(ms)-(emimx(2,n,k)-emimx(1,n,k))**2/4.+swlt0(ms)
	wsw(n,k)=swlt(ms)+swlt0(ms)+swltmx(ms)
        wsw(n,k)=wsw(n,k)-(emimx(2,n,k)-emimx(1,n,k))**2/4.
	wsw(n,k)=1./wsw(n,k)
        ENDDO
	enddo
	return
	end

        



