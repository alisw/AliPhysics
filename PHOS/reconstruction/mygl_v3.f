	subroutine mygl
	common /error/ier

	include 'comgl.for'
	include 'comgam.for'
	include 'comalwl.for'
	include 'comwlgen.for'
	include 'comgamx.for'


	call zsavegm
	call zclust
	call cfgln(1)
	
c clasters search
	call wlfn(1)
	call clust(1,0.01)
	call clrwln(0)
c end clusters search

	do ic=1,ncls(1)
c new cluster
                call loadclgl(ic,1)
        	kg=0
        	call wlfn(1)
        	call gamsimpl(1)	
c        	call gamaddtosimpl(1)	
c        	call gamma1(1)	
c                call obmanchik(1,ngggn,ic)
        	call clrwln(0)
c save gamma information
        	call savegm
c goto next cluster if it is
	enddo

	call unsavegm
*	write (*,*) ' Reconstructed gamma = ',kg

	return
	end

        subroutine loadclgl(ic,iw)
	include 'comgl.for'
	include 'comalwl.for'
	include 'comwlgen.for'

	ipwcc=ipwlcl(iw)
	lc=lgclst(ic+ipwcc)
	ipgcc=ipclst(ic+ipwcc)
	ipgl=lc
	nglw(iw)=lc
		do ni=1,lc
		nic=ipgcc+ni
		mgl(ni)=mclst(nic)
		egl(ni)=eclst(nic)
c	write (*,*) ' m,e ',mgl(ni),egl(ni)
		enddo
        return
        end

	subroutine gamsimpl(iw)
        common /error/ier
        common /comkey/keykey
	include 'comgamx.for'
	include 'comalwl.for'
	include 'comgam.for'
	include 'comgl.for'
	include 'comwlgen.for'


	if(iw.le.0.or.iw.gt.nwall) return
	ip=ipwgl(iw)+1
	ng=nglw(iw)

        keykey=0
	call gamma0(ng,mgl(ip),egl(ip),0.07,3.)
	if(ier.ne.0) goto 2
	if(kg.eq.0) goto 2

	call filet
	call filwt
	call diven
	call ecgam
c        write (*,*) ' egam0 ',kg,(e(k),k=1,kg)
	call clrwt

        kgold=kg
        call killnear
        if(kg.eq.kgold) then
                call gamaddtosimpl(1)	
        endif
2	continue
        if(kg.eq.0) then
	call gamma0(ng,mgl(ip),egl(ip),0.035,1.)
	level=-1
c	level=0
	call fitgaml(ng,mgl(ip),egl(ip),hi,himax,mmax,level)
        endif
c        write (*,*) ' fitgaml in ',kg
c        rrxx=sqrt((xwx(1)-xwx(2))**2+(ywx(1)-ywx(2))**2)
c        if(kg.eq.2) then
c        if(ex(1).gt.0.3.and.ex(2).gt.0.3.and.rrxx.lt.65.) then
c        write (*,*) ' kg,kgx ',kg,kgx
c	level=-1
c	level=0
c	level=1
c        keykey=1
c        keykey=0
c	call fitgaml(ng,mgl(ip),egl(ip),hi,himax,mmax,level)
c        keykey=0
c        endif
c        endif
c        write (*,*) ' fitgaml out ',kg

c        keykey=1
c	call filet
c	call filwt
c	call diven
c	call ecgam
c        call simplfite0
c	call clrwt

	return
	end

        subroutine killnear

	include 'comgl.for'
	include 'comgam.for'
	include 'comalwl.for'
	include 'comwlgen.for'

        
        if(kg.le.1) return

        do k1=1,kg-1
        e1=e(k1)
        do k2=k1+1,kg
        e2=e(k2)

        if(e1+e2.gt.0.2) then

        xm=(e1*xw(k1)+e2*xw(k2))/(e1+e2)
        ym=(e1*yw(k1)+e2*yw(k2))/(e1+e2)
        dx1=xw(k1)-xm
        dy1=yw(k1)-ym
        dsc1=dx1*raxay(4,k1)+dy1*raxay(5,k1)
        dx2=xw(k2)-xm
        dy2=yw(k2)-ym
        dsc2=dx2*raxay(4,k1)+dy2*raxay(5,k1)
        dd=sqrt((xw(k1)-xw(k2))**2+(yw(k1)-yw(k2))**2)

        if(dd.lt.30.) then
        if(abs(dsc1).gt.0.1*dd.or.abs(dsc2).gt.0.1*dd) then
                if(e1.gt.e2) then
                       call concatgm(k2,k1)
                else
                       call concatgm(k1,k2)
                endif
                return
        endif
        endif

        endif

        enddo
        enddo

        return
        end

	subroutine gamaddtosimpl(iw)
        common /error/ier
        common /comkey/keykey
	include 'comalwl.for'
	include 'comgam.for'
	include 'comgl.for'
	include 'comwlgen.for'

	if(ifl.eq.0) then

	ifl=1
	endif

	ier=0
	iaddflag=0
	keykey=0
	iadflmax=0

	if(iw.le.0.or.iw.gt.nwall) return
	ip=ipwgl(iw)+1
	ng=nglw(iw)

        keykey=1
c        keykey=0
	call filet
c	call filwt
	call filwltgl(ng,mgl(ip))
	call hiwlnew(ng,mgl(ip),himin,himax,prb,mmax)

	if(kg.gt.0.and.himax.gt.3.85) then
	kgdo=kg
c		call gmplus(ng,mgl(ip),egl(ip),0.07,0.13,0.)
		call gmplus(ng,mgl(ip),egl(ip),0.07,0.13,0.)

c	write (*,*) ' kgdo ,kg ',kgdo,kg

c		if(kg.eq.kgdo) then
cc  hi is bad but gamma not found !		
c        	kgold=kg
c        	call gmpluss(ng,mgl(ip),egl(ip),0.085,0.13,40.)
c	write (*,*) ' vtoroi kgold ,kg ',kgold,kg
c                endif

		if(kg.gt.kgdo) then
c       	level=0
        	level=-1
        	call fitgaml(ng,mgl(ip),egl(ip),hi,himax,mmax,level)
                endif
	endif	



c	call clrwt
	call clrwltgl(ng,mgl(ip))
        keykey=0
        return
        end
        

