        SUBROUTINE CLUST(iwl,gamthr)
C
C      9.12.92  Kolosov V.N.
C        GLASS CLUSTERS SEARCH
C
        common /error/ier

	include 'comgl.for'
	include 'comalwl.for'
	include 'comwlgen.for'

	integer mark(nglal),icl(ialw)


	miwl=1
	mawl=nwall
c	write (*,*) ' nwall ',nwall
	if(iwl.gt.0.and.iwl.le.nwall) then
		miwl=iwl
		mawl=iwl
	endif

	do iw=miwl,mawl
	ip=ipwgl(iw)
	ng=nglw(iw)
c	call wlf(wl,ng,mgl(ip+1),egl(ip+1))
	ncls(iw)=0
	ipwlcl(iw)=ntclst

	nc=ntclst
	do k=ip+1,ip+ng
	m=mgl(k)
	e=egl(k)
c	write (*,*) ' m,e ',m,e
	if(icl(m).eq.0.and.e.gt.gamthr) then
c New cluster
		if(nc.eq.nmxcl) then
		write (*,*) ' max cluster err !!! ',nmxcl
		ier=36
		return
		endif
		nc=nc+1
co	write (*,*) ' new cluster ',nc,ipcl
		icl(m)=nc
		ipclst(nc)=ipcl
		ncl=1
		ncld=1
		DO WHILE (NCLD.gt.0)
		NCLD=0
		DO 2 K1=ip+1,ip+NG
co	write (*,*) ' neib? #,m,icl',k1-ip,mgl(k1),icl(mgl(k1))
		IF(MARK(K1).NE.0) GOTO 2
		M1=MGl(K1)
		IF(ICL(M1).NE.NC) GOTO 2
co	write (*,*) ' search ns,wl',ns(m1),wl(m1)
			DO L=1,NS(M1)
			M1S=KS(L,M1)
co	write (*,*) ' sosed ',l,m1s,wl(m1s)
			IF(WL(M1S).NE.0..AND.ICL(M1S).EQ.0) THEN
				ICL(M1S)=NC
				NCLD=NCLD+1
			ENDIF
			ENDDO
		MARK(K1)=NC
c Save here
		if(ipcl.eq.nglal) then
		write (*,*) ' over cluster buff !!! ',ipcl
		ier=37
		return
		endif
		eclst(ipcl+1)=egl(k1)
		mclst(ipcl+1)=mgl(k1)
		ipcl=ipcl+1
2		CONTINUE
		ncl=ncl+ncld
		ENDDO
		lgclst(nc)=ncl
C   END  NEW  CLUSTER		
	endif
	enddo
	ncls(iw)=nc-ntclst
	ntclst=ntclst+nc
c  zero marks
	DO K=ip+1,ip+NG
	m=mgl(k)
	MARK(K)=0
	icl(m)=0	
	ENDDO
c New wall
	enddo

	return

	entry zclust

	ipcl=0
	ntclst=0
	do i=1,nalw
	ncls(i)=0
	enddo
	do i=1,ialw
	icl(i)=0
	enddo	
	do i=1,nglal
	mark(i)=0
	enddo	

	return
	end

        SUBROUTINE cCLUST(iwl)
        common /error/ier

	include 'comgl.for'
	include 'comalwl.for'
	include 'comwlgen.for'

	real xt(3),xmi(3),xma(3)

	miwl=1
	mawl=nwall
	if(iwl.gt.0.and.iwl.le.nwall) then
		miwl=iwl
		mawl=iwl
	endif
	do iw=miwl,mawl
c new wall
	ipwcc=ipwlcl(iw)
	etwc=0.
	nck(iw)=ncls(iw)

	do i=1,ncls(iw)
c new cluster
	lc=lgclst(i+ipwcc)
	ipgcc=ipclst(i+ipwcc)
	et=0.
	emax=0.
	call vzero(xt,3)
	call ucopy(x(1,mclst(ipgcc+1)),xmi,3)
	call ucopy(x(1,mclst(ipgcc+1)),xma,3)
		do ni=1,lc
c new glass in cluster
		nic=ipgcc+ni
		mgc=mclst(nic)
		egc=eclst(nic)
		if(emax.lt.egc) then
			emax=egc
			mmax=mgc
		endif
		if(xmi(1).gt.x(1,mgc)) xmi(1)=x(1,mgc) 
		if(xmi(2).gt.x(2,mgc)) xmi(2)=x(2,mgc) 
		if(xmi(3).gt.x(3,mgc)) xmi(3)=x(3,mgc) 
		if(xma(1).lt.x(1,mgc)) xma(1)=x(1,mgc) 
		if(xma(2).lt.x(2,mgc)) xma(2)=x(2,mgc) 
		if(xma(3).lt.x(3,mgc)) xma(3)=x(3,mgc) 

		xt(1)=xt(1)+x(1,mgc)*egc
		xt(2)=xt(2)+x(2,mgc)*egc
		xt(3)=xt(3)+x(3,mgc)*egc
		et=et+egc
c end glass in cluster
		enddo
	if(et.eq.0.) then
		write (*,*) ' error!! energy cluster=0 '
		write (*,*) ' wall,nc,lc,ipgcc ',iw,i,lc,ipgcc
		ier=33
		return
	endif
	xt(1)=xt(1)/et
	xt(2)=xt(2)/et
	xt(3)=xt(3)/et
	
	eck(i,iw)=et
	emxck(i,iw)=emax
	call ucopy(xt,xavck(1,i,iw),3)
	call ucopy(xmi,xmick(1,i,iw),3)
	call ucopy(xma,xmack(1,i,iw),3)

	etwc=etwc+et
c end cluster
	enddo
	etck(iw)=etwc
c end wall
	enddo

	return
	end
