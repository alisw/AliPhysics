        SUBROUTINE GAMMA0(ng,mg,eg,glthr,ff)

C
C K+  12.90  KOLOSOV V.N.
C     14.03.92
C              SIMPLE GAMMA SEARCH
C
	dimension mg(*),eg(*)
        common /error/ier
	include 'comwlgen.for'
	include 'comalwl.for'
	include 'comgam.for'

	dimension e3na3(3,3),ijdd(2)

	data kgmax /1000/
	
        IF(NG.EQ.0) RETURN
        DO 1 N=1,NG
        EM=EG(N)
	
        M=MG(N)
        IF(EM.LT.GLTHR) GOTO 1
c        IF(EM.LT.GmTHR(m)) GOTO 1
	if(iflgmthr.eq.0) then
	iflgmthr=1
c	write (*,*) ' Warning !! no glass thresholds now !'
	endif

c piks search
	call findpik(m,iokpik,e3na3,ff)



	if(iokpik.eq.0) goto 1

	ijdd(1)=1
	ijdd(2)=1
	if(e3na3(1,2).gt.e3na3(3,2)) ijdd(1)=-1
	if(e3na3(2,1).gt.e3na3(2,3)) ijdd(2)=-1
	
c new gamma
          IF(KG.GT.kgmax) THEN
            IER=32
            RETURN
          ENDIF

        xx=X(1,M)
        yy=x(2,M)
        ee=WL(M)
	see=ee/4.
	sxx=d(1,m)/2.
	syy=d(2,m)/2.
	call adgmnew(m,ee,xx,yy,see,sxx,syy,ijdd)

1       CONTINUE

        RETURN
        END

	subroutine savegm

	include 'comgam.for'
        COMMON /SVGAMMA/KGsv,MWsv(ngp),IDsv(ngp),JDsv(ngp),svE(ngp),
     ,  svE4(ngp),
     ,  svXW(ngp),svYW(ngp),svES(nsps,ngp),svET(nsps,ngp),
     ,  ISsdsv(ngp),
     ,  IGDEVsv(ngp),svZGDEV(ngp),svsigexy(3,ngp),
     ,  svEmimx(2,nsps,ngp),
     ,  kgfixsv,svigfix(ngp),svcgfix(3,ngp),svsgfix(3,ngp),
     ,  svhiw(ngp),
     ,  svwsw(nsps,ngp),svh1w(ngp),svh0w(ngp)

	if(kg.eq.0) return

		do k=1,kg
		mwsv(kgsv+k)=mw(k)
		idsv(kgsv+k)=id(k)
		jdsv(kgsv+k)=jd(k)
		sve(kgsv+k)=e(k)
		sve4(kgsv+k)=e4(k)
		svxw(kgsv+k)=xw(k)
		svyw(kgsv+k)=yw(k)

		do kn=1,nsps
			sves(kn,kgsv+k)=es(kn,k)
			svet(kn,kgsv+k)=et(kn,k)
			svemimx(1,kn,kgsv+k)=emimx(1,kn,k)
			svemimx(2,kn,kgsv+k)=emimx(2,kn,k)
			svwsw(kn,kgsv+k)=wsw(kn,k)
		enddo
		issdsv(kgsv+k)=issd(k)
		igdevsv(kgsv+k)=igdev(k)
		svzgdev(kgsv+k)=zgdev(k)
		call ucopy(sigexy(1,k),svsigexy(1,kgsv+k),3)
		enddo
		kgsv=kgsv+kg
	return
	
	entry zsavegm
		kgsv=0
	return

	entry unsavegm
		kg=kgsv
		do k=1,kgsv
		mw(k)=mwsv(k)
		id(k)=idsv(k)
		jd(k)=jdsv(k)
		e(k)=sve(k)
		e4(k)=sve4(k)
		xw(k)=svxw(k)
		yw(k)=svyw(k)

		do kn=1,nsps
			es(kn,k)=sves(kn,k)
			et(kn,k)=svet(kn,k)
			emimx(1,kn,k)=svemimx(1,kn,k)
			emimx(2,kn,k)=svemimx(2,kn,k)
			wsw(kn,k)=svwsw(kn,k)
		enddo
		issd(k)=issdsv(k)
		igdev(k)=igdevsv(k)
		zgdev(k)=svzgdev(k)
		call ucopy(sigexy(1,k),svsigexy(1,k),3)
		enddo
	return
	end
        SUBROUTINE findpik(m,iokpik,e3na3,ffls)

        common /error/ier
	include 'comwlgen.for'
	include 'comalwl.for'
	include 'comgam.for'

	real e3na3(3,3)

	
	
	ff=ffls
	if(ffls.le.0.) ff=1. 
        em=wl(m)

c ALL arround less then em
c	do k=2,ns(m)
c	ms=ks(k,m)
c	if(em.lt.wl(ms)) goto 1
c	enddo
	
	call vzero(e3na3,9)
	
        IF(NKS(M).EQ.0) THEN
c simple case
          ILS=KS(2,M)
          IF(EM.LE.WL(ILS)) GOTO 1
          EL=WL(ILS)
          IRS=KS(3,M)
          IF(EM.LT.WL(IRS)) GOTO 1
          ER=WL(IRS)
          IUS=KS(4,M)
          IF(EM.LT.WL(IUS)) GOTO 1
          EU=WL(IUS)
          IDS=KS(5,M)
          IF(EM.LE.WL(IDS)) GOTO 1
          ED=WL(IDS)
c L-sosedi
          ILUS=KS(6,M)
          IF(EM.LE.WL(ILS)/ff) GOTO 1
          ELU=WL(ILUS)
          IRUS=KS(7,M)
          IF(EM.LT.WL(IRUS)/ff) GOTO 1
          ERU=WL(IRUS)
          ILDS=KS(8,M)
          IF(EM.LT.WL(ILDS)/ff) GOTO 1
          ELD=WL(ILDS)
          IRDS=KS(9,M)
          IF(EM.LE.WL(IRDS)/ff) GOTO 1
          ERD=WL(IRDS)


c new gamma in standart case
         ELSE
c complex case
          NLS=IBITS(NKS(M),0,2)
          NRS=IBITS(NKS(M),2,2)
          NUS=IBITS(NKS(M),4,2)
          NDS=IBITS(NKS(M),6,2)
          NLUS=IBITS(NKS(M),8,1)
          NRUS=IBITS(NKS(M),9,1)
          NLDS=IBITS(NKS(M),10,1)
          NRDS=IBITS(NKS(M),11,1)
          L=1
          EL=0.
          DO K=1,NLS
          L=L+1
          ILS=KS(L,M)
          IF(EM.LE.WL(ILS)) GOTO 1
          EL=EL+WL(ILS)
          ENDDO
          ER=0.
          DO K=1,NRS
          L=L+1
          IRS=KS(L,M)
          IF(EM.LT.WL(IRS)) GOTO 1
          ER=ER+WL(IRS)
          ENDDO
          EU=0.
          DO K=1,NUS
          L=L+1
          IUS=KS(L,M)
          IF(EM.LE.WL(IUS)) GOTO 1
          EU=EU+WL(IUS)
          ENDDO
          ED=0.
          DO K=1,NDS
          L=L+1
          IDS=KS(L,M)
          IF(EM.LE.WL(IDS)) GOTO 1
          ED=ED+WL(IDS)
          ENDDO
c L-sosedi
          ELU=0.
          DO K=1,NLUS
          L=L+1
          ILUS=KS(L,M)
          IF(EM.LE.WL(ILUS)/ff) GOTO 1
          ELU=ELU+WL(ILUS)
          ENDDO
          ERU=0.
          DO K=1,NRUS
          L=L+1
          IRUS=KS(L,M)
          IF(EM.LT.WL(IRUS)/ff) GOTO 1
          ERU=ERU+WL(IRUS)
          ENDDO
          ELD=0.
          DO K=1,NLDS
          L=L+1
          ILDS=KS(L,M)
          IF(EM.LE.WL(ILDS)/ff) GOTO 1
          ELD=ELD+WL(ILDS)
          ENDDO
          ERD=0.
          DO K=1,NRDS
          L=L+1
          IRDS=KS(L,M)
          IF(EM.LE.WL(IRDS)/ff) GOTO 1
          ERD=ERD+WL(IRDS)
          ENDDO
c end of both cases
        ENDIF

	e3na3(2,2)=em
	e3na3(1,2)=el
	e3na3(3,2)=er
	e3na3(2,1)=ed
	e3na3(2,3)=eu
	e3na3(1,3)=elu
	e3na3(3,3)=eru
	e3na3(1,1)=eld
	e3na3(3,1)=erd

	iokpik=1
	return

1       iokpik=0
        RETURN
        END
c
	subroutine adgmnew(m,ee,xx,yy,see,sxx,syy,ijdd)
        common /error/ier
	include 'comwlgen.for'
	include 'comalwl.for'
	include 'comgam.for'

	dimension ijdd(2)

	kg=kg+1
	mw(kg)=m
	id(kg)=ijdd(1)
	jd(kg)=ijdd(2)
c set wall marker to gamma
	igdev(kg)=0
	do iw=1,nwall
	if(m.gt.madwl(iw).and.m.le.madwl(iw+1)) then
		igdev(kg)=iw
		goto 555
	endif
	enddo
555	continue
c gamma coord och. grubo
        XW(KG)=xx
        YW(KG)=yy
        E(KG)=ee
        zgdev(kg)=zmidlshower(ee,m)
	sigexy(1,kg)=see
	sigexy(2,kg)=sxx
	sigexy(3,kg)=syy
        DO K=1,NS(M)
        LS=KS(K,M)
        ES(K,KG)=WL(LS)
        ENDDO
	call filraxay(kg)
c end new gamma
	return
	end
C
	subroutine filraxay(k)
	include 'comalwl.for'
	include 'comgam.for'
	include 'comwlgen.for'

        data ifl/0/
        if(ifl.eq.0) then
c        write (*,*) ' FILRAXAY WARNING dz=4600. eto SMESHNO !!! '
        ifl=1
        endif
        
	dz=4600.
c	dz=xyzwall(3,iw)-xyzvtx(3)

	m=mw(k)

	raxay(1,k)=xw(k)
	raxay(2,k)=yw(k)
	raxay(3,k)=zgdev(k)
c	raxay(3,k)=46.
	raxay(4,k)=-x(1,m)/dz
	raxay(5,k)=-x(2,m)/dz
c	raxay(4,k)=xw(k)/dz
c	raxay(5,k)=yw(k)/dz

	return
	end

        SUBROUTINE FILWT0
	include 'comalwl.for'
	include 'comgam.for'

        common /comdevpar/sigmaph,sigmapd,sigphsq,sigpdsq
	data ifl/0/
c	write (*,*) ' filwt0 '

	if(ifl.eq.0) then
		do n=1,nt
		swlt0(n)=sigpdsq
		enddo
	write (*,*) ' filwt0 OK  !!! '
	ifl=1
        else
	 write (*,*) ' filwt0 ne nado !!! '
        
	endif

        return
        end
        
        SUBROUTINE FILWT
	include 'comalwl.for'
	include 'comgam.for'

        common /comdevpar/sigmaph,sigmapd,sigphsq,sigpdsq
	data ifl/0/
C

c	write (*,*) ' filwt '

	if(ifl.eq.0) then
	sigmaped=0.01

c        call filwt0
	ifl=1
	endif


        IF(KG.EQ.0) RETURN
        DO K=1,KG
        M=MW(K)

        DO N=1,NS(M)
        MS=KS(N,M)
        WLT(MS)=WLT(MS)+ET(N,K)
c	swlt(ms)=swlt(ms)+(emimx(2,n,k)-emimx(1,n,k))**2/4.
	swlt(ms)=swlt(ms)+dispeces(n,k)
        ENDDO
        ENDDO

        RETURN
        END
C
        SUBROUTINE CLRWT
	include 'comalwl.for'
	include 'comgam.for'
C
        IF(KG.EQ.0) RETURN
        DO K=1,KG
        M=MW(K)

        DO N=1,NS(M)
        MS=KS(N,M)
        WLT(MS)=0.
        SWLT(MS)=0.
        swltmx(ms)=0.
        ENDDO
        ENDDO

	call vzero(wlt,nt)

        RETURN
        END
C

        SUBROUTINE DIVEN

	include 'comalwl.for'
	include 'comgam.for'

        IF(KG.LE.1) RETURN
        DO K=1,KG
	h1w(k)=0.
	h0w(k)=0.
        MK=MW(K)
          DO N=1,NS(MK)
          MN=KS(N,MK)
	          ES(N,K)=WL(MN)*ET(N,K)/WLT(MN)
          ENDDO

        ENDDO

        RETURN
        END

        subroutine simplfite0
	include 'comalwl.for'
	include 'comgam.for'
        do k=1,kg
        mk=mw(k)
        sum1=0.
        sum2=0.
        do n=1,ns(mk)
        ms=ks(n,mk)
        en=et(n,k)/e(k)
        disp=swlt(ms)
        sum1=sum1+es(n,k)*en/disp
        sum2=sum2+en*en/disp
        enddo
        enew=sum1/sum2
c        write (*,*) ' k, e_old, e_new ',k,e(k),enew
        enddo
        return
        end


        SUBROUTINE ECGAM
C
C K+  12.90  KOLOSOV V.N.
C     14.03.92
C             E,X,Y   GAMMA CALCULATION
C
	include 'comalwl.for'
	include 'comgam.for'
        COMMON /COORF/CFUN(200)
        DIMENSION EEE(9),NNN(9),KEY(3,3)
        DATA KEY  /   8,5,9,
     ,                2,1,3,
     ,                6,4,7/
C
        IF(KG.EQ.0) RETURN
        DO 1 N=1,KG
        M=MW(N)
        IF(NKS(M).EQ.0) THEN
          DO I=1,9
          EEE(I)=ES(I,N)
          ENDDO
         ELSE
          NNN(1)=1
          NNN(2)=IBITS(NKS(M),0,2)
          NNN(3)=IBITS(NKS(M),2,2)
          NNN(4)=IBITS(NKS(M),4,2)
          NNN(5)=IBITS(NKS(M),6,2)
          NNN(6)=IBITS(NKS(M),8,1)
          NNN(7)=IBITS(NKS(M),9,1)
          NNN(8)=IBITS(NKS(M),10,1)
          NNN(9)=IBITS(NKS(M),11,1)
          L=1
          DO I=1,9
          EEE(I)=0.
          DO K=1,NNN(I)
          EEE(I)=EEE(I)+ES(L,N)
          L=L+1
          ENDDO
          ENDDO
        ENDIF
C
        E(N)=0.
        DO K=1,9
        E(N)=E(N)+EEE(K)
        ENDDO
        E4(N)=0.
C
        DO I=2,2+ID(N),ID(N)
        DO J=2,2+JD(N),JD(N)
        K=KEY(I,J)
        E4(N)=E4(N)+EEE(K)
        ENDDO
        ENDDO

        IF(ID(N).EQ.1) THEN
          IL=2
         ELSE
          IL=1
        ENDIF

        IF(JD(N).EQ.1) THEN
          JL=2
         ELSE
          JL=1
        ENDIF

        EL=0.
        DO J=2,2+JD(N),JD(N)
        K=KEY(IL,J)
        EL=EL+EEE(K)
        ENDDO

        ED=0.
        DO I=2,2+ID(N),ID(N)
        K=KEY(I,JL)
        ED=ED+EEE(K)
        ENDDO
	if(el.eq.0.and.ed.eq.0) then
		write (*,*) ' ECGAM giagnostic EG=0 '
		write(*,*) ' KG,ng,id,jd ',kg,n,id(n),jd(n)
		write(*,*) ' egam',(es(i,n),i=1,9)
	endif
        DHx=D(1,M)/2.
        DHy=D(2,M)/2.

        X0=X(1,M)+ID(N)*DHx
        Y0=x(2,M)+JD(N)*DHy

        ax=raxay(4,n)
        ay=raxay(5,n)

c        write (*,*) ' ax,ay ',ax,ay
	xx=coornew(el,e4(n)-el,ax,m)
	yy=coornew(ed,e4(n)-ed,ay,m)

c        write (*,*) ' xx,yy ',xx,yy
c put cell center if too far
        IF(XX.LT.-DHx) XX=-DHx
        IF(XX.GT. DHx) XX= DHx
        IF(YY.LT.-DHy) YY=-DHy
        IF(YY.GT. DHy) YY= DHy

        
        XW(N)=X0+XX
        YW(N)=Y0+YY

        dxm=xw(n)-x(1,mw(n))
        if(ax.le.0.) then
        xxc=dxm
        corr=1.076+0.163*xxc-0.00955*xxc**2-0.0013215*xxc**3
        else
        xxc=-dxm
        corr=1.076+0.163*xxc-0.00955*xxc**2-0.0013215*xxc**3
        corr=-corr
        endif
        xw(n)=xw(n)-corr
        
        dxm=yw(n)-x(2,mw(n))
        if(ay.le.0.) then
        xxc=dxm
        corr=1.076+0.163*xxc-0.00955*xxc**2-0.0013215*xxc**3
        else
        xxc=-dxm
        corr=1.076+0.163*xxc-0.00955*xxc**2-0.0013215*xxc**3
        corr=-corr
        endif
        yw(n)=yw(n)-corr
        
 
	call filraxay(n)
1       CONTINUE

c        call simplfite0
        call betterthanfit


        RETURN
        END

        subroutine betterthanfit
	include 'comalwl.for'
	include 'comgam.for'
	include 'comgl.for'
	include 'comwlgen.for'
c ADD energy out 3*3 cells for crystalls only
	do k=1,kg
	e(k)=1.065*e(k)
        dxm=xw(k)-x(1,mw(k))
        dym=yw(k)-x(2,mw(k))
        rm=sqrt(dxm**2+dym**2)
        e(k)=e(k)*(1.+0.01*rm/15.)
        ee=e(k)
        e(k)=ee*(1.+0.03*exp(-ee))
	enddo
        return
        end

	subroutine filet
	include 'comalwl.for'
	include 'comgam.for'
        common /comdevpar/sigmaph,sigmapd,sigphsq,sigpdsq
	real rxxyy(5)

	if(kg.eq.0) return

	do k=1,kg
	m=mw(k)

	xg=xw(k)
	yg=yw(k)

	do n=1,ns(m)
	ms=ks(n,m)
	xs=x(1,ms)
	ys=x(2,ms)

	et(n,k)=e(k)*ampcelnew(e(k),raxay(1,k),x(1,ms),d(1,ms),ms)

	if(et(n,k).le.0) then

		write (*,*) ' et=0!! kg,k,m,n,e(k) ',kg,k,m,n,e(k)
		et(n,k)=e(k)*0.0001
	endif

	emimx(1,n,k)=et(n,k)
	emimx(2,n,k)=et(n,k)

        dispeces(n,k)=
     =  (sigwlgam0(et(n,k),e(k)))**2+sigmphsq*et(n,k)

	sigmaes0(n,k)=sqrt(dispeces(n,k))

	enddo

	enddo

        RETURN
        END

	subroutine killgm(kkg)

        common /error/ier
	include 'comgam.for'
	include 'comalwl.for'

	M=MW(Kkg)


	if(kkg.eq.kg) then
		kg=kg-1
	else 
		kg=kg-1
		do k=kkg,kg
		mw(k)=mw(k+1)
		id(k)=id(k+1)
		jd(k)=jd(k+1)
		e(k)=e(k+1)
		e4(k)=e4(k+1)
		xw(k)=xw(k+1)
		yw(k)=yw(k+1)

		do kn=1,nsps
			es(kn,k)=es(kn,k+1)
			et(kn,k)=et(kn,k+1)
			emimx(1,kn,k)=emimx(1,kn,k+1)
			emimx(2,kn,k)=emimx(2,kn,k+1)
			wsw(kn,k)=wsw(kn,k+1)
		enddo
		issd(k)=issd(k+1)
		igdev(k)=igdev(k+1)
		zgdev(k)=zgdev(k+1)
		call ucopy(sigexy(1,k+1),sigexy(1,k),3)
		call ucopy(raxay(1,k+1),raxay(1,k),5)
		enddo
	endif
	return
	end

	subroutine concatgm(k1,k2)

        common /error/ier
	include 'comgam.for'
	include 'comalwl.for'

        if(k1.gt.kg.or.k1.le.0) return
        if(k2.gt.kg.or.k2.le.0) return
        enew=e(k1)+e(k2)
        xnew=(xw(k1)*e(k1)+xw(k2)*e(k2))/enew
        ynew=(yw(k1)*e(k1)+yw(k2)*e(k2))/enew
        e(k2)=enew
        xw(k2)=xnew
        yw(k2)=ynew
        raxay(1,k2)=xnew
        raxay(2,k2)=ynew
        id(k2)=1
        jd(k2)=1
        if(xnew.lt.x(1,mw(k2))) id(k2)=-1
        if(ynew.lt.x(2,mw(k2))) jd(k2)=-1

        call killgm(k1)

	call filet
	call filwt
	call diven
	call ecgam
	call clrwt

	return
	end



