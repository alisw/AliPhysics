	function acl3(e,vx,vc,vb)

	real vx(5),vc(3),vb(3)
	common /comsav/nzsav,vxsav(5),esav,ezsav(100)

	if(e.le.0.) then 
c		write (*,*) ' acl3 ',e
		acl3=0.
		return
	endif
	acl3=0.
c esli cell daleko to nechego schitat'!	
c
cold	ax=vx(4)
cold	ay=vx(5)
cnew?
	ax=-vx(4)
	ay=-vx(5)
	dcx=vb(1)
	dcy=vb(2)
	dcz=vb(3)
c	write (*,*) ' ax,ay,vx,vy,vz',ax,ay,dcx,dcy,dcz
c	write (*,*) ' x,y,z ',vx(1),vx(2),vx(3)
c	write (*,*) ' xc,yc,zc ',vc(1),vc(2),vc(3)

c	nz=10
	nz=25
c	nz=3
	
	dz=200./nz
	nz=vb(3)/dz
	dz=vb(3)/nz
	dl=dz

c?????	z0=vc(3)-dcz/2.
	z0=0.
	xz0=vx(1)-vc(1)+ax*(z0-vx(3))
	yz0=vx(2)-vc(2)+ay*(z0-vx(3))

	adx=abs(xz0)
	ady=abs(yz0)
	iflpr=0
	if(adx.gt.80..or.ady.gt.80.) then
		acl3=0.0001
		return
c		write (*,*) ' adx,ady ',adx,ady
c		iflpr=1
	endif


	if(e.ne.esav.or.nzsav.ne.nz) then
		esav=e
		nzsav=nz
		call ucopy(vx(1),vxsav(1),5)
		do iz=1,nz
		z=z0+dz*(float(iz)-0.5)
		cl=z
		ezsav(iz)=azlngd(e,cl,dz)
		enddo
	endif


	do iz=1,nz
	z=z0+dz*(float(iz)-0.5)
	cl=z
c	ez=azlngd(e,cl,dz)
	ez=ezsav(iz)
c xz,yz gamma coord from cell center at z level
	xz=xz0+ax*z
	yz=yz0+ay*z
	eez=aclz(xz,yz,dcx,dcy,z)
	acl3=acl3+ez*eez
c	write (*,*) ' iz,ez,eez,cl ',iz,ez,eez,cl
	enddo

c	if(iflpr.eq.1) then
c		iflpr=0
c		write(*,*) ' acl3= ',acl3
c	endif
c		write(*,*) ' acl3= ',acl3
	return
	end
c along z
	function aznorm(e,dz)
c	--------------------
	data zmax/500./
	en=0.
	do iz=1,10000 
	z=dz*(float(iz)-0.5)
	if(z.gt.zmax) goto 1
	ead=azlng(e,z)
	en=en+ead
	enddo
1	aznorm=en
	return
	end

	function azlngd(e,z,dz)
	et=aznorm(e,dz)
	azlngd=azlng(e,z)/et
	return
	end

c longditude distribution 
	function azlng(e,z)
c	--------------------
	common/shdat/ rl,zmxr,er,tmx
c	data rl,er/11.0,0.025/
c	data rl,er/17.6,0.025/
	data rl,er/9.2,0.025/
	if(z.lt.0.) then
		azlng=0.
		return
	endif
	if(e.le.0) then
		write (*,*) ' e,z ',e,z
		azlng=0
		return
	endif
c for electron ? 
c	tmx=alog(e/er)-0.5
c for gamma?
	tmx=alog(e/er)+0.5
	tm=tmx
c directly
	azlng=azlnt(tm,z)
c by tables
c	azlng=azltt(tm,z)
	return
	end
c ---------------------------------------------------------
	function azlnt(tm,z)
c	--------------------
	common/shdat/ rl,zmxr,er,tmx

	t=z/rl
c	write (*,*) ' tmx ',tmx
	x=t/2.
	xm=tm/2.
c	if(x.gt.10.) write (*,*) ' x ',x
	if(x.gt.40.) then
		azlnt=0.
	else
		azlnt=exp(xm*alog(x))*exp(-x)
	endif
	return
	end
c =========
c And the same by tables   ******************************
	function azltt(tmx,z)
	common /comaltt/nx,ny,xmi,xma,ymi,yma,att(100,200)

	if(ifl.eq.0) then
	open(20,file='att.ffun',form='formatted',status='unknown')
	read (20,*) nx,ny,xmi,xma,ymi,yma
	read (20,*) att
	close (20)
	ifl=1
	endif
	azltt=f2tab(att,nx,tmx,xmi,xma,nx,z,ymi,yma,ny)
	return
	end
c ---------------------------------------------------------
c   2-dim functions
        FUNCTION ACLz(XG,YG,dax,day,z)
c      -------------------------------      
C       XG,YG GAMMA COOR. X=0,Y=0 - CENTER OF THE CELL
C       dx,dy - SIZE OF THE CELL (MM)
C
        EXTERNAL A2Fz
        X=-ABS(XG)
        Y=-ABS(YG)
        Dx=DAx/2.
        Dy=DAy/2.
        XP=X+Dx
        YP=Y+Dy
        XM=X-Dx
        YM=Y-Dy
        ACLz=A2Fz(XP,YP,z)+A2Fz(XM,YM,z)-A2Fz(XP,YM,z)-A2Fz(XM,YP,z)
        RETURN
        END

c =============================================================

	function a2fz(xc,yc,z)
	a2fz=a2fzo(xc,yc,z)
	return
	end

c =============================================================
        FUNCTION a2FzO(xc,yc,z)
	ax=afz(xc,z)
	ay=afz(yc,z)
c directly
	a2fzo=aa2fz(ax,ay)
c by tables
c	a2fzo=aa2fzt(ax,ay)
        RETURN
	END
c ------------------------------------------------------------
	function aa2fz(ax,ay)
C        DATA  NPAR,C0,C1,C2/3,0.34,0.00,11.00/
C  LEDNEV'S FUNCTION
        DATA  NPAR,C0,C1,C2/3,0.82,-2.2,24.30/
        DATA  CPI/6.283186/

	X=0.5-AX
	Y=0.5-AY
	X2=2.*X
	Y2=2.*Y
	PNX=X*(1.-X2*X2)
	PNY=Y*(1.-Y2*Y2)
        PXY=C0+C1*(X**2+Y**2)+C2*X**2*Y**2
c        ESIN=SIN(CPI*AX)*SIN(CPI*AY)*PXY
	aa2fz=AX*AY+PNX*PNY*PXY
	return
	end

c And the same by tables  *********************
        function aa2fzt(ax,ay)
	common /comaa2/nx,ny,xmi,xma,ymi,yma,aa2(100,100)
	if(ifl.eq.0) then
	open(20,file='aa2.ffun',form='formatted',status='unknown')
	read (20,*) nx,ny,xmi,xma,ymi,yma
	read (20,*) aa2
	close (20)

	ifl=1
	endif
	aa2fzt=f2tab(aa2,nx,ax,xmi,xma,nx,ay,ymi,yma,ny)
	return
	end
c ==============================================
c 1-dim functions
        FUNCTION AFz(X,z)
C       -------------
C    SHOWER ENERGY FRACTION RIGHT FROM BOUNDARY
C    BOUNDARY COOR X=0.
	if(x.gt.0) then
		afz=1.-cumprf(x,z)
	else
		afz=cumprf(-x,z)
	endif
	if(afz.lt.0) then
		write (*,*) ' afz ',afz
	endif
	return
	end
c -----
	function cumprf(cx,cl)
	common/shdat/ rl,zmxr,er,tmx

c        write (*,*) 'rl,zmxr,er,tmx',rl,zmxr,er,tmx
	clmx=rl*(tmx+2.)

	dc=clmx-cl

c	dcarg=0.9*dc/clmx
c	dcarg=3.5*dc/clmx
	dcarg=0.9*dc/clmx
c	write (*,*) ' dcarg ',dcarg
	if(dcarg.lt.-40.) then
		cfcl=0.
	elseif(dcarg.lt.0.) then
c		cfcl=0.71*1.15*exp(dcarg)	
		cfcl=0.71*1.15*exp(dcarg)	
	else
c		cfcl=0.71*1.15*exp(0.5*dcarg**2)	
		cfcl=0.71*1.15*exp(dcarg)	
	endif
c	write (*,*) ' dcarg ',dcarg,cfcl
	cumprf=cumcryex(cx*cfcl)
c	if(cumprf.lt.0.) write (*,*) ' cumprf ',cumprf
	if(cumprf.lt.0.) cumprf=0.
	if(cumprf.gt.1.) write (*,*) ' cumprf ',cumprf
	return
	end
c ==================================================================
	function cumcryex(xc)
	x=abs(xc)
c directly
	if(x.gt.80.) then
	a=0.
	else
	a=acumex(x)
	endif
c by tables
c	a=acumext(x)
	if(xc.ge.0.) then
		cumcryex=a
	else
		cumcryex=1.-a
	endif
	return
	end
c ==================================================================
        function acumex(x)
c	data p1,p2 / 2.526, 1.104/
	data p1,p2 /1.480, 1.600/
	a=0.5*exp((p1-sqrt(p1**2+4.*p2*x))/2./p2)
	acumex=a
	return
	end
c And the same by tables   ****************************
        function acumext(x)
	common /comext/nx,xmi,xma,acux(2000)
	if(ifl.eq.0) then
	open(20,file='ac1.ffun',form='formatted',status='unknown')
	read (20,*) nx,xmi,xma
	read (20,*) acux
	close (20)

	ifl=1
	endif
	acumext=f1tab(acux,x,xmi,xma,nx)
	return
	end
c ==================================================================
C             SHOWER FUNCTIONES DERIVATIVES
C
C                                               23.11.90  KOLOSOV V.N.
C                                                  ??        OTLADKA
c                                               18.03.94  continue buging
c                                               july.96 
	function dacl3(nd,e,vx,vc,vb)
	real vx(5),vc(3),vb(3)
	real vxp(5),vcp(3),vbp(3)
	dacl3=0.


	if(e.le.0.) then 
c		write (*,*) ' dacl3 energy',e
		dacl3=0.
		return
	endif
	if(nd.lt.1.or.nd.gt.4) return
	if(nd.eq.1) then
		dacl3=dxacl3(e,vx,vc,vb)
	elseif(nd.eq.2) then
		call ucopy(vx,vxp,5)
		call ucopy(vc,vcp,3)
		call ucopy(vb,vbp,3)
		call chng(vxp,1,2)
		call chng(vxp,4,5)
		call chng(vcp,1,2)
		call chng(vbp,1,2)
		dacl3=dxacl3(e,vxp,vcp,vbp)
	endif
	return
	end
c ---------------------------------------------------------
	function dxacl3(e,vx,vc,vb)
	common /comsav/nzsav,vxsav(5),esav,ezsav(100)

	common /acflall/ precx,precy
	real vx(5),vc(3),vb(3)
	data afin/0.01/
	dxacl3=0.
	ax=vx(4)
	ay=vx(5)
	dcx=vb(1)
	dcy=vb(2)
	dcz=vb(3)

	nz=10
c	nz=3



c	do iz=1,nz
c	z=z0+dz*(float(iz)-0.5)
c	cl=z
c	ez=azlngd(e,cl,dz)
cc xz,yz gamma coord from cell center at z level
c	xz=xz0+ax*z
c	yz=yz0+ay*z
c	dddd=dxaclz(xz,yz,dcx,dcy,z)
c	dxacl3=dxacl3+ez*dddd
c	enddo


	dz=200./nz
	nz=vb(3)/dz
	dz=vb(3)/nz
	dl=dz

c?????	z0=vc(3)-dcz/2.
	z0=0.
	xz0=vx(1)-vc(1)+ax*(z0-vx(3))
	yz0=vx(2)-vc(2)+ay*(z0-vx(3))

	adx=abs(xz0)
	ady=abs(yz0)
	iflpr=0
	if(adx.gt.80..or.ady.gt.80.) then
		dxacl3=0.0
		return
c		write (*,*) ' adx,ady ',adx,ady
c		iflpr=1
	endif

	if(e.ne.esav.or.nzsav.ne.nz) then
		esav=e
		nzsav=nz
		do iz=1,nz
		z=z0+dz*(float(iz)-0.5)
		cl=z
		ezsav(iz)=azlngd(e,cl,dz)
		enddo
	endif


	do iz=1,nz
	z=z0+dz*(float(iz)-0.5)
	cl=z
c	ez=azlngd(e,cl,dz)
	ez=ezsav(iz)
c xz,yz gamma coord from cell center at z level
	xz=xz0+ax*z
	yz=yz0+ay*z
	dddd=dxaclz(xz,yz,dcx,dcy,z)
	dxacl3=dxacl3+ez*dddd
	enddo

	return
	end
c ===========================================================

        FUNCTION dxACLz(XG,YG,dax,day,z)
c      -------------------------------      
C       XG,YG GAMMA COOR. X=0,Y=0 - CENTER OF THE CELL
C       dx,dy - SIZE OF THE CELL (MM)
C
        EXTERNAL dxA2Fz
        X=-ABS(XG)
        Y=-ABS(YG)
        Dx=DAx/2.
        Dy=DAy/2.
        XP=X+Dx
        YP=Y+Dy
        XM=X-Dx
        YM=Y-Dy
        dxACLz=dxA2Fz(XP,YP,z)+dxA2Fz(XM,YM,z)-
     -	dxA2Fz(XP,YM,z)-dxA2Fz(XM,YP,z)
	if(xg.gt.0.) dxaclz=-dxaclz
        RETURN
        END
C ----------------------------
        FUNCTION dxa2Fz(xc,yc,z)
	ax=afz(xc,z)
	ay=afz(yc,z)
	dxax=dxafz(xc,z)
c directly
	daxa2=daa2fz(ax,ay)
c by tables
c	daxa2=daa2fzt(ax,ay)
	dxa2fz=daxa2*dxax
        RETURN
	END
c ----------------------------------------------------------
	function daa2fz(ax,ay)
C      ----------------------
C        DATA  NPAR,C0,C1,C2/3,0.34,0.00,11.00/
C  LEDNEV'S FUNCTION
        DATA  NPAR,C0,C1,C2/3,0.82,-2.2,24.30/
        DATA  CPI/6.283186/

	X=0.5-AX
	Y=0.5-AY
	PNX=X*(1.-4.*X*X)
	dxpnx=1.-12.*x*x
	PNY=Y*(1.-4.*Y*Y)
        PXY=C0+C1*(X**2+Y**2)+C2*X**2*Y**2
	dxpxy=2.*(c1+c2*y**2)*x
c	a2fz=AX*AY+PNX*PNY*PXY
	daa2fz=ay-pny*(dxpnx*pxy+pnx*dxpxy)
	return
	end
c And the same by tables  ***********************************
        function daa2fzt(ax,ay)
	common /comdaa2/nx,ny,xmi,xma,ymi,yma,daa2(100,100)
	if(ifl.eq.0) then
c read from the file
	open(20,file='daa2.ffun',form='formatted',status='unknown')
	read (20,*) nx,ny,xmi,xma,ymi,yma
	read (20,*) daa2
	close (20)

	ifl=1
	endif
	daa2fzt=f2tab(daa2,nx,ax,xmi,xma,nx,ay,ymi,yma,ny)
	return
	end
c ------------------------------------------------------------
        FUNCTION dxAFz(X,z)
C       -------------
C    SHOWER ENERGY FRACTION RIGHT FROM BOUNDARY
C    BOUNDARY COOR X=0.
	if(x.gt.0) then
		dxafz=-dxcumprf(x,z)
	else
		dxafz=-dxcumprf(-x,z)
	endif
	return
	end
c ---------------------------------------------------------
	function dxcumprf(cx,cl)
	common/shdat/ rl,zmxr,er,tmx
	clmx=rl*(tmx+2.)
	dc=clmx-cl
	dcarg=0.9*dc/clmx
c	write (*,*) ' dcarg ',dcarg
	if(dcarg.lt.-10.) then
		cfcl=0.
	elseif(dcarg.gt.10.) then
		cfcl=0.
	else
		cfcl=0.71*1.15*exp(dcarg)	
	endif
c	cfcl=1.15*exp(0.9*dc/clmx)	
c	cfcl=0.71*cfcl
	dxcumprf=cfcl*dccryex(cx*cfcl)
	return
	end
c ---------------------------------------------------------
	function dccryex(xc)
	x=abs(xc)
c directly
	dccryex=dacumex(x)
c by tables
c	dccryex=dacumext(x)
	return
	end
c ---------------------------------------------------------
	function dacumex(x)
c	data p1,p2 / 2.526, 1.104/
	data p1,p2 /1.480, 1.600/
	dacumex=0.5*exp((p1-sqrt(p1**2+4.*p2*x))/2./p2)
	dacumex=dacumex/2./p2
	dacumex=-0.5*dacumex/sqrt(p1**2+4.*p2*x)*4.*p2
	return
	end
c And the same by tables
	function dacumext(x)
	common /comdext/nx,xmi,xma,dacux(2000)
	if(ifl.eq.0) then
c read from the file
	open(20,file='dac1.ffun',form='formatted',status='unknown')
	read (20,*) nx,xmi,xma
	read (20,*) dacux
	close (20)
	ifl=1
	endif
	dacumext=f1tab(dacux,x,xmi,xma,nx)
	return
	end
c ----------------------------------------------------------
	subroutine chng(v,i,j)
	real v(10)
	c=v(i)
	v(i)=v(j)
	v(j)=c
	return
	end
C
C    FUNCTIONS EVAL. LIBRARY 
C
	FUNCTION F1TAB(F,X,XL,XR,NX)
C     GET FUNCTION FROM TABLE 
C     F - FUNCTION ARRAY DIM NX. X - ARGUMENT. XL,XR- ARG.ARRAY 
C     LINEAR INTERPOLATION
	REAL F(NX)
C
	RN=NX*(X-XL)/(XR-XL)+0.5
	N=RN
	IF(N.LT.1) THEN
		N=1
	ELSE IF(N.GT.NX-1) THEN
		N=NX-1
	ENDIF
	DX=RN-N
	F1TAB=F(N)+(F(N+1)-F(N))*DX
	RETURN
	END
C
	FUNCTION F2TAB(F,NXMX,X,XL,XR,NX,Y,YL,YR,MY)
C     GET FUNCTION FROM TABLE 
C     F - FUNCTION ARRAY DIM NXMX(NX)*MY, X,Y - ARGUMENTS
C     XL,XR,YL,YR - CLEAR 
C     LINEAR INTERPOLATION
	REAL F(NXMX,MY)
C
	if(x.lt.xl) x=xl
	if(x.gt.xr) x=xr
	RN=NX*(X-XL)/(XR-XL)+0.5
c	if(abs(rn).gt.32000.) then
c		write (*,*) ' rn ',rn
c	endif
	
	N=RN
	IF(N.LT.1) THEN
		N=1
	ELSE IF(N.GT.NX-1) THEN
		N=NX-1
	ENDIF
	DX=RN-N
C
	RM=MY*(Y-YL)/(YR-YL)+0.5

	if(abs(rm).gt.32000.) then
		write (*,*) ' rm ',rm
	endif

	M=RM
	IF(M.LT.1) THEN
		M=1
	ELSE IF(M.GT.MY-1) THEN
		M=MY-1
	ENDIF
	DY=RM-M
c	write (*,*) 'n,rn,m,rm',n,rn,m,rm
C
	F2TAB=F(N,M)+(F(N+1,M)-F(N,M))*DX+(F(N,M+1)-F(N,M))*DY
	RETURN
	END
