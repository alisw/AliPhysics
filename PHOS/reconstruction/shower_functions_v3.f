c ==================================================================

	function coornew(el,er,angl,mmm)
	mmnew=mmconvert(mmm)	
	esum=el+er
	ex=el/esum

	if(esum.le.0.) then
		write (*,*) ' coornew diagnostic bad esum=',esum
		return
	endif
	
	itip=newtipcell(mmnew)
	if(itip.eq.1) then
		xx=xcryexgl(ex)
	elseif(itip.eq.2) then
c		xx=xcryexcr(ex)
		xx=xcryexcrangl(ex,angl)
	endif
	coornew=xx
	return
	end	
	
	function xcryexgl(ar)
c FOR LEAD GLASS
	data p1,p2 /1.480, 1.600/
c ar - energy fraction right from boundary
	

	if(ar.gt.0.5) then
		a=1.-ar
		s=-1.
	else
		a=ar
		s=1.
	endif
	if(a.le.0.) a=0.0001

	y=alog(2.*a)
	x=p2*y**2-p1*y
c lead glass for my geant
c	x=x/0.55
c lead glass for my geant 2*2 cells
	x=x/0.65
c corrections for crystalls 2*2, systematic is better than 0.25 mm for 3Gev geant showers
c94	x=x*0.54-0.25
c end corrections for crystalls
	xcryexgl=s*x
	return
	end

	function xcryexcr(ar)
c FOR CRYSTALS
c	data p1,p2 /1.480, 1.600/
c  very new 26.02.98
	data p1,p2 /1.700, 1.150/
c ar - energy fraction right from boundary
	if(ar.gt.0.5) then
		a=1.-ar
		s=-1.
	else
		a=ar
		s=1.
	endif
	if(a.le.0.) a=0.0001
	y=alog(2.*a)
	x=p2*y**2-p1*y
	xcryexcr=s*x
	return
	end

	function xcryexcrangl(ar,angl)
c FOR CRYSTALS ALICE parametrisation 23.02.98
c	data p1,p2 /1.480, 1.600/
	data p1,p2 /1.700, 1.150/
c ar - energy fraction right from boundary
        aangl=abs(angl)
c        write (*,*) ' ar,angl ',ar,angl
        if(angl.ge.0.) then
                pp1=p1
                pp2=p2-0.7*aangl
        else
                pp1=p1+7.6*aangl
                pp2=p2-0.9*aangl
        endif
c        if(angl.ge.0.) then
c                pp1=p1+7.6*aangl
c                pp2=p2-0.9*aangl
c        else
c                pp1=p1
c                pp2=p2-0.7*aangl
c        endif
	if(ar.gt.0.5) then
		a=1.-ar
		s=-1.
	else
		a=ar
		s=1.
	endif
	if(a.le.0.) a=0.0001
	y=alog(2.*a)
	x=pp2*y**2-pp1*y
	xcryexcrangl=s*x
	return
	end
c ==============================================================
	real function cumulnew(xc,itip)
        real xc
        integer itip

	if(itip.eq.1) then
		a=cucryexgl(xc)
	elseif(itip.eq.2) then
		a=cucryexcr(xc)
	endif
        cumulnew=a
        return
        end
        
	real function dcumulnew(xc,itip)
        real xc
        integer itip

	if(itip.eq.1) then
		a=dcucryexgl(xc)
	elseif(itip.eq.2) then
		a=dcucryexcr(xc)
	endif
        dcumulnew=a
        return
        end
        
	real function cucryexcr(xc)
c fit for ALICE 94 RUN DATA
c	data p1,p2 /2.526, 1.104/
c94	data p1,p2 /2.300, 2.044/
c crystals for my geant
c	data p1,p2 /1.480, 1.600/
c  very new 26.02.98
	data p1,p2 /1.700, 1.150/
	x=abs(xc)
	a=0.5*exp((p1-sqrt(p1**2+4.*p2*x))/2./p2)
	if(xc.ge.0.) then
		cucryexcr=1.-a
	else
		cucryexcr=a
	endif
	return
	end

	real function dcucryexcr(xc)
c fit for ALICE 94 RUN DATA
c	data p1,p2 /2.526, 1.104/
c	data p1,p2 /2.300, 2.044/
c crystals for my geant
c	data p1,p2 /1.480, 1.600/
c  very new 26.02.98
	data p1,p2 /1.700, 1.150/
	x=abs(xc)

	rr=sqrt(p1**2+4.*p2*x)
	da=0.5*exp((p1-rr)/2./p2)/rr
	dcucryexcr=da
	
	return
	end

	real function cucryexgl(xc)
c fit for ALICE 94 RUN DATA
c	data p1,p2 /2.526, 1.104/
c94	data p1,p2 /2.300, 2.044/
c crystals for my geant
	data p1,p2 /1.480, 1.600/
	x=abs(xc)
c lead glass for my geant
	x=x*0.55

	a=0.5*exp((p1-sqrt(p1**2+4.*p2*x))/2./p2)
	if(xc.ge.0.) then
		cucryexgl=1.-a
	else
		cucryexgl=a
	endif
	return
	end

	real function dcucryexgl(xc)
c fit for ALICE 94 RUN DATA
c	data p1,p2 /2.526, 1.104/
c	data p1,p2 /2.300, 2.044/
c crystals for my geant
	data p1,p2 /1.480, 1.600/
	x=abs(xc)
c lead glass for my geant
	x=x*0.55

	rr=sqrt(p1**2+4.*p2*x)
	da=0.5*exp((p1-rr)/2./p2)/rr
c lead glass for my geant
	da=da*0.55
	dcucryexgl=da
	
	return
	end
c ==============================================================

	real function ampcelnew(e,vx,vc,vb,mmm)
	common /comkey/ keykey
	real vx(5),vc(3),vb(3)
        integer mmm

	mmnew=mmconvert(mmm)	
	itip=newtipcell(mmnew)
	xg=vx(1)
	yg=vx(2)
	xc=vc(1)
	yc=vc(2)
	if(keykey.eq.0) then
        	a=aclnew(xg-xc,yg-yc,vb,itip)
	elseif(keykey.eq.1) then
c new with angle
        	a=acl3(e,vx,vc,vb)
	else
		write (*,*) ' Net takoi keykey !! ',keykey
                a=0.
	endif	
        ampcelnew=a
        return
        end

	real function dampcelnew(nd,e,vx,vc,vb,mmm)
	common /comkey/keykey
	real vx(5),vc(3),vb(3)
        integer mmm

	mmnew=mmconvert(mmm)	
	itip=newtipcell(mmnew)
c old no angle
	if(keykey.eq.0) then
		xg=vx(1)
		yg=vx(2)
		xc=vc(1)
		yc=vc(2)
	   if(nd.eq.1) then
		da=dxaclnew(xg-xc,yg-yc,vb,itip)
	   elseif(nd.eq.2) then
		da=dxaclnew(yg-yc,xg-xc,vb,itip)
	   else
		write (*,*) ' net takoi nd ',nd
	   endif
	elseif(keykey.eq.1) then
c new with angle
		da=dacl3(nd,e,vx,vc,vb)
	else
		write (*,*) ' net takoi keykey ',keykey
                da=0.
	endif
        dampcelnew=da

	return
	end

        FUNCTION aclnew(XG,YG,DA,mmn)
C
C     SHOWER FRACTION IN THE CELL
C       XG,YG GAMMA COOR. X=0,Y=0 - CENTER OF THE CELL
C       DA - SIZE OF THE CELL (MM)
C
	real da(3)
        integer mmn
        
        EXTERNAL a2fnew
C       SX=SIGN(1.,XG)
C       SY=SIGN(1.,YG)    FOR ANGLE
        X=-ABS(XG)
        Y=-ABS(YG)
        Dx=DA(1)/2.
	dy=da(2)/2.
        XP=X+Dx
        YP=Y+Dy
        XM=X-Dx
        YM=Y-Dy
        aclnew=a2fnew(XP,YP,mmn)+a2fnew(XM,YM,mmn)-
     -         a2fnew(XP,YM,mmn)-a2fnew(XM,YP,mmn)
c	write (*,*) ' xg,yg,dx,dy,acl ',xg,yg,dx,dy,acl
        RETURN
        END
C
        FUNCTION dxaclnew(XG,YG,DA,mmn)
C
C       DXACL=D(ACL(X,Y,DA))/D(X)
C
        integer mmn
        EXTERNAL dxa2fnew
        SX=SIGN(1.,XG)
C       SY=SIGN(1.,YG)    FOR ANGLE
        X=-ABS(XG)
        Y=-ABS(YG)
        D=DA/2.
        XP=X+D
        YP=Y+D
        XM=X-D
        YM=Y-D
        A=dxa2fnew(XP,YP,mmn)+dxa2fnew(XM,YM,mmn)-
     -    dxa2fnew(XP,YM,mmn)-dxa2fnew(XM,YP,mmn)
        dxaclnew=-SX*A
        RETURN
        END

        FUNCTION a2fnew(X,Y,mmn)
        integer mmn
        EXTERNAL afnew
        AX=afnew(X,mmn)
        AY=afnew(Y,mmn)
	a2fnew=aa2fz(ax,ay)
	return
        END

        FUNCTION dxa2fnew(X,Y,mmn)
        integer mmn
        EXTERNAL afnew,dafnew
        AX=afnew(X,mmn)
        AY=afnew(Y,mmn)
	dxa2fnew=daa2fz(ax,ay)*dafnew(x,mmn)
	return
        END

        FUNCTION afnew(x,mmn)
	afnew=cumulnew(x,mmn)
	return
        end
        
        FUNCTION dafnew(x,mmn)
	dafnew=dcumulnew(x,mmn)
	return
        end

	function sigwlgam0(ewl,egam)
c        write (*,*) ' ewl,egam ',ewl,egam 
	if(egam.le.0.) then
		sigwlgam0=0.
	else
		sqwl=sqrt(ewl)
		fr=ewl/egam
		if(fr.gt.0.8) fr=0.8
c		sigwlgam0=0.36*(1.-fr/0.88)*sqwl
        xp=abs(fr-0.05)
        sfun=0.5*sin(3.14/0.7*xp)+0.5*sin(3.14/0.76*sqrt(xp*0.7))
c        sfun=2.5*egam*sfun+0.13-0.13*xp
c        fm=(1.+exp(-egam/0.2))
c        if(fm.lt.1.) fm=1.
c        if(fm.gt.2.) fm=2.
        fm=2.5
        sfun=egam*fm*sfun+0.13-0.13*xp
c        sigwlgam0=sqrt(egam)*sfun*sqwl
        sigwlgam0=0.1/sqrt(egam)*sfun*sqwl
c        write (*,*) ' fr sfun ',fr,sfun 
c        write (*,*) ' sigwlgam0 ',sigwlgam0,ewl,egam 

	endif
	return
	end	

	function dispwlm(ewl,m)

        common /comsigma/sigmaph,sigmapd,sigphsq,sigpdsq


	dispph=sigphsq*ewl
	disppd=sigpdsq	
	dispwlm=dispph+disppd

	return
	end

        real function zmidlshower(ee,m)
	include 'comgam.for'
        zmidlshower = 46.
c        zmidlshower = 50.
        return
        end


