        
        subroutine hisgamnew

	include 'comgam.for'

        nh0=110000
	if(ifl.eq.0) then
        ngmax=200
        egmax=5.
        etmax=150.
        ptmax=5.
        smax=1500.
        fngmax=float(ngmax)
c		call hbook1(nh0+1,' ngam$',ngmax,0.,fngmax,0)
c		call hbook1(nh0+2,' egam$',100,0.,egmax,0)
c		call hbook1(nh0+3,' etot$',150,0.,etmax,0)
c		call hbook1(nh0+4,'  pt  $',100,0.,ptmax,0)
c	call hbook2(nh0+5,'  x*y  $',50,-smax,smax,50,-smax,smax,0)
c	        call hbook1(nh0+6,'  x  $',500,-smax,smax,0)
c	        call hbook1(nh0+7,'  y  $',500,-smax,smax,0)
c		call hbook1(nh0+8,'  ptalice  $',100,0.,ptmax,0)
	ifl=1
	endif

c        call hf1(nh0+1,float(kg),1.)
        etot=0.
        do k=1,kg
             etot=etot+e(k)
        pt=sqrt(p4(1,k)**2+p4(2,k)**2)
        ptalice=sqrt(p4(1,k)**2+p4(3,k)**2)
c             call hf1(nh0+2,e(k),1.)
c             call hf2(nh0+5,xw(k),yw(k),1.)
c             call hf1(nh0+6,xw(k),1.)
c             call hf1(nh0+7,yw(k),1.)
c             call hf1(nh0+4,pt,1.)
c             call hf1(nh0+8,ptalice,1.)
         enddo
c        call hf1(nh0+3,etot,1.)

        return
        end

        subroutine p4gamnew
        common /comgeom/igeomflag

	include 'comgam.for'
	include 'comwlgen.for'
        external xgamvert,ygamvert,zgamvert

        do k=1,kg
        iv=igamvert(k)
        xv=xgamvert(iv)
        yv=ygamvert(iv)
        zv=zgamvert(iv)
        igd=igdev(k)
        if(igeomflag.eq.2) then
                rg=zgdev(k)+xyzwall(3,igd)
                rr=4600.
                fi=xw(k)/rr
                xg=rg*sin(fi)-xv
                zg=rg*cos(fi)-zv
                yg=yw(k)-yv
        else
        xg=xw(k)-xv
        yg=yw(k)-yv
c        write (*,*) ' k zgdev xyzwall zv',k,zgdev(k),xyzwall(3,igd),zv
        zg=zgdev(k)+xyzwall(3,igd)-zv
        endif
        rr=sqrt(xg**2+yg**2+zg**2)
        p4(1,k)=e(k)*xg/rr
        p4(2,k)=e(k)*yg/rr
        p4(3,k)=e(k)*zg/rr
        p4(4,k)=e(k)
        enddo

        return
        end

        real function xgamvert(ivert)
	include 'comgam.for'
        if(ivert.le.0) then
                xgamvert=0.
        elseif(ivert.gt.ntotvertx) then
                xgamvert=0.
                write (*,*) ' wrong vertex number '
        else
                xgamvert=gamvertex(1,ivert)
        endif
        return
        end

        real function ygamvert(ivert)
	include 'comgam.for'
        if(ivert.le.0) then
                ygamvert=0.
        elseif(ivert.gt.ntotvertx) then
                ygamvert=0.
                write (*,*) ' wrong vertex number '
        else
                ygamvert=gamvertex(2,ivert)
        endif
        return
        end

        real function zgamvert(ivert)
	include 'comgam.for'
        if(ivert.le.0) then
                zgamvert=0.
        elseif(ivert.gt.ntotvertx) then
                zgamvert=0.
                write (*,*) ' wrong vertex number '
        else
                zgamvert=gamvertex(3,ivert)
        endif
        return
        end

        subroutine clrcompart
	include 'compart.for'
        nbpart=0
        return
        end
        
        subroutine addp4tocompart(ncrad)

        common /comgeom/igeomflag
        include 'event_format.inc'
	include 'comgam.for'
	include 'compart.for'

        data ifl/0/
        real fir(10)
        data fir/10*0./

        if(ifl.eq.0) then
        fi0=22.*44./4600.
        fi0=atan(fi0)
        nmmmax=crystals_matrix_amount_PHOS
        do i=1,nmmmax
cnew
        fir(i)=fi0*float(nmmmax-1)-2.*fi0*float(i-1)
        enddo
        ifl=1
        endif

        do k=1,kg
        if(nbpart+1.gt.nbpartmax) then
                write (*,*) ' incr nbpartmax ',nbpartmax
                return
        endif
        nbpart=nbpart+1

        if(igeomflag.eq.1) then
                call rotrefr(p4(1,k),ppart(1,nbpart),fir(ncrad))
        else
                call rotrefc(p4(1,k),ppart(1,nbpart))
        endif

        enddo
        return
        end

        subroutine extrp4wrf(ncrad)

        common /comgeom/igeomflag
        include 'event_format.inc'
	include 'comgam.for'
	include 'compart.for'
	include 'comggen.for'

        data ifl/0/
        real fir(10)
        data fir/10*0./
	integer ijdd(2)
        real p4temp(4)


        if(ifl.eq.0) then
        fi0=22.*44./4600.
        fi0=atan(fi0)
        nmmmax=crystals_matrix_amount_PHOS
        do i=1,nmmmax
        fir(i)=fi0*float(nmmmax-1)-2.*fi0*float(i-1)
        enddo
        ifl=1
        endif

c        do n=1,nbpart
c        if(igeomflag.eq.1) then
c                call rotrefr_a(p4temp,ppart(1,n),fir(ncrad))
c        else
c                call rotrefc_a(p4temp,ppart(1,n))
c        endif
c
c
c        call xydetwlnew(xx,yy,p4temp)
c	call mdetwlnew(mmww,xx,yy,ijdd,1)
c
c        enddo

        do n=1,ngamgen
        if(igeomflag.eq.1) then
                call rotrefr_a(p4temp,g(1,n),fir(ncrad))
        else
                call rotrefc_a(p4temp,g(1,n))
        endif


        call xydetwlnew(xx,yy,p4temp)
        write (*,*) ' xx,yy ',xx,yy,xgamgen(n),ygamgen(n)
c	call mdetwlnew(mmww,xx,yy,ijdd,1)

       enddo

        return
        end

        subroutine calcgamgennew(ncrad)

        common /comgeom/igeomflag
        include 'event_format.inc'
	include 'comgam.for'
	include 'compart.for'
	include 'comggen.for'

        data ifl/0/
        real fir(10)
        data fir/10*0./
	integer ijdd(2)
        real p4temp(4)


        if(ifl.eq.0) then
        fi0=22.*44./4600.
        fi0=atan(fi0)
        nmmmax=crystals_matrix_amount_PHOS
        do i=1,nmmmax
        fir(i)=fi0*float(nmmmax-1)-2.*fi0*float(i-1)
        enddo
        ifl=1
        endif

        do n=1,ngamgen
	call ucopy(g(1,n),p4ggen(1,n),3)
	p4ggen(4,n)=sqrt(g(1,n)**2+g(2,n)**2+g(3,n)**2)
	egamgen(n)=p4ggen(4,n)
        if(igeomflag.eq.1) then
                call rotrefr_a(p4temp,p4ggen(1,n),fir(ncrad))
        else
                call rotrefc_a(p4temp,p4ggen(1,n))
        endif


        call xydetwlnew(xx,yy,p4temp)
        xgamgen(n)=xx
        ygamgen(n)=yy

       enddo

	do n=1,ngamgen
	rmingen(n)=10000.
	kmingen(n)=n
	enddo

	do n1=1,ngamgen-1
	do n2=n1+1,ngamgen
	rx=xgamgen(n1)-xgamgen(n2)
	ry=ygamgen(n1)-ygamgen(n2)
	rr=sqrt(rx**2+ry**2)
	if(rr.lt.rmingen(n1)) then
		rmingen(n1)=rr
		kmingen(n1)=n2
	endif	
	if(rr.lt.rmingen(n2)) then
		rmingen(n2)=rr
		kmingen(n2)=n1
	endif	

	enddo
	enddo

        return
        end

        subroutine xydetwlnew(x,y,p4)
        common /comgeom/igeomflag
        real x,y,p4(4)

c	zc=4600.+46.
	zc=4600.+40.
        if(igeomflag.eq.2) then
c cylindric case
                rcentr=4600.
        	fi1=atan(p4(1)/p4(3))
                g12=sqrt(p4(3)**2+p4(1)**2)
        	x=fi1*rcentr
        	y=p4(2)/g12*zc
        elseif(igeomflag.eq.1) then
c rect case
        	x=p4(1)/p4(3)*zc
        	y=p4(2)/p4(3)*zc
        else
                write (*,*) ' check geometry flag '
                stop
        endif
        return
        end



        subroutine rotrefr(p1,p2,fi)
        real p1(4),p2(4)

        cosfi=cos(fi)
        sinfi=sin(fi)
        xp=p1(3)
        yp=p1(1)
        p2(1)=xp*cosfi+yp*sinfi
        p2(2)=yp*cosfi-xp*sinfi
        p2(3)=p1(2)
        p2(4)=p1(4)
        return
        end

        subroutine rotrefc(p1,p2)
        real p1(4),p2(4)
        p2(1)=p1(3)
        p2(2)=p1(1)
        p2(3)=p1(2)
        p2(4)=p1(4)
        return
        end
        
        subroutine rotrefr_a(p1,p2,fi)
        real p1(4),p2(4)

        cosfi=cos(fi)
        sinfi=sin(fi)

        p1(3)=p2(1)*cosfi-p2(2)*sinfi
        p1(1)=p2(2)*cosfi+p2(1)*sinfi

        p1(2)=p2(3)
        p1(4)=p2(4)
        return
        end

        subroutine rotrefc_a(p1,p2)
        real p1(4),p2(4)
        p1(3)=p2(1)
        p1(1)=p2(2)
        p1(2)=p2(3)
        p1(4)=p2(4)
        return
        end
        
        
        
