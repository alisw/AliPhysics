c *******************     INIWL          ********************************

        SUBROUTINE iniwl(s1gev,ampshum)

        common /comgeom/igeomflag

        lun=99
        open (lun,FILE='phos_geom_r_new.dat',FORM='FORMATTED',STATUS
     $       ='OLD')
        call rgeomnew(lun)

c   general walls information & geometry
	call rwlgen
	call rallwl
	call addwsh
c coefficients
	iw=1
	cfdum=1.
	call rcoef(iw,nbcff,cfdum)

c Initialisation PHOS parameters
        call inisigma(s1gev,ampshum)

	return
        end

	subroutine rwlgen
        common /comgeom/igeomflag
	include 'comwlgen.for'

	madwl(1)=0

	nwall=1

	idw1=104
        if(igeomflag.eq.2) then
        	jdw1=352
        else
        	jdw1=88
        endif

	nwll(1)=idw1*jdw1
	madwl(1+1)=madwl(1)+nwll(1)	

	idimw(1)=idw1
	jdimw(1)=jdw1
	nspace(1)=0	

	dimblc(1,1)=22.
	dimblc(2,1)=22.
	dimblc(3,1)=200.

	xyzwall(1,1)=0.
	xyzwall(2,1)=0.
	xyzwall(3,1)=0.

	nhwl(1)=0

	dwx=dimblc(1,1)*float(idimw(1))
	dwy=dimblc(2,1)*float(jdimw(1))
	xlbon(1,1)=-dwx/2.+xyzwall(1,1)
	xrbon(1,1)=dwx/2.+xyzwall(1,1)
	xlbon(2,1)=-dwy/2.+xyzwall(2,1)
	xrbon(2,1)=dwy/2.+xyzwall(2,1)

	return
	
	end

        SUBROUTINE addwsh
        common /comgeom/igeomflag

	include 'comalwl.for'
	include 'comwlgen.for'

	xyzwall(1,1)=0.
	xyzwall(2,1)=0.
	xyzwall(3,1)=4600.

        if(igeomflag.eq.1) then
c	nw=1
c        DO I=1,NT
c	if(mm(i).gt.madwl(nw+1)) nw=nw+1
c        x(1,I)=x(1,I)+XYZWALL(1,NW)
c        x(2,I)=x(2,I)+XYZWALL(2,NW)
c        x(3,I)=x(3,I)+XYZWALL(3,NW)
c        ENDDO
        else
        endif
        
	return
	end

        SUBROUTINE rallwl

        common /comgeom/igeomflag
	include 'comalwl.for'
	include 'comwlgen.for'

	lun=20

        OPEN (lun,FILE='geomf_old.dat',FORM='FORMATTED',status='old')
c	WRITE (*,*) 'READ WALL INF from geomf_old.dat'

        READ  (lun,*) NT
c	write (*,*) ' ngl ',nt
        READ  (lun,*) (MM(I),I=1,NT)
        READ  (lun,*) (D(1,I),I=1,NT)
        READ  (lun,*) (D(2,I),I=1,NT)
        READ  (lun,*) (D(3,I),I=1,NT)
        READ  (lun,*) (X(1,I),I=1,NT)
        READ  (lun,*) (x(2,I),I=1,NT)
        READ  (lun,*) (x(3,I),I=1,NT)
        READ  (lun,*) (NS(I),I=1,NT)
        READ  (lun,*) (NKS(I),I=1,NT)
        READ  (lun,*) ((KS(J,I),J=1,10),I=1,NT)
             CLOSE (lun)
        RETURN
        END


        SUBROUTINE RCOEF(jw,ng,dum)

        common /coefil/ filcff(4)

	include 'comalwl.for'
	include 'comwlgen.for'

	real dum(10)

c	WRITE (*,*) ' IN THIS VERSION ALL COEFFICIENTS SET TO 1. '

	iw=jw
	i0w=madwl(iw)
	nggw=madwl(iw+1)-madwl(iw)
	do i=1,nggw
	cf(i0w+i)=dum(iw)
	enddo		

	return
	end

        subroutine inisigma(sigph,sigpd)

        common /comdevpar/sigmaph,sigmapd,sigphsq,sigpdsq

        sigmaph=sigph
        sigmapd=sigpd

	sigphsq=sigmaph**2
	sigpdsq=sigmapd**2

c        write (*,*) ' phe statistic ',sigph,' pedestal noise ',sigpd
        
        return
        end


        subroutine rgeomnew(lun)
	include 'comarray.for'

        read (lun,*) ncells

        do n=1,ncells
        read (lun,*) idcells(n),idcelmat(n),(cellsize(k,n),k=1,3)
        enddo 

        read (lun,*) narray

        do n=1,narray
        read (lun,*) idarray(n),idcellar(n),isizarray(n),nholes(n),
     ,                (ijarray(k,n),k=1,2),(sizarray(k,n),k=1,3)

        do nh=1,nholes(n)
        read (lun,*) (poshole(k,nh,n),k=1,3),(sizhole(k,nh,n),k=1,3)
        enddo
        enddo

        read (lun,*) ncompose

        do n=1,ncompose
        read (lun,*) ncomparr(n),mmincomp(n),mmaxcomp(n),
     ,                (poscomp(k,n),k=1,3)
        do na=1,ncomparr(n)
        read (lun,*) idcomparr(na,n),mmaxcompar(na,n),
     ,  mmincompar(na,n),(poscompar(k,na,n),k=1,3)              
        enddo
        enddo

        read (lun,*) ipmmm

        do n=1,ipmmm
        read (lun,*) mmcells(n),neibortyp(n),
     ,                (cellpos(k,n),k=1,3)
        enddo
        
        close (lun)
        return
        end

        



