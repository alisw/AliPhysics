      subroutine USERGRIDevolve(xin,qin,pdf)
      implicit real*8 (a-h,o-z)
      include 'parmsetup.inc'
      double precision parm(nopmax)
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      CHARACTER*80 LINE
      dimension pdf(-6:6)
      integer init,set,i,j,k,l,nset,iset
      parameter(nhess=0)
      double precision fgrid(0:nhess,201,201,-6:6),grid(402)
      double precision agrid(201),alfas,Qalfa
      integer nacross,iparton(13)
      logical NPARTON(-6:6)
      common/ugrid/NUMX,NUMQ2,fgrid
      double precision up,upv,dn,dnv,usea,dsea,str,chm,bot,glu
      double precision qq(5),yntmp(5)
      save 
      x=xin
      q2=qin*qin 
      call getnset(iset)
      upv  = 0.0d0
      dnv  = 0.0d0
      usea = 0.0d0
      dsea = 0.0d0
      str  = 0.0d0
      sbar = 0.0d0
      chm  = 0.0d0
      cbar = 0.0d0
      bot  = 0.0d0
      bbar = 0.0d0
      top  = 0.0d0
      tbar = 0.0d0
      glu = 0.0d0
      if(name(iset)(1:len_trim(name(iset))).eq.'USERGRIDQ4') then
        ! HERE FOR QUARTIC POLYNOMIAL INTERPOLATION
          do np = 1,npartons
              if(iparton(np).eq.-6) tbar = USERGRIDQ4(x,Q2,grid,-6)
              if(iparton(np).eq.-5) bbar = USERGRIDQ4(x,Q2,grid,-5)
              if(iparton(np).eq.-4) cbar = USERGRIDQ4(x,Q2,grid,-4)
              if(iparton(np).eq.-3) sbar = USERGRIDQ4(x,Q2,grid,-3)
              if(iparton(np).eq.-2) usea = USERGRIDQ4(x,Q2,grid,-2)
              if(iparton(np).eq.-1) dsea = USERGRIDQ4(x,Q2,grid,-1)
              if(iparton(np).eq.0) glu = USERGRIDQ4(x,Q2,grid,0)
              if(iparton(np).eq.1) dnv = USERGRIDQ4(x,Q2,grid,1)
              if(iparton(np).eq.2) upv = USERGRIDQ4(x,Q2,grid,2)
              if(iparton(np).eq.3) str = USERGRIDQ4(x,Q2,grid,3)
              if(iparton(np).eq.4) chm = USERGRIDQ4(x,Q2,grid,4)
              if(iparton(np).eq.5) bot = USERGRIDQ4(x,Q2,grid,5)
              if(iparton(np).eq.6) top = USERGRIDQ4(x,Q2,grid,6)
          enddo
      else if(name(iset)(1:len_trim(name(iset))).eq.'USERGRIDQ3') then
        ! HERE FOR CUBIC POLYNOMIAL INTERPOLATION
          do np = 1,npartons
              if(iparton(np).eq.-6) tbar = USERGRIDQ3(x,Q2,grid,-6)
              if(iparton(np).eq.-5) bbar = USERGRIDQ3(x,Q2,grid,-5)
              if(iparton(np).eq.-4) cbar = USERGRIDQ3(x,Q2,grid,-4)
              if(iparton(np).eq.-3) sbar = USERGRIDQ3(x,Q2,grid,-3)
              if(iparton(np).eq.-2) usea = USERGRIDQ3(x,Q2,grid,-2)
              if(iparton(np).eq.-1) dsea = USERGRIDQ3(x,Q2,grid,-1)
              if(iparton(np).eq.0) glu = USERGRIDQ3(x,Q2,grid,0)
              if(iparton(np).eq.1) dnv = USERGRIDQ3(x,Q2,grid,1)
              if(iparton(np).eq.2) upv = USERGRIDQ3(x,Q2,grid,2)
              if(iparton(np).eq.3) str = USERGRIDQ3(x,Q2,grid,3)
              if(iparton(np).eq.4) chm = USERGRIDQ3(x,Q2,grid,4)
              if(iparton(np).eq.5) bot = USERGRIDQ3(x,Q2,grid,5)
              if(iparton(np).eq.6) top = USERGRIDQ3(x,Q2,grid,6)
          enddo  
      else if(name(iset)(1:len_trim(name(iset))).eq.'USERGRIDQ2') then
        ! HERE FOR QUADRATIC POLYNOMIAL INTERPOLATION
          do np = 1,npartons
              if(iparton(np).eq.-6) tbar = USERGRIDQ2(x,Q2,grid,-6)
              if(iparton(np).eq.-5) bbar = USERGRIDQ2(x,Q2,grid,-5)
              if(iparton(np).eq.-4) cbar = USERGRIDQ2(x,Q2,grid,-4)
              if(iparton(np).eq.-3) sbar = USERGRIDQ2(x,Q2,grid,-3)
              if(iparton(np).eq.-2) usea = USERGRIDQ2(x,Q2,grid,-2)
              if(iparton(np).eq.-1) dsea = USERGRIDQ2(x,Q2,grid,-1)
              if(iparton(np).eq.0) glu = USERGRIDQ2(x,Q2,grid,0)
              if(iparton(np).eq.1) dnv = USERGRIDQ2(x,Q2,grid,1)
              if(iparton(np).eq.2) upv = USERGRIDQ2(x,Q2,grid,2)
              if(iparton(np).eq.3) str = USERGRIDQ2(x,Q2,grid,3)
              if(iparton(np).eq.4) chm = USERGRIDQ2(x,Q2,grid,4)
              if(iparton(np).eq.5) bot = USERGRIDQ2(x,Q2,grid,5)
              if(iparton(np).eq.6) top = USERGRIDQ2(x,Q2,grid,6)
          enddo  
      else
          print *,'Unknown interpolation method ',   &
     & name(iset)(1:len_trim(name(iset))),' called for!'
          stop
      endif
      up = upv + usea
      dn = dnv + dsea  

      if(NPARTON(-6)) then
          pdf(-6) = tbar
      else 
          pdf(-6) = top
      endif
      pdf(6) = top

      if(NPARTON(-5)) then
          pdf(-5) = bbar
      else
          pdf(-5) = bot
      endif
      pdf(5) = bot

      if(NPARTON(-4)) then
          pdf(-4) = cbar
      else
          pdf(-4) = chm
      endif
      pdf(4) = chm

      if(NPARTON(-3)) then
          pdf(-3) = sbar
      else
          pdf(-3) = str
      endif
      pdf(3) = str
      pdf(-2) = usea
      pdf(2)  = up
      pdf(-1) = dsea
      pdf(1)  = dn
      pdf(0)  = glu
      return
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry USERGRIDread(nset)
!     
      call getnmem(nset,imem) 

      read(1,*)nmem(nset),ndef(nset),NUMX,NUMQ2,nacross,npartons,(iparton(k),k=1,npartons)
      
      do np = -6,6
          NPARTON(np) = .FALSE.
      enddo
      
      nnx = numx/nacross
      nnq2 = numq2/nacross
      
      do iq=1,nnq2
         read(1,*) (grid(NUMX+(iq-1)*nacross+ii),ii=1,nacross)
      enddo

      do jx=1,nnx
         read(1,*) (grid((jx-1)*nacross+ii),ii=1,nacross)
      enddo

    !read in to alphas grid  
      do iq=1,nnq2
         read(1,*) (agrid((iq-1)*nacross+ii),ii=1,nacross)
      enddo

      do ns=0,nmem(nset)
        do k = 1,npartons
          NPARTON(iparton(k)) = .TRUE.
          do IQ=1,NUMQ2
            do JX=1,nnx
	      read(1,*)(fgrid(ns,(jx-1)*nacross+ii,iq,iparton(k)),ii=1,nacross)
            enddo
          enddo
        enddo
      enddo
      return
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry USERGRIDalfa(alfas,Qalfa)
      q2 = Qalfa*Qalfa
      nqlow = -1
      nqhi = NUMQ2+1
      do while (nqhi-nqlow.gt.1)
         nqmid = (nqlow+nqhi)/2
	 if(q2.ge.grid(NUMX+nqmid)) then
	    nqlow = nqmid
	 else 
	    nqhi = nqmid
	 endif
      enddo
      ! do a linear interpolation on the alphas grid
      npt = 2 
      call userpolint(grid(NUMX+nqlow),agrid(nqlow),npt,q2,alfas,dy)
      call getnset(iset)
      call getnmem(iset,imem)
      call listPDF(iset,imem,parm)
      alfas=alfas*parm(1)/0.1176d0
      return
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry USERGRIDinit(Eorder,Q2fit)
      return
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry USERGRIDpdf(mem)
      imem = mem
      call getnset(iset)
      call setnmem(iset,imem)
      return
!
 1000 format(5e13.5)
      end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! attempt at polynomial (order 4) interpolation based on  polint
      double precision function USERGRIDQ4(x,Q2,grid,n)
      implicit real*8(a-h,o-z)
      integer iset,imem,nhess
      parameter(nhess=0)
      double precision grid(402),x,Q2
      double precision fgrid(0:nhess,201,201,-6:6)
      double precision ya(5,5),yntmp(5),ymtmp(5)
      common/ugrid/NUMX,NUMQ2,fgrid
      
      npt = 5
      call getnset(iset)
      call getnmem(iset,imem)      
      
! find the x bins around x
      nxlow = -1
      nxhi = NUMX+1
      do while (nxhi-nxlow.gt.1)
         nxmid = (nxlow+nxhi)/2
	 if(x.ge.grid(nxmid)) then
	    nxlow = nxmid
	 else 
	    nxhi = nxmid
	 endif
      enddo

! find the q2 bins around q2
      nqlow = -1
      nqhi = NUMQ2+1
      do while (nqhi-nqlow.gt.1)
         nqmid = (nqlow+nqhi)/2
	 if(q2.ge.grid(NUMX+nqmid)) then
	    nqlow = nqmid
	 else 
	    nqhi = nqmid
	 endif
      enddo
      
! fill the temp 4x4 funtion array (allowing for endpoints and extrapolation)
	if(nxlow.le.0) nxlow=1
	if(nxlow.ge.NUMX) nxlow=NUMX-1
        if(nxlow.eq.1) then
	   nxbot = 1
	elseif(nxlow.eq.2) then
	   nxbot = 2
	else
	   if(nxlow.eq.NUMX-1) then
	      nxbot = 4
	   else
	      nxbot = 3
	   endif
	endif
	if(nqlow.le.0) nqlow=1
	if(nqlow.ge.NUMQ2) nqlow=NUMQ2-1
        if(nqlow.eq.1) then
	   nqbot = 1
	elseif(nqlow.eq.2) then
	   nqbot = 2
	else
	   if(nqlow.eq.NUMQ2-1) then
	      nqbot = 4
	   else   
	      nqbot = 3
	   endif
	endif
	
        do nx=1,5
          do nq=1,5 
	      ya(nx,nq) = fgrid(imem,nxlow+nx-nxbot,nqlow+nq-nqbot,n)
	  enddo   
        enddo      
      
      do j=1,5
        do k=1,5
	      yntmp(k)=ya(j,k)
	    enddo
	    call userpolint(grid(NUMX+nqlow-nqbot+1),yntmp,npt,q2,ymtmp(j),dy)
      enddo 
      call userpolint(grid(nxlow-nxbot+1),ymtmp,npt,x,y,dy)
      
      usergridq4=y
      
      return
      end     
!=========================================================
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! attempt at polynomial (order 3) interpolation based on  polint
      double precision function USERGRIDQ3(x,Q2,grid,n)
      implicit real*8(a-h,o-z)
      integer iset,imem,nhess
      parameter(nhess=0)
      double precision grid(402),x,Q2
      double precision fgrid(0:nhess,201,201,-6:6)
      double precision ya(4,4),yntmp(4),ymtmp(4)
      common/ugrid/NUMX,NUMQ2,fgrid
      
      npt = 4
      call getnset(iset)
      call getnmem(iset,imem)      
      
! find the x bins around x
      nxlow = -1
      nxhi = NUMX+1
      do while (nxhi-nxlow.gt.1)
         nxmid = (nxlow+nxhi)/2
	 if(x.ge.grid(nxmid)) then
	    nxlow = nxmid
	 else 
	    nxhi = nxmid
	 endif
      enddo

! find the q2 bins around q2
      nqlow = -1
      nqhi = NUMQ2+1
      do while (nqhi-nqlow.gt.1)
         nqmid = (nqlow+nqhi)/2
	 if(q2.ge.grid(NUMX+nqmid)) then
	    nqlow = nqmid
	 else 
	    nqhi = nqmid
	 endif
      enddo
      
! fill the temp 4x4 funtion array (allowing for endpoints and extrapolation)
	if(nxlow.le.0) nxlow=1
	if(nxlow.ge.NUMX) nxlow=NUMX-1
    if(nxlow.eq.1) then
	   nxbot = 1
	else
	   if(nxlow.eq.NUMX-1) then
	      nxbot = 3
	   else
	      nxbot = 2
	   endif
	endif
	if(nqlow.le.0) nqlow=1
	if(nqlow.ge.NUMQ2) nqlow=NUMQ2-1
    if(nqlow.eq.1) then
	   nqbot = 1
	else
	   if(nqlow.eq.NUMQ2-1) then
	      nqbot = 3
	   else   
	      nqbot = 2
	   endif
	endif
	
        do nx=1,4
          do nq=1,4 
	      ya(nx,nq) = fgrid(imem,nxlow+nx-nxbot,nqlow+nq-nqbot,n)
	  enddo   
        enddo      
      
      do j=1,4
        do k=1,4
	      yntmp(k)=ya(j,k)
	    enddo
	    call userpolint(grid(NUMX+nqlow-nqbot+1),yntmp,npt,q2,ymtmp(j),dy)
      enddo 
      call userpolint(grid(nxlow-nxbot+1),ymtmp,npt,x,y,dy)
      
      usergridq3=y
      
      return
      end     
!=========================================================
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! attempt at polynomial (order 2) interpolation based on  polint
      double precision function USERGRIDQ2(x,Q2,grid,n)
      implicit real*8(a-h,o-z)
      integer iset,imem,nhess
      parameter(nhess=0)
      double precision grid(402),x,Q2
      double precision fgrid(0:nhess,201,201,-6:6)
      double precision ya(3,3),yntmp(3),ymtmp(3)
      common/ugrid/NUMX,NUMQ2,fgrid
      
      npt = 3
      call getnset(iset)
      call getnmem(iset,imem)      
      
! find the x bins around x
      nxlow = -1
      nxhi = NUMX+1
      do while (nxhi-nxlow.gt.1)
         nxmid = (nxlow+nxhi)/2
	 if(x.ge.grid(nxmid)) then
	    nxlow = nxmid
	 else 
	    nxhi = nxmid
	 endif
      enddo

! find the q2 bins around q2
      nqlow = -1
      nqhi = NUMQ2+1
      do while (nqhi-nqlow.gt.1)
         nqmid = (nqlow+nqhi)/2
	 if(q2.ge.grid(NUMX+nqmid)) then
	    nqlow = nqmid
	 else 
	    nqhi = nqmid
	 endif
      enddo
      
! fill the temp 4x4 funtion array (allowing for endpoints and extrapolation)
	if(nxlow.le.0) nxlow=1
	if(nxlow.ge.NUMX) nxlow=NUMX-1
    if(nxlow.eq.1) then
	   nxbot = 1
	else
	   if(nxlow.eq.NUMX-1) then
	      nxbot = 2
	   else
	      nxbot = 1
	   endif
	endif
	if(nqlow.le.0) nqlow=1
	if(nqlow.ge.NUMQ2) nqlow=NUMQ2-1
    if(nqlow.eq.1) then
	   nqbot = 1
	else
	   if(nqlow.eq.NUMQ2-1) then
	      nqbot = 2
	   else   
	      nqbot = 1
	   endif
	endif
	
        do nx=1,3
          do nq=1,3 
	      ya(nx,nq) = fgrid(imem,nxlow+nx-nxbot,nqlow+nq-nqbot,n)
	  enddo   
        enddo      
      
      do j=1,3
        do k=1,3
	      yntmp(k)=ya(j,k)
	    enddo
	    call userpolint(grid(NUMX+nqlow-nqbot+1),yntmp,npt,q2,ymtmp(j),dy)
      enddo 
      call userpolint(grid(nxlow-nxbot+1),ymtmp,npt,x,y,dy)
      
      usergridq2=y
      
      return
      end     
!=========================================================
      
         SUBROUTINE USERPOLINT (XA,YA,N,X,Y,DY)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!                                        Adapted from "Numerical Recipes"
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) RETURN
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
!=========================================================
