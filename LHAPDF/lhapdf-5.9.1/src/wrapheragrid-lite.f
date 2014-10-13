      subroutine HERAGRIDevolve(xin,qin,pdf)
      implicit real*8 (a-h,o-z)
      include 'parmsetup.inc'
      double precision parm(nopmax)
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      double precision gridx(nmxgridx),gridq(nmxgridq)
      integer ngridx,ngridq,jx,jq
      CHARACTER*80 LINE
      character*512 setpath
      dimension pdf(-6:6)
      integer init,set,i,j,k,l,nset,iset
      parameter(nhess=0)
      double precision fgrid(0:nhess,161,161,0:7),grid(322)
      double precision agrid(161),alfas,Qalfa
      common/fgridz/fgrid,grid
      double precision up,dn,usea,dsea,str,chm,bot,glu
      double precision heragrid
      double precision qq(5),yntmp(5)
      real*8 mc,mc2,mb,mb2,mt,mt2,mz,mz2,alfa0,scale0
      COMMON/QCCONS/                                                    &
     &PI,PROTON,EUTRON,UCLEON,UDSCBT(6),AAM2H,BBM2H,AAM2L,BBM2L,        &
     &AAAR2,BBBR2,FL_FAC,CBMSTF(4:7),CHARGE(4:7),                       &
     &C1S3,C2S3,C4S3,C5S3,C8S3,C11S3,C14S3,C16S3,C20S3,C22S3,C28S3,     &
     &C38S3,C40S3,C44S3,C52S3,C136S3,C11S6,C2S9,C4S9,C10S9,C14S9,C16S9, &
     &C40S9,C44S9,C62S9,C112S9,C182S9,C11S12,C35S18,C61S12,C215S1,      &
     &C29S12,CPI2S3,CPIA,CPIB,CPIC,CPID,CPIE,CPIF,CCA,CCF,CTF,CATF,CFTF 
      save 
      x=xin
      q2=qin*qin 
      call getnset(iset)
      up =   HERAGRID(x,Q2,1)
      dn =   HERAGRID(x,Q2,2)
      usea = HERAGRID(x,Q2,3)
      dsea = HERAGRID(x,Q2,4)
      if(name(iset)(1:10).eq.'HERAGRID10') then
        up = up + usea
        dn = dn + dsea  
      endif
      str =  HERAGRID(x,Q2,5)
      chm =  HERAGRID(x,Q2,6)
      bot =  HERAGRID(x,Q2,7)
      glu =  HERAGRID(x,Q2,0)
      pdf(-6) = 0.0d0
       pdf(6) = 0.0d0
      pdf(-5) = bot
       pdf(5) = bot
      pdf(-4) = chm
       pdf(4) = chm
      pdf(-3) = str
       pdf(3) = str
      pdf(-2) = usea
       pdf(2) = up
      pdf(-1) = dsea
       pdf(1) = dn
       pdf(0) = glu
      return
!      
      entry HERAGRIDgetgrid(nset,ngridx,ngridq,gridx,gridq)

      ngridx=161
      do jx=1,ngridx
          gridx(jx)=grid(jx)
      enddo
      ngridq=161
      do jq=1,ngridq
          gridq(jq)=grid(jq+161)
      enddo
      return        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry HERAGRIDread(nset)
!     
      call getnmem(nset,imem) 

      read(1,*)nmem(nset),ndef(nset)
      
      do iq=1,23
         read(1,*) (grid(161+(iq-1)*7+ii),ii=1,7)
      enddo

      do jx=1,23
         read(1,*) (grid((jx-1)*7+ii),ii=1,7)
      enddo

    !read in to alphas grid  
      do iq=1,23
         read(1,*) (agrid((iq-1)*7+ii),ii=1,7)
      enddo

 
 ! - dummy read in to get to End: (stream 1 is still open) 

      do ns=0,nmem(nset)
        do k = 0,7
          do IQ=1,161
            do JX=1,23
	      read(1,'(a)')line
            enddo
          enddo
        enddo
      enddo
      
      return
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry HERAGRIDalfa(alfas,Qalfa)

       call getnset(iset)
       q2 = Qalfa*Qalfa
       nf=6
       if (Q2.lt.mt2) nf=5
       if (Q2.lt.mb2) nf=4
       if (Q2.lt.mc2) nf=3
       call GetOrderAsM(iset,iord)
       iord=iord+1
       call listPDF(iset,imem,parm)
       alfa0=parm(1)
       alfas=A0TOA1(q2,mz2,alfa0,iord,nf,ierr)      

      return
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry HERAGRIDinit(Eorder,Q2fit)

      call getnset(iset)
      mz = 91.187d0
      mz2=mz*mz
      call getQmassM(iset,4,mc)
      call getQmassM(iset,5,mb)
      call getQmassM(iset,6,mt)
      mc2=mc*mc
      mb2=mb*mb
      mt2=mt*mt
      UDSCBT(1) = 0.005 
      UDSCBT(2) = 0.01 
      UDSCBT(3) = 0.3 
      UDSCBT(4) = mc 
      UDSCBT(5) = mb 
      UDSCBT(6) = mt 

      return
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry HERAGRIDpdf(mem)
      imem = mem
      call getnset(iset)
      call setnmem(iset,imem)

  ! have to reopen stream 1 
  call getsetpath(setpath)
  open(1,file=setpath(1:len_trim(setpath)),action='READ')

  line = ''
  do while (line(2:11).ne.'Evolution:') 
     read(1,'(a)'),line
  enddo
  read(1,'(a)'),line
  read(1,'(a)'),line

  read(1,*)nmem(iset),ndef(iset)

  do iq=1,23
     read(1,*) (grid(161+(iq-1)*7+ii),ii=1,7)
  enddo

  do jx=1,23
     read(1,*) (grid((jx-1)*7+ii),ii=1,7)
  enddo

  !read in to alphas grid  
  do iq=1,23
     read(1,*) (agrid((iq-1)*7+ii),ii=1,7)
  enddo
 ! - dummy read in to get to End: (stream 1 is still open) 

   do ns=0,mem-1
     do k = 0,7
       do IQ=1,161
   	 do JX=1,23
     	   read(1,'(a)')line
   	 enddo
       enddo
     enddo
   enddo

! - read in the requested member 
   do k = 0,7
     do IQ=1,161
       do JX=1,23
   	 read(1,*)(fgrid(0,(jx-1)*7+ii,iq,k),ii=1,7)
       enddo
     enddo
   enddo
  
   close(1)
   return
!
 1000 format(5e13.5)
      end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! attempt at polynomial (order 4) interpolation based on  polint
      double precision function HERAGRID(x,Q2,n)
      implicit real*8(a-h,o-z)
      integer iset,imem,nhess
      parameter(nhess=0)
      double precision grid(322),x,Q2
      double precision fgrid(0:nhess,161,161,0:7)
      double precision ya(5,5),yntmp(5),ymtmp(5)
      common/fgridz/fgrid,grid
      
      npt = 5
      call getnset(iset)
      call getnmem(iset,imem)
      
! find the x bins around x
      nxlow = -1
      nxhi = 162
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
      nqhi = 162
      do while (nqhi-nqlow.gt.1)
         nqmid = (nqlow+nqhi)/2
	 if(q2.ge.grid(161+nqmid)) then
	    nqlow = nqmid
	 else 
	    nqhi = nqmid
	 endif
      enddo
      
! fill the temp 4x4 funtion array
	if(nxlow.le.0) nxlow=1
	if(nxlow.ge.161) nxlow=160
        if(nxlow.eq.1) then
	   nxbot = 1
	elseif(nxlow.eq.2) then
	   nxbot = 2
	else
	   if(nxlow.eq.160) then
	      nxbot = 4
	   else
	      nxbot = 3
	   endif
	endif
	if(nqlow.le.0) nqlow=1
	if(nqlow.ge.161) nqlow=160
        if(nqlow.eq.1) then
	   nqbot = 1
	elseif(nqlow.eq.2) then
	   nqbot = 2
	else
	   if(nqlow.eq.160) then
	      nqbot = 4
	   else   
	      nqbot = 3
	   endif
	endif
        do nx=1,5
          do nq=1,5 
	    ya(nx,nq) = fgrid(0,nxlow+nx-nxbot,nqlow+nq-nqbot,n)
	  enddo   
        enddo      
      
      do j=1,5
        do k=1,5
	  yntmp(k)=ya(j,k)
	enddo
	call herapolint(grid(161+nqlow-nqbot+1),yntmp,npt,q2,ymtmp(j),dy)
      enddo 
      call herapolint(grid(nxlow-nxbot+1),ymtmp,npt,x,y,dy)
      
      heragrid=y
      
      return
      end     
!=========================================================
         SUBROUTINE HERAPOLINT (XA,YA,N,X,Y,DY)

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
