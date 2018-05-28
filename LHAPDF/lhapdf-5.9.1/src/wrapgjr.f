      subroutine GJRevolve(xin,qin,pdf)
      implicit real*8 (a-h,o-z)
      include 'parmsetup.inc'
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      double precision gridx(nmxgridx),gridq(nmxgridq)
      integer ngridx,ngridq,jx,jq
      CHARACTER*80 LINE
      dimension pdf(-6:6)
      integer ng(9),init,set,i,j,k,l,nset,iset
      double precision fgrid(118,99,-5:3,0:26),grid(217)
      !double precision fgrid(118,99,-5:3,0:0),grid(217)
!      common/fgridc/fgrid
      double precision upv,dnv,usea,dsea,str,chm,bot,glu
      double precision arg(9)
      double precision lha_dfint
      double precision lha_gjr08
      data ng /118,99,0,0,0,0,0,0,0/

      data grid &
     & /1d-9,1.25d-9,1.6d-9,2d-9,2.5d-9,3.16d-9,4d-9,5d-9,6.3d-9,8d-9, &
     &  1d-8,1.25d-8,1.6d-8,2d-8,2.5d-8,3.16d-8,4d-8,5d-8,6.3d-8,8d-8, &
     &  1d-7,1.25d-7,1.6d-7,2d-7,2.5d-7,3.16d-7,4d-7,5d-7,6.3d-7,8d-7, &
     &  1d-6,1.25d-6,1.6d-6,2d-6,2.5d-6,3.16d-6,4d-6,5d-6,6.3d-6,8d-6, &
     &  1d-5,1.25d-5,1.6d-5,2d-5,2.5d-5,3.16d-5,4d-5,5d-5,6.3d-5,8d-5, &
     &  1d-4,1.25d-4,1.6d-4,2d-4,2.5d-4,3.16d-4,4d-4,5d-4,6.3d-4,8d-4, &
     &  1d-3,1.25d-3,1.6d-3,2d-3,2.5d-3,3.16d-3,4d-3,5d-3,6.3d-3,8d-3, &
     &  1d-2,1.25d-2,1.6d-2,2d-2,2.5d-2,3.16d-2,4d-2,5d-2,6.3d-2,8d-2, &
     &  0.10d0,0.125d0,0.15d0,0.175d0,0.20d0,0.225d0,0.25d0,0.275d0, &
     &  0.30d0,0.325d0,0.35d0,0.375d0,0.40d0,0.425d0,0.45d0,0.475d0, &
     &  0.50d0,0.525d0,0.55d0,0.575d0,0.60d0,0.625d0,0.65d0,0.675d0, &
     &  0.70d0,0.725d0,0.75d0,0.775d0,0.80d0,0.825d0,0.85d0,0.875d0, &
     &  0.9d0,0.920d0,0.94d0,0.960d0,0.98d0,1d0, &
     &  0.3d0,0.31d0,0.35d0,0.375d0,0.4d0,0.45d0,0.5d0,0.51d0,0.525d0, &
     &  0.55d0,0.575d0,0.6d0,0.65d0,0.7d0,0.75d0,0.8d0,0.85d0,0.9d0, &
     &  1d0,1.25d0,1.6d0,2d0,2.5d0,3.16d0,4d0,5d0,6.3d0,8d0, &
     &  1d1,1.25d1,1.6d1,2d1,2.5d1,3.16d1,4d1,5d1,6.3d1,8d1, &
     &  1d2,1.25d2,1.6d2,2d2,2.5d2,3.16d2,4d2,5d2,6.3d2,8d2, &
     &  1d3,1.25d3,1.6d3,2d3,2.5d3,3.16d3,4d3,5d3,6.3d3,8d3, &
     &  1d4,1.25d4,1.6d4,2d4,2.5d4,3.16d4,4d4,5d4,6.3d4,8d4, &
     &  1d5,1.25d5,1.6d5,2d5,2.5d5,3.16d5,4d5,5d5,6.3d5,8d5, &
     &  1d6,1.25d6,1.6d6,2d6,2.5d6,3.16d6,4d6,5d6,6.3d6,8d6, &
     &  1d7,1.25d7,1.6d7,2d7,2.5d7,3.16d7,4d7,5d7,6.3d7,8d7,1d8/

      save 
      x=xin
      q2=qin*qin 
      call getnset(iset)
      call getnmem(iset,imem)
       upv =  LHA_GJR08(x,Q2,grid,fgrid,ng,1,imem)
      dnv =  LHA_GJR08(x,Q2,grid,fgrid,ng,2,imem)
      usea = LHA_GJR08(x,Q2,grid,fgrid,ng,-1,imem)
      dsea = LHA_GJR08(x,Q2,grid,fgrid,ng,-2,imem)
      str =  LHA_GJR08(x,Q2,grid,fgrid,ng,-3,imem)
      glu =  LHA_GJR08(x,Q2,grid,fgrid,ng,0,imem)
      pdf(-6) = 0.0d0
       pdf(6) = 0.0d0
      pdf(-5) = 0.0d0
       pdf(5) = 0.0d0
      pdf(-4) = 0.0d0
       pdf(4) = 0.0d0
      pdf(-3) = str
       pdf(3) = str
      pdf(-2) = usea
       pdf(2) = upv+usea
      pdf(-1) = dsea
       pdf(1) = dnv+dsea
       pdf(0) = glu
      if(name(iset)(1:7).eq.'GJR08VF'.or. &
     &   name(iset)(1:7).eq.'GJR08LO') then
        chm =  LHA_GJR08(x,Q2,grid,fgrid,ng,-4,imem)
       bot =  LHA_GJR08(x,Q2,grid,fgrid,ng,-5,imem)
       pdf(-5) = bot
        pdf(5) = bot
       pdf(-4) = chm
        pdf(4) = chm
      endif
      return
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      
      entry GJRgetgrid(nset,ngridx,ngridq,gridx,gridq)
     
      ngridx=118
      do jx=1,118
          gridx(jx)=grid(jx)
      enddo

      ngridq=99      
      do jq=1,99
          gridq(jq)=grid(118+jq)
      enddo 
       
      return
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GJRread(nset)
!     
      call getnmem(nset,imem) 
      read(1,*)nmem(nset),ndef(nset)
      lstart = -3
      if(name(nset)(1:7).eq.'GJR08VF') lstart=-5
      do i=0,nmem(nset)
      !do ii=0,nmem(nset)
!       i = 0
       do j=1,118
        do k=1,99
         if(name(nset)(1:7).eq.'GJR08VF'.or. &
     &      name(nset)(1:7).eq.'GJR08LO') then
           read(1,*) fgrid(j,k,-5,i),fgrid(j,k,-4,i), &
     &  	   fgrid(j,k,-3,i),fgrid(j,k,-2,i),fgrid(j,k,-1,i), &
     &  	   fgrid(j,k,0,i),fgrid(j,k,1,i),fgrid(j,k,2,i), &
     &  	   fgrid(j,k,3,i)
         else
          read(1,*) fgrid(j,k,-3,i),fgrid(j,k,-2,i),fgrid(j,k,-1,i), &
     &  	   fgrid(j,k,0,i),fgrid(j,k,1,i),fgrid(j,k,2,i), &
     &  	    fgrid(j,k,3,i)
         endif
        if(name(nset)(1:7).ne.'GJR08LO') then
          do  l=-lstart,3
            if (grid(118+k) < 0.5d0) then
              fgrid(j,k,l,i)=0d0
              fgrid(j,k,l,i)=0d0/fgrid(j,k,l,i)
            endif
          enddo
         endif
        enddo
       enddo
      enddo

      return
!

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GJRalfa(alfas,qalfa)
      call getnset(iset)
      call getnmem(iset,imem)
      arg(1) = 1d-9
      arg(2) = qalfa*qalfa
!      imem = 0
      alfas = lha_dfint(9,arg,ng,grid,fgrid(1,1,3,imem))
      return
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GJRinit(Eorder,Q2fit)
      return
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GJRpdf(mem)
      imem = mem
      call getnset(iset)
       call setnmem(iset,imem)
      return
!
 1000 format(5e13.5)
      end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function LHA_GJR08(x,Q2,grid,fgrid,ng,n,set)
      implicit none
      integer ng(9),n,set
      double precision grid(217),arg(9),x,Q2
      double precision lha_dfint
      double precision fgrid(118,99,-5:3,0:26)
!      common/fgridc/fgrid
      arg(1) = x
      arg(2) = Q2
       LHA_GJR08 = lha_dfint(9,arg,ng,grid,fgrid(1,1,n,set))
      end
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 
!! CERNLIB E104 modified to be used with GJR08 GRIDS:
!! Name changed from fint to dfint.
!! Name changed from dfint to lha_dfint.
!! Real variables changed to double precision.
!! External references to CERNLIB (error handling) routines removed.
          DOUBLE PRECISION FUNCTION LHA_DFINT(NARG,ARG,NENT,ENT,TABLE)
          INTEGER   NENT(9), INDEX(32)
          DOUBLE PRECISION ARG(9),   ENT(9),   TABLE(9), WEIGHT(32)
          LHA_DFINT  =  0d0
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1
          INDEX(1)  =  1
          WEIGHT(1) =  1d0
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0d0)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
!             IF(X-ENT(LOCC))  12, 20, 13
             IF(X-ENT(LOCC).lt.0)  goto 12
             IF(X-ENT(LOCC).eq.0)  goto 20
             IF(X-ENT(LOCC).gt.0)  goto 13            
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT
  22            CONTINUE
             GOTO 90
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             LHA_DFINT  =  LHA_DFINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
          RETURN
 300      WRITE(*,1000) NARG
          STOP
1000      FORMAT( 7X, 24HFUNCTION DFINT... NARG =,I6, &
     &              17H NOT WITHIN RANGE)
          END

