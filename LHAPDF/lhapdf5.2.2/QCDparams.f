      subroutine GetLam4(mem,lam4)
      implicit none
      integer mem,nset
      real*8 lam4,lam5
      nset = 1
      call GetLam4M(nset,mem,lam4)
      return
      
      entry GetLam5(mem,lam5)
      nset = 1
      call GetLam5M(nset,mem,lam5)
      return

      end
      
      subroutine initQCDparams(nset)
      implicit real*8(a-h,o-z)
      include 'parmsetup.inc'
      real*8 parmQCD(nmxset,0:noemax,2),lam4,lam5
      integer nset
c      integer iset,imem
c      common/SET/iset,imem
      save
      read(1,*)nmem,nrep
      if(nrep.eq.0) then
        do i=0,nmem
          read(1,*) (parmQCD(nset,i,j),j=1,2)
        enddo
      else
        read(1,*) (parmQCD(nset,0,j),j=1,2)
	do i=1,nmem
	  do j=1,2
	    parmQCD(nset,i,j)=parmQCD(nset,0,j)
	  enddo
        enddo
      endif
      return

      entry GetLam4M(nset,mem,lam4)
      lam4 = parmQCD(nset,mem,1)
      return
c
      entry GetLam5M(nset,mem,lam5)
      lam5 = parmQCD(nset,mem,2)
      return
c        
      end
      
      
      
