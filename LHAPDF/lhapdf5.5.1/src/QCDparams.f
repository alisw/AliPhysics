! -*- F90 -*-


subroutine GetLam4(mem,lam4)
  implicit none
  integer mem,nset
  double precision lam4,lam5
  nset = 1
  call GetLam4M(nset,mem,lam4)
  return
  
  entry GetLam5(mem,lam5)
  nset = 1
  call GetLam5M(nset,mem,lam5)
  return
end subroutine GetLam4


subroutine GetXmin(mem,xmin)
  implicit none
  integer mem,nset
  double precision xmin,xmax,q2min,q2max
  nset = 1
  call GetXminM(nset,mem,xmin)
  return
  
  entry GetXmax(mem,xmax)
  nset = 1
  call GetXmaxM(nset,mem,xmax)
  return
  
  entry GetQ2min(mem,q2min)
  nset = 1
  call GetQ2minM(nset,mem,q2min)
  return
  
  entry GetQ2max(mem,q2max)
  nset = 1
  call GetQ2maxM(nset,mem,q2max)
  return
  
  entry GetMinMax(mem,xmin,xmax,q2min,q2max)
  nset = 1
  call GetMinMaxM(nset,mem,xmin,xmax,q2min,q2max)
  return
end subroutine GetXmin

      
subroutine initQCDparams(nset)
  implicit double precision(a-h,o-z)
  include 'parmsetup.inc'
  double precision parmQCD(nmxset,0:noemax,2),lam4,lam5
  integer nset
  ! integer iset,imem
  ! common/SET/iset,imem
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
  
  entry GetLam5M(nset,mem,lam5)
  lam5 = parmQCD(nset,mem,2)
  return
end subroutine initQCDparams

      
subroutine initMinMax(nset)
  implicit double precision(a-h,o-z)
  include 'parmsetup.inc'
  double precision parmXmin(nmxset,0:noemax),xmin
  double precision parmXmax(nmxset,0:noemax),xmax
  double precision parmQ2min(nmxset,0:noemax),q2min
  double precision parmQ2max(nmxset,0:noemax),q2max
  integer nset
  save
  read(1,*)nmem,nrep
  if(nrep.eq.0) then
     do i=0,nmem
        read(1,*) &
             parmXmin(nset,i),&
             parmXmax(nset,i),&
             parmQ2min(nset,i),&
             parmQ2max(nset,i)
     enddo
  else
     read(1,*) &
          parmXmin(nset,0),&
          parmXmax(nset,0),&
          parmQ2min(nset,0),&
          parmQ2max(nset,0)
     do i=1,nmem
        parmXmin(nset,i) = parmXmin(nset,0)
        parmXmax(nset,i) = parmXmax(nset,0)
        parmQ2min(nset,i) = parmQ2min(nset,0)
        parmQ2max(nset,i) = parmQ2max(nset,0)
     enddo
  endif
  return

  entry GetXminM(nset,mem,xmin)
  xmin = parmXmin(nset,mem)
  return
  
  entry GetXmaxM(nset,mem,xmax)
  xmax = parmXmax(nset,mem)
  return
  
  entry GetQ2minM(nset,mem,q2min)
  q2min = parmQ2min(nset,mem)
  return
  
  entry GetQ2maxM(nset,mem,q2max)
  q2max = parmQ2max(nset,mem)
  return
  
  entry GetMinMaxM(nset,mem,xmin,xmax,q2min,q2max)
  xmin = parmXmin(nset,mem)
  xmax = parmXmax(nset,mem)
  q2min = parmQ2min(nset,mem)
  q2max = parmQ2max(nset,mem)
  return
end subroutine initMinMax
