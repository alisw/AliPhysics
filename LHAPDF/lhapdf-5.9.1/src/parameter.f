! -*- F90 -*-


subroutine numberPDF(noe)
  implicit none
  integer nset,noe
  nset = 1
  call numberPDFM(nset,noe)
  return
end subroutine numberPDF

      
subroutine listPDF(nset, imem, parm)
  implicit none
  include 'parmsetup.inc'
  character*16 name(nmxset)
  integer nmem(nmxset),ndef(nmxset)
  common/NAME/name,nmem,ndef,mem
  character*16 s1
  integer i,j,mem,imem,noe,nop,listN(nmxset),listP(nmxset),type
  double precision parmL(nmxset,0:noemax,nopmax),parm(nopmax)
  integer nset
  save listN,listP,parmL
#ifdef NNPDF
! nnpdf variables
  integer MXPDF
  parameter(MXPDF=13)
  integer NTOTPDF
  parameter(NTOTPDF=5)      
  integer MXPAR
  parameter(MXPAR=2e2)
  integer MXL,MXN	   
  parameter(MXL=5,MXN=10)
  integer IPAR,IPDF,IDUM
  integer NPAR(MXPDF)
  integer MXREP
  parameter(MXREP=1e3)
  real*8 PARTR(MXPAR,0:MXREP,MXPDF)
  real*8 ANORMTR(MXPDF,0:MXREP)
  common/nnpdf10CPARTR/PARTR,ANORMTR
!
#endif
  mem=imem
  if (mem.gt.listN(nset)) then
     write(*,*) 'Maximum number of PDFs in list exceeded: ', mem,' > ',listN
     write(*,*) 'Returning most likely PDF'
     mem=0
  endif
  if (mem.lt.0) then
     write(*,*) 'Negative PDF member requested: ', mem
     write(*,*) 'Returning most likely PDF'
     mem=0
  endif
  do i=1,listP(nset)
     parm(i)=parmL(nset,mem,i)
  enddo
  return

  entry nopPDF(nset,nop)
  nop=listP(nset)
  return
  
  entry numberPDFM(nset, noe)
  if (      NAME(nset).eq.'MRSTgrid'    .or. NAME(nset).eq.'MRST98grid'  &
       .or. NAME(nset).eq.'A02'         .or. NAME(nset).eq.'A02M'        &
       .or. NAME(nset).eq.'ABKM09'      .or. NAME(nset).eq.'ABM11'       &
       .or. NAME(nset).eq.'CTEQ5grid'   .or. NAME(nset).eq.'CTEQ6grid'   &
       .or. NAME(nset).eq.'CTEQ65grid'  .or. NAME(nset).eq.'CTEQ66grid'  &
       .or. NAME(nset).eq.'CTEQ65cgrid' .or. NAME(nset).eq.'CTEQ6ABgrid' &
       .or. NAME(nset).eq.'SASG'        .or. NAME(nset).eq.'GRVG'        &
       .or. NAME(nset).eq.'DOG'         .or. NAME(nset).eq.'DGG'         &
       .or. NAME(nset).eq.'LACG'        .or. NAME(nset).eq.'GSG'         &
       .or. NAME(nset).eq.'GSG96'       .or. NAME(nset).eq.'ACFGP'       &
       .or. NAME(nset).eq.'WHITG'       .or. NAME(nset).eq.'OWP'         &
       .or. NAME(nset).eq.'SMRSP'       .or. NAME(nset).eq.'GRVP'        &
       .or. NAME(nset).eq.'ABFKWP'      .or. NAME(nset).eq.'ZEUSGRID'    &
       .or. NAME(nset)(1:5).eq.'GJR08'  .or. NAME(nset).eq.'NNPDFint'    &
       .or. NAME(nset).eq.'HKNgrid'                                      &
       .or. NAME(nset).eq.'CT12grid' & 
       !.or. NAME.eq.'H12000' &
       !.or. NAME.eq.'GRV' &
       ) then
     noe = nmem(nset)
  else
     noe=listN(nset)
  endif
  return
  
  entry InitListPDF(nset)
  type=-1
  read(1,*) s1,listN(nset),listP(nset)
  !print *,s1,listN(nset),listP(nset)
  if (index(s1,'list').eq.1) then
     type=1
     do i=0,listN(nset)
        read(1,*) (parmL(nset,i,j),j=1,listP(nset))
       !print *,i,(parmL(nset,i,j),j=1,listP(nset))
     enddo
      !print *,parmL(nset,0,1)
!
#ifdef NNPDF
  elseif (index(s1,'nnpar').eq.1) then
     type=2 ! is type=2 ok?

     do i=0,listN(nset)
  	read(1,*)idum
  	if(idum.ne.i)write(*,*)"error in InitListPDF"
  	do ipdf = 1,listP(nset)
  	   read(1,*) npar(ipdf)
  	   read(1,*) anormtr(ipdf,i)
  	   do ipar = 1,npar(ipdf)
  	      read(1,*) partr(ipar,i,ipdf)
  	   enddo
  	enddo
     enddo
#endif
  endif
  if (type.lt.0) then
     write(*,*) 'File description error:'
     write(*,*) 'Unknown parameter list type ',s1
     stop
  endif
  return
end subroutine listPDF
