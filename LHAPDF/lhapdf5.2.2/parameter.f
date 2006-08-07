      subroutine numberPDF(noe)
      implicit none
      integer nset,noe
      nset = 1
      call numberPDFM(nset,noe)
      return
      end
c      
      subroutine listPDF(nset,imem,parm)
      implicit none
      include 'parmsetup.inc'
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset)
      common/NAME/name,nmem,ndef,mem
      character*16 s1
      integer i,j,mem,imem,noe,nop,listN(nmxset),listP(nmxset),type
      real*8 parmL(nmxset,0:noemax,nopmax),parm(nopmax)
      integer nset
      save listN,listP,parmL
*
       mem=imem
      if (mem.gt.listN(nset)) then
         write(*,*) 'Maximum number of PDFs in list exceeded: ',
     .              mem,' > ',listN
         write(*,*) 'Returning most likely PDF'
         mem=0
      endif
      if (mem.lt.0) then
         write(*,*) 'Negative PDF member requested: ',mem
         write(*,*) 'Returning most likely PDF'
         mem=0
      endif
      do i=1,listP(nset)
         parm(i)=parmL(nset,mem,i)
      enddo
      return
*
      entry nopPDF(nset,nop)
      nop=listP(nset)
      return
*
      entry numberPDFM(nset,noe)
      if(NAME(nset).eq.'MRSTgrid'
     +.or.NAME(nset).eq.'MRST98grid'
     +.or.NAME(nset).eq.'A02'
     +.or.NAME(nset).eq.'A02M'
     +.or.NAME(nset).eq.'CTEQ5grid'
     +.or.NAME(nset).eq.'CTEQ6grid'
     +.or.NAME(nset).eq.'CTEQ6ABgrid'		!*** added ***
     +.or.NAME(nset).eq.'SASG'
     +.or.NAME(nset).eq.'GRVG'
     +.or.NAME(nset).eq.'DOG'
     +.or.NAME(nset).eq.'DGG'
     +.or.NAME(nset).eq.'LACG'
     +.or.NAME(nset).eq.'GSG'
     +.or.NAME(nset).eq.'GSG96'
     +.or.NAME(nset).eq.'ACFGP'
     +.or.NAME(nset).eq.'WHITG'
     +.or.NAME(nset).eq.'OWP'
     +.or.NAME(nset).eq.'SMRSP'
     +.or.NAME(nset).eq.'GRVP'
     +.or.NAME(nset).eq.'ABFKWP'
c     +.or.NAME.eq.'H12000'
c     +.or.NAME.eq.'GRV'
     +) then
        noe = nmem(nset)
      else
        noe=listN(nset)
      endif
     
      return
*
      entry InitListPDF(nset)
      type=-1
      read(1,*) s1,listN(nset),listP(nset)
c      print *,s1,listN,listP
      if (index(s1,'list').eq.1) then
         type=1
         do i=0,listN(nset)
            read(1,*) (parmL(nset,i,j),j=1,listP(nset))
         enddo
c	 print *,parmL(0,1)
      endif
      if (type.lt.0) then
         write(*,*) 'File description error:'
         write(*,*) 'Unknown parameter list type ',s1
         stop
      endif
      return
*
      end
      
