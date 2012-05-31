      subroutine descriptionPDF(nset,id)
      implicit none
      include 'parmsetup.inc'
      integer nset
      integer id,token,l,nline
      character*64 string
      character*64 desc(nmxset,linemax)
      integer lhasilent
      common/lhasilent/lhasilent
      save nline,desc
*
      l=0
      if(lhasilent.eq.0) 
     + write(*,*) '>>>>>> PDF description: <<<<<<'
 1    read(1,*) string
      id=token(string)
      if (id.eq.0) then
         if(lhasilent.eq.0) then
          write(*,*) string
         endif
         l=l+1
         if (l.gt.linemax) then
            write(*,*) 'Too many lines in PDF description to store.'
            write(*,*) 'Increase linemax variable in parmsetup.inc'
            write(*,*) 'Ignoring additional description lines'
            l=linemax
         else
            desc(nset,l)=string
            goto 1
         endif
      endif
      nline=l
      if(lhasilent.eq.0) then
       write(*,*)'>>>>>>                   <<<<<<'
       write(*,*)
      endif
      return
*      
      entry GetDescM(nset)
      do l=1,nline
         write(*,*) desc(nset,l)
      enddo
      return
*     
      end
*     
      subroutine GetDesc()
      implicit none
      integer nset
      nset = 1
      call GetDescM(nset)
      return
      end
