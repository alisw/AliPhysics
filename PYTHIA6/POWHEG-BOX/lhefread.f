c...lhefheader(nlf)
c...reads initialization information from a les houches events file on unit nlf. 
      subroutine lhefreadhdr(nlf)
      implicit none
      integer nlf
      character * 100 string
      integer ipr
      include 'LesHouches.h'
 1    read(nlf,fmt='(a)',err=998,end=998) string
      if(string(1:5).eq.'<init') then
         read(nlf,*) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         do ipr=1,nprup
            read(nlf,*) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &           lprup(ipr)
         enddo
         goto 999
      else
         goto 1
      endif
 998  write(*,*) 'lhefreadhdr: could not find <init> data'
      call exit(1)
 999  end


c...reads event information from a les houches events file on unit nlf. 
      subroutine lhefreadev(nlf)
      implicit none
      integer nlf
      character * 100 string
      include 'LesHouches.h'
      integer i,j
 1    continue
      string=' '
      read(nlf,fmt='(a)',err=777,end=666) string
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      if(string(1:6).eq.'<event') then
c on error try next event. The error may be cause by merging
c truncated event files. On EOF return with no event found
         read(nlf,*,end=998,err=1)nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            read(nlf,*,end=998,err=1) idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         goto 999
      else
         goto 1
      endif
c no event found:
 777  continue
      print *,"Error in reading"
      print *,string
      stop
 666  continue
      print *,"reached EOF"
      print *,string
      stop
 998  continue
      print *,"read </LesHouchesEvents>"
      nup=0      
 999  end
