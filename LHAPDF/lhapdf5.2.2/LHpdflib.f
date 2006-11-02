c     Automatically determine the path to the system's PDF set
c     collection using the lhapdf-config utility (which must be
c     in the user's execution path
c     ---------------------------------------------------------
      subroutine InitPDFsetByCodes(code1, code2, code3)
      write(*,*) "Not implemented yet: this will move the 'glue' interf
     +ace to LHAPDF proper and use the InitPDFsetByName function to get
     +the path automatically."
      return
      end
c     ---------------------------------------------------------


c     Automatically determine the path to the system's PDF set
c     collection using the lhapdf-config utility (which must be
c     in the user's execution path
c     --------------------------------------------------------
      subroutine InitPDFsetByName(setname)
      implicit none
      character setname*(*)
      integer nset
      nset = 1
      call InitPDFsetByNameM(nset,setname)
      return
      end
      
      subroutine InitPDFsetByNameM(nset,setname)
      implicit none
      INTEGER LNBLNK
      character setname*(*)
      integer nset
c      integer :: ierror
      integer n, dirpathlength, setnamelength
      character*512 dirpath, setpath

      INTEGER LNROOT
      CHARACTER*1000 CHROOT
      CHROOT=' '

c check enviromental variable LHAPATH
      call getenv('LHAPATH',dirpath)
      if (dirpath.eq.'') then
C     Take the data from $ALICE_ROOT/LHAPDF/PDFsets
         CALL GETENV('ALICE_ROOT',CHROOT)
         LNROOT = LNBLNK(CHROOT)
         IF(LNROOT.LE.0) THEN
            dirpath='PDFsets'   ! Default value
         ELSE
            dirpath=CHROOT(1:LNROOT)//'/LHAPDF/PDFsets'
         ENDIF
      endif

c     Now do some mangling to get the right path length from the 
c     (hopefully) over-long string read in from the file
      n = 512
      do while (dirpath(n:n) .eq. ' ' .and. n .gt. 0)
         n = n-1
      enddo
      dirpathlength = n

c     How long is 'name', really?
      n = len(setname)
      do while (setname(n:n) .eq. ' ' .and. n .gt. 0)
         n = n-1
      enddo
      setnamelength = n

c     Combine the set directory path and the set name
      setpath(1:dirpathlength) = dirpath(1:dirpathlength)
      setpath(dirpathlength+1:dirpathlength+1) = "/"
      setpath(dirpathlength+2:dirpathlength+setnamelength+1) = setname(1
     $:setnamelength)
c      setpath(dirpathlength+setnamelength+2:dirpathlength+setnamelength+
c     $2) = ":"
c      write(*,*) setpath(1:dirpathlength+setnamelength+2)

      call InitPDFsetM(nset,setpath(1:dirpathlength+setnamelength+1))
      return
      end
c     ---------------------------------------------------------

      subroutine InitPDFset(setpath)
      implicit none
      integer nset
      character setpath*(*)
      nset = 1
      call InitPDFsetM(nset,setpath)
      return
      end      
c
      subroutine InitPDFsetM(nset,setpath)
      implicit none
      include 'parmsetup.inc'
      character setpath*(*)
      character*64 string
      character*16 s1,s2
      integer id,token,Ctoken
      integer lhasilent
      common/lhasilent/lhasilent
      integer nset,imem
c
      call setnset(nset)
c      
      open(unit=1,file=setpath,status='old')
      read(1,*) s1,s2
      if ((index(s2,'1.0').ne.1)
     +.and.(index(s2,'1.1').ne.1)
     +.and.(index(s2,'2.0').ne.1)
     +.and.(index(s2,'2.1').ne.1)
     +.and.(index(s2,'3.0').ne.1) 
     +.and.(index(s2,'3.1').ne.1)
     +.and.(index(s2,'4.0').ne.1)
     +.and.(index(s2,'5.0').ne.1))then
         write(*,*) 
     .        'Version ',s2,' not supported by this version of LHAPDF'
         stop
      else  
       if(lhasilent.eq.0) then
         write(*,*) '*************************************'
         write(*,*) '*       LHAPDF Version 5.2.2          *'
         write(*,*) '*************************************'
         write(*,*)
       endif
      endif
      id=Ctoken()
 1    read(1,*) string
      id=token(string)
c      print *,'id = ',id,string
      if (id.eq.0) then
         write(*,*) 'File description error:'
         write(*,*) 'Command not understood: ',string
         stop
      endif
      if (id.eq.1) call descriptionPDF(nset,id)
c      print *,'1/2'
      if (id.eq.2) call initEvolve(nset)
c      print *,'2/3'
      if (id.eq.3) call initAlphasPDF(nset)
c      print *,'3/4'
      if (id.eq.4) call initInputPDF(nset)
c      print *,'4/5'
      if (id.eq.5) call initListPDF(nset)
c      print *,'5/6'
      if (id.eq.6) call initQCDparams(nset)
c      print *,'6/7'
      if (id.ne.7) goto 1
      close(1)
c      print *,'calling InitEvolveCode',nset
      call InitEvolveCode(nset)
*
      return
      end
*     
      integer function token(s)
      implicit none
      character*16 s
      integer not,i,Ctoken
      parameter(not=7)
      character*16 t(not)
      data t/'Description:','Evolution:','Alphas:',
     .                    'Parametrization:','Parameterlist:',
     .                    'QCDparams:',
     .                    'End:'/
      integer count(not)
      save count
*
      token=0
      do i=1,not
         if (s.eq.t(i)) token=i
      enddo
      if (token.ne.0) then
         count(token)=count(token)+1
         if (count(token).eq.2) then
            write(*,*) 'File description error:'
            write(*,*) 'Second definition of entry: ',s
            stop
         endif
      endif
      return
*
      entry Ctoken()
      do i=1,not
         count(i)=0
      enddo
      Ctoken=0
      return
*     
      end
c
      subroutine LHAprint(iprint)
      implicit none
      integer lhasilent,iprint
      common/lhasilent/lhasilent
      lhasilent=iprint
c      print *,'lhasilent',lhasilent
      return
      end
