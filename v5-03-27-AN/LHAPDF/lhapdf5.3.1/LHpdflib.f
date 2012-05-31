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
      include 'parmsetup.inc'
      include 'pathsetup.inc'
      integer LNBLNK
      common/LHAPDFC/lhapath
      character*20 lhaparm(20)
      real*8 lhavalue(20)
      common/LHACONTROL/lhaparm,lhavalue
      character setname*(*)
      integer nset
c      integer :: ierror
      integer n, dirpathlength, setnamelength
      character*512 dirpath, homepath, dotlhapath, cachepath, setpath
c pr      
      INTEGER LNROOT
      INTEGER LEN_TRIM
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

c     Now build the path to the PDF set
      setpath = dirpath(:LEN_TRIM(dirpath)) // "/" //
     +setname(:LEN_TRIM(setname))

c     Initialize using the detected PDF set
      call InitPDFsetM(nset, setpath(:LEN_TRIM(setpath)))
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
     +.and.(index(s2,'5.0').ne.1)
     +.and.(index(s2,'5.3').ne.1))then
         write(*,*) 
     .        'Version ',s2,' not supported by this version of LHAPDF'
         stop
      else  
       if(lhasilent.eq.0) then
         write(*,*) '*************************************'
         write(*,*) '*       LHAPDF Version 5.3.0          *'
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
      if (id.eq.7) call initMinMax(nset)
c      print *,'7/8'
      if (id.ne.8) goto 1
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
      parameter(not=8)
      character*16 t(not)
      data t/'Description:','Evolution:','Alphas:',
     .                    'Parametrization:','Parameterlist:',
     .                    'QCDparams:','MinMax:',
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
