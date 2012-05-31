! -*- F90 -*-


!   LHAGLUE Interface to LHAPDF library of modern parton
!    density functions (PDF) with uncertainties
! 
! Authors for v4: Dimitri Bourilkov, Craig Group, Mike Whalley
! 
! Authors for v3: Dimitri Bourilkov, Craig Group, Mike Whalley
! 
! Author for v1 and v2: Dimitri Bourilkov  bourilkov@mailaps.org
!                       University of Florida
! 
! HERWIG interface by Dimitri Bourilkov and Craig Group
! 
! New numbering scheme and upgrade for LHAPDF v2.1
! by Dimitri Bourilkov and Mike Whalley
! 
! For more information, or when you cite this interface, currently
! the official reference is:
! D.Bourilkov, "Study of Parton Density Function Uncertainties with
! LHAPDF and PYTHIA at LHC", hep-ph/0305126.
! 
! The official LHAPDF page is:
! 
!    http://durpdg.dur.ac.uk/lhapdf/index.html 
! 
! The interface contains four subroutines (similar to PDFLIB).
! It can be used seamlessly by Monte Carlo generators 
! interfaced to PDFLIB or in stand-alone mode.
! 
!     For initialization (called once)
! 
!     PDFSET(PARM,VALUE)
! 
!     For the proton/pion structure functions
! 
!     STRUCTM(X,Q,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
! 
!     For the photon structure functions
! 
!     STRUCTP(X,Q2,P2,IP2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
! 
!     For statistics ON structure functions (under/over-flows)
! 
!     PDFSTA
! 
! This interface can be invoked in 3 ways depending
! on the value of PARM(1) provided by the user when
! calling PDFSET(PARM,VALUE):
! 
!     For PYTHIA:         PARM(1).EQ.'NPTYPE'
!       (this is set automatically by PYTHIA)
! 
!     For HERWIG:         PARM(1).EQ.'HWLHAPDF'
!       (set by the USER e.g. in the main program like this:
!           AUTPDF(1) = 'HWLHAPDF'
!           AUTPDF(2) = 'HWLHAPDF'                         )
! 
!     For Stand-alone:    PARM(1).EQ.'DEFAULT'
!       (can be used for PDF studies or when interfacing
!        new generators)
! 
! The LHAPDF set/member is selected depending on the value of:
! 
!         PYTHIA:   ABS(MSTP(51)) - proton
!                   ABS(MSTP(53)) - pion
!                   ABS(MSTP(55)) - photon
! 
!         HERWIG:   ABS(INT(VALUE(1)))
! 
!    STAND-ALONE:   ABS(INT(VALUE(1)))
! 
! 
!         CONTROL switches:
!        ==================
! 
!      THE LOCATION OF THE LHAPDF LIBRARY HAS TO BE SPECIFIED
!      AS DESCRIBED BELOW (the rest is optional)
! 
!      if the user does nothing, sensible defaults
!      are active; to change the behaviour, the corresponding
!      values of LHAPARM() should be set to the values given below
! 
!   Location of the LHAPDF library of PDFs (pathname):
!      uses common block /LHAPDFC/
! 
!    If the user does nothing => default = subdir PDFsets of the 
!                               current directory (can be real subdir
!                               OR a soft link to the real location)
!    If the user sets LHAPATH => supplied by the USER who defines the
!                         path in common block COMMON/LHAPDFC/LHAPATH
!                         BEFORE calling PDFSET
! 
!    Other controls:
!    ===============
!      use common block /LHACONTROL/
! 
!   Collect statistics on under/over-flow requests for PDFs
!   outside their validity ranges in X and Q**2
!   (call PDFSTA at end of run to print it out)
! 
!       LHAPARM(16).EQ.'NOSTAT' => No statistics (faster)
!       LHAPARM(16).NE.'NOSTAT' => Default: collect statistics
! 
!   Option to use the values for the strong coupling alpha_s
!   as computed in LHAPDF in the MC generator
!   (to ensure uniformity between the MC generator and the PDF set)
!   WARNING: implemented ONLY for PYTHIA in LHAPDFv4
! 
!       LHAPARM(17).EQ.'LHAPDF' => Use alpha_s from LHAPDF
!       LHAPARM(17).NE.'LHAPDF' => Default (same as LHAPDF v1/v3)
! 
!   Extrapolation of PDFs outside LHAPDF validity range given by
!   [Xmin,Xmax] and [Q2min,Q2max]; DEFAULT => PDFs "frozen" at the
!   boundaries
! 
!       LHAPARM(18).EQ.'EXTRAPOLATE' => Extrapolate PDFs on OWN RISK
!                           WARNING: Crazy values can be returned
! 
!   Printout of initialization information in PDFSET (by default)
! 
!       LHAPARM(19).EQ.'SILENT' => No printout (silent mode)
!       LHAPARM(19).EQ.'LOWKEY' => Print 5 times (almost silent mode)
! 
!
!*********************************************************************
!
! $Id: lhaglue.f 356 2008-08-28 15:58:02Z buckley $
!
! $Log$
! Revision 1.7  2005/12/02 14:50:54  whalley
! Changes for new CTEQ code/AB sets
!
! Revision 1.6  2005/10/18 15:35:48  whalley
! fix to allow LHAPATH to be user defined as well as lhapdf-config
!
! Revision 1.5  2005/10/18 11:47:48  whalley
! Change to only set LHAPATH once per run
!
! Revision 1.1.1.2  1996/10/30 08:29:06  cernlib
! Version 7.04
!
! Revision 1.1.1.1  1996/04/12 15:29:26  plothow
! Version 7.01
!
! v5.0  06-Oct-2005  Major change to allow multiset-initializations 
! v4.0  28-Apr-2005  PDFSTA routine; option to use Alfa_s from LHAPDF
! v4.0  21-Mar-2005  Photon/pion/new p PDFs, updated for LHAPDF v4
! v3.1  26-Apr-2004  New numbering scheme, updated for LHAPDF v2/v3
! v3.0  23-Jan-2004  HERWIG interface added
! v2.0  20-Sep-2003  PDFLIB style adopted
! v1.0  05-Mar-2003  First working version from PYTHIA to LHAPDF v1
! 
! interface to LHAPDF library
!*********************************************************************


! PDFSET
! Initialization for use of parton distributions
!  according to the LHAPDF interface.
! 
! v4.0  28-Apr-2005  Option to use Alfa_s from LHAPDF
! v4.0  21-Mar-2005  Photon/pion/new p PDFs, updated for LHAPDF v4
! v3.1  26-Apr-2004  New numbering scheme
! v3.0  23-Jan-2004  HERWIG interface added
! 
! Interface to LHAPDF library
 subroutine pdfset(parm,value, mstu11, mstp51, mstp53, mstp55, qcdl4, qcdl5, axmin, axmax, aq2min, aq2max)
 entry pdfset_herwig(parm, value) 
  ! Double precision and integer declarations.
  implicit double precision(a-h, o-z)
  implicit integer(i-n)
  ! Common blocks
  include 'commonlhapdf.inc'
  include 'commonlhasets.inc'
  include 'commonlhapdfc.inc'
  include 'commonlhacontrol.inc'
  include 'commonlhaglsta.inc'
  ! additions for multiset use
  double precision xxmin(nmxset),xxmax(nmxset),qq2min(nmxset),qq2max(nmxset)
  save xxmin,xxmax,qq2min,qq2max
  ! Interface to LHAPDFLIB.
  double precision qcdlha4, qcdlha5
  integer nfllha
  common/lhapdfr/qcdlha4, qcdlha5, nfllha
  save /lhapdfr/
  integer lhaextrp
  common/lhapdfe/lhaextrp
  save /lhapdfe/
  integer lhasilent
  common/lhasilent/lhasilent
  save /lhasilent/
  ! Interface to PDFLIB.
  common/w50511/ nptypepdfl,ngrouppdfl,nsetpdfl,modepdfl,nflpdfl,lopdfl,tmaspdfl
  save /w50511/
  double precision tmaspdfl
  double precision qcdl4,qcdl5
  double precision axmin, axmax, aq2min, aq2max
  ! Interface to PDFLIB.
  common/w50513/xmin,xmax,q2min,q2max
  save /w50513/
  double precision xmin,xmax,q2min,q2max
  ! Local arrays and character variables (NOT USED here DB)
  character*20 parm(20)
  double precision value(20)
  integer lhapathlen
  integer :: lhainput = 1
  !integer lhaselect
  integer lhaprint
  integer lhaonce
  integer lhafive
  save lhaonce
  save lhafive
  data lhaonce/0/
  data lhafive/0/
  logical first
  data first/.true./
  character*512 dirpath
  save first
  character*1000 chroot
  chroot=' '
  ! Initialise common blocks  
  call commoninit()

!  if (first) then
!     call getdirpath(dirpath)
!     first = .FALSE.
!  endif
  if(first .AND. (LHAPARM(20).NE.'LHAPATH')) then
!...overide the default PDFsets path
! ... check first if the environmental variable LHAPATH is set
         call getenv('LHAPATH',lhapath)
         if(lhapath.eq.'') then
!     The environment variable LHAPATH is not set.
!     Take the data from $ALICE_ROOT/LHAPDF/PDFsets
            CALL GETENV('ALICE_ROOT',CHROOT)
            LNROOT = LNBLNK(CHROOT)
            IF(LNROOT.LE.0) THEN
               LHAPATH='PDFsets' ! Default value
            ELSE
               LHAPATH=CHROOT(1:LNROOT)//'/LHAPDF/PDFsets'
            ENDIF
         endif
      first=.FALSE.
      endif


  ! Init
  lhaextrp = 0
  if(lhaparm(18).EQ.'EXTRAPOLATE') then  ! Extrapolate PDFs on own risk
     lhaextrp = 1
  endif
  lhasilent = 0
  if (lhaparm(19).EQ.'SILENT') then    !  No printout (silent MODE)
     lhasilent = 1
  elseif (lhaparm(19).EQ.'LOWKEY') then ! print 5 times (lowkey mode)
     if (lhafive .lt. 6) then
        lhafive = lhafive + 1
     else
        lhasilent = 1
     endif
  endif
  if (parm(1).EQ.'NPTYPE') then        !  pythia
     lhaprint = mstu11
     if(value(1) .eq. 1) then         !   nucleon
        lhainput = abs(mstp51)
     elseif(value(1) .eq. 2) then     !   pion
        lhainput = abs(mstp53)
     elseif(value(1) .eq. 3) then     !   photon
        lhainput = abs(mstp55)
     endif
     if (lhasilent .ne. 1) print *,'==== PYTHIA WILL USE LHAPDF ===='
  elseif(parm(1).EQ.'HWLHAPDF') then  !  herwig
     lhainput = abs(int(value(1)))
     if(lhaonce.eq.lhainput) return
     if(lhasilent .ne. 1) print *,'==== HERWIG WILL USE LHAPDF ===='
     lhaprint = 6
     lhaonce = lhainput
  elseif(parm(1).EQ.'DEFAULT') then  !  stand-alone
     lhainput = abs(int(value(1)))
     if(lhaonce.eq.lhainput) return
     if(lhasilent .ne. 1) print *,'==== STAND-ALONE LHAGLUE MODE TO USE LHAPDF ===='
     lhaprint = 6
     lhaonce = lhainput
  else
     print *,'== UNKNOWN LHAPDF INTERFACE CALL! STOP EXECUTION! =='
     stop
  endif
  ! Initialize parton distributions: LHAPDFLIB.
  lhapathlen=index(lhapath,' ') - 1
  lhaset = lhainput
  xmin = 1.0D-6      ! X_min for current PDF set
  xmax = 1.0D0       ! X_max for current PDF set
  q2min = 1.0D0**2   ! Q**2_min scale for current PDF set [GeV]
  q2max = 1.0D5**2   ! Q**2_max scale for current PDF set [GeV]
  ! 
  ! Protons
  ! 
  ! CTEQ Family
  if ((lhainput .ge. 10000) .and. (lhainput .le. 19999)) then
     q2max = 1.0d08
     if ((lhainput .ge. 10000) .and. (lhainput .le. 10040)) then
        lhaset = 10000
        lhaname=lhapath(1:lhapathlen)//'/cteq6.LHpdf'
        q2min = 1.69d0
     elseif((lhainput .ge. 10041) .and. (lhainput .le. 10041)) then
        lhaset = 10041
        lhaname=lhapath(1:lhapathlen)//'/cteq6l.LHpdf'
        q2min = 1.69d0
     elseif((lhainput .ge. 10042) .and. (lhainput .le. 10042)) then
        lhaset = 10042
        lhaname=lhapath(1:lhapathlen)//'/cteq6ll.LHpdf'
        q2min = 1.69d0
     elseif((lhainput .ge. 10050) .and. (lhainput .le. 10090)) then
        lhaset = 10050
        lhaname=lhapath(1:lhapathlen)//'/cteq6mE.LHgrid'
        q2min = 1.69d0
     elseif((lhainput .ge. 10100) .and. (lhainput .le. 10140)) then
        lhaset = 10100
        lhaname=lhapath(1:lhapathlen)//'/cteq61.LHpdf'
        q2min = 1.69d0
     elseif((lhainput .ge. 10150) .and. (lhainput .le. 10190)) then
        lhaset = 10150
        lhaname=lhapath(1:lhapathlen)//'/cteq61.LHgrid'
        q2min = 1.69d0
     elseif((lhainput .ge. 10250) .and. (lhainput .le. 10269)) then
        lhaset = 10250
        lhaname=lhapath(1:lhapathlen)//'/cteq6AB.LHgrid'
        q2min = 1.69d0
     elseif((lhainput .ge. 10350) .and. (lhainput .le. 10390)) then
        lhaset = 10350
        lhaname=lhapath(1:lhapathlen)//'/cteq65.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-7
     elseif((lhainput .ge. 10450) .and. (lhainput .le. 10456)) then
        lhaset = 10450
        lhaname=lhapath(1:lhapathlen)//'/cteq65c.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-7
     elseif((lhainput .ge. 10460) .and. (lhainput .le. 10467)) then
        lhaset = 10460
        lhaname=lhapath(1:lhapathlen)//'/cteq65s.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-7
     elseif((lhainput .ge. 10550) .and. (lhainput .le. 10594)) then
        lhaset = 10550
        lhaname=lhapath(1:lhapathlen)//'/cteq66.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10650) .and. (lhainput .le. 10653)) then
        lhaset = 10650
        lhaname=lhapath(1:lhapathlen)//'/cteq66c.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10660) .and. (lhainput .le. 10663)) then
        lhaset = 10660
        lhaname=lhapath(1:lhapathlen)//'/cteq66a.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10670) .and. (lhainput .le. 10677)) then
        lhaset = 10670
        lhaname=lhapath(1:lhapathlen)//'/cteq6lg.LHgrid'
        q2min = 1.69d0
     elseif((lhainput .ge. 19050) .and. (lhainput .le. 19050)) then
        lhaset = 19050
        lhaname=lhapath(1:lhapathlen)//'/cteq5m.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19051) .and. (lhainput .le. 19051)) then
        lhaset = 19051
        lhaname=lhapath(1:lhapathlen)//'/cteq5m1.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19053) .and. (lhainput .le. 19053)) then
        lhaset = 19053
        lhaname=lhapath(1:lhapathlen)//'/cteq5f3.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19054) .and. (lhainput .le. 19054)) then
        lhaset = 19054
        lhaname=lhapath(1:lhapathlen)//'/cteq5f4.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19060) .and. (lhainput .le. 19060)) then
        lhaset = 19060
        lhaname=lhapath(1:lhapathlen)//'/cteq5d.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19070) .and. (lhainput .le. 19070)) then
        lhaset = 19070
        lhaname=lhapath(1:lhapathlen)//'/cteq5l.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19150) .and. (lhainput .le. 19150)) then
        lhaset = 19150
        lhaname=lhapath(1:lhapathlen)//'/cteq4m.LHgrid'
        q2min = 2.56d0
        xmin=1.0d-5
     elseif((lhainput .ge. 19160) .and. (lhainput .le. 19160)) then
        lhaset = 19160
        lhaname=lhapath(1:lhapathlen)//'/cteq4d.LHgrid'
        q2min = 2.56d0
        xmin=1.0d-5
     elseif((lhainput .ge. 19170) .and. (lhainput .le. 19170)) then
        lhaset = 19170
        lhaname=lhapath(1:lhapathlen)//'/cteq4l.LHgrid'
        q2min = 2.56d0
        xmin=1.0d-5
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! MRST Family
  elseif((lhainput .ge. 20000) .and. (lhainput .le. 29999)) then
     q2min = 1.25d0
     q2max = 1.0d07
     xmin = 1.0d-5
     if((lhainput .ge. 20000) .and. (lhainput .le. 20004)) then
        lhaset = 20000
        lhaname=lhapath(1:lhapathlen)//'/MRST2001nlo.LHpdf'
     elseif((lhainput .ge. 20050) .and. (lhainput .le. 20054)) then
        lhaset = 20050
        lhaname=lhapath(1:lhapathlen)//'/MRST2001nlo.LHgrid'
     elseif((lhainput .ge. 20060) .and. (lhainput .le. 20061)) then
        lhaset = 20060
        lhaname=lhapath(1:lhapathlen)//'/MRST2001lo.LHgrid'
     elseif((lhainput .ge. 20070) .and. (lhainput .le. 20074)) then
        lhaset = 20070
        lhaname=lhapath(1:lhapathlen)//'/MRST2001nnlo.LHgrid'
     elseif((lhainput .ge. 20100) .and. (lhainput .le. 20130)) then
        lhaset = 20100
        lhaname=lhapath(1:lhapathlen)//'/MRST2001E.LHpdf'
     elseif((lhainput .ge. 20150) .and. (lhainput .le. 20180)) then
        lhaset = 20150
        lhaname=lhapath(1:lhapathlen)//'/MRST2001E.LHgrid'
     elseif((lhainput .ge. 20200) .and. (lhainput .le. 20201)) then
        lhaset = 20200
        lhaname=lhapath(1:lhapathlen)//'/MRST2002nlo.LHpdf'
     elseif((lhainput .ge. 20250) .and. (lhainput .le. 20251)) then
        lhaset = 20250
        lhaname=lhapath(1:lhapathlen)//'/MRST2002nlo.LHgrid'
     elseif((lhainput .ge. 20270) .and. (lhainput .le. 20271)) then
        lhaset = 20270
        lhaname=lhapath(1:lhapathlen)//'/MRST2002nnlo.LHgrid'
     elseif((lhainput .ge. 20300) .and. (lhainput .le. 20301)) then
        lhaset = 20300
        lhaname=lhapath(1:lhapathlen)//'/MRST2003cnlo.LHpdf'
        q2min = 10.0d0
        xmin = 1.0d-3
     elseif((lhainput .ge. 20350) .and. (lhainput .le. 20351)) then
        lhaset = 20350
        lhaname=lhapath(1:lhapathlen)//'/MRST2003cnlo.LHgrid'
        q2min = 10.0d0
        xmin = 1.0d-3
     elseif((lhainput .ge. 20370) .and. (lhainput .le. 20371)) then
        lhaset = 20370
        lhaname=lhapath(1:lhapathlen)//'/MRST2003cnnlo.LHgrid'
        q2min = 7.0d0
        xmin = 1.0d-3
     elseif((lhainput .ge. 20400) .and. (lhainput .le. 20401)) then
        lhaset = 20400
        lhaname=lhapath(1:lhapathlen)//'/MRST2004nlo.LHpdf'
     elseif((lhainput .ge. 20406) .and. (lhainput .le. 20407)) then
        lhaset = 20406
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF3nlo.LHpdf'
     elseif((lhainput .ge. 20408) .and. (lhainput .le. 20409)) then
        lhaset = 20408
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF4nlo.LHpdf'
     elseif((lhainput .ge. 20450) .and. (lhainput .le. 20451)) then
        lhaset = 20450
        lhaname=lhapath(1:lhapathlen)//'/MRST2004nlo.LHgrid'
     elseif((lhainput .ge. 20452) .and. (lhainput .le. 20453)) then
        lhaset = 20452
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF3lo.LHgrid'
     elseif((lhainput .ge. 20454) .and. (lhainput .le. 20455)) then
        lhaset = 20454
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF4lo.LHgrid'
     elseif((lhainput .ge. 20456) .and. (lhainput .le. 20457)) then
        lhaset = 20456
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF3nlo.LHgrid'
     elseif((lhainput .ge. 20458) .and. (lhainput .le. 20459)) then
        lhaset = 20458
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF4nlo.LHgrid'
     elseif((lhainput .ge. 20460) .and. (lhainput .le. 20462)) then
        lhaset = 20460
        lhaname=lhapath(1:lhapathlen)//'/MRST2004qed.LHgrid'
     elseif((lhainput .ge. 20470) .and. (lhainput .le. 20471)) then
        lhaset = 20470
        lhaname=lhapath(1:lhapathlen)//'/MRST2004nnlo.LHgrid'
     elseif((lhainput .ge. 20550) .and. (lhainput .le. 20580)) then
        lhaset = 20550
        lhaname=lhapath(1:lhapathlen)//'/MRST2006nnlo.LHgrid'
        q2min = 1.0d0
        q2max = 1.0d09
        xmin = 1.0d-6
     elseif((lhainput .ge. 20650) .and. (lhainput .le. 20650)) then
        lhaset = 20650
        lhaname=lhapath(1:lhapathlen)//'/MRST2007lomod.LHgrid'
     elseif((lhainput .ge. 20651) .and. (lhainput .le. 20651)) then
        lhaset = 20651
        lhaname=lhapath(1:lhapathlen)//'/MRSTMCal.LHgrid'
     elseif((lhainput .ge. 29000) .and. (lhainput .le. 29003)) then
        lhaset = 29000
        lhaname=lhapath(1:lhapathlen)//'/MRST98.LHpdf'
     elseif((lhainput .ge. 29040) .and. (lhainput .le. 29045)) then
        lhaset = 29040
        lhaname=lhapath(1:lhapathlen)//'/MRST98lo.LHgrid'
     elseif((lhainput .ge. 29050) .and. (lhainput .le. 29055)) then
        lhaset = 29050
        lhaname=lhapath(1:lhapathlen)//'/MRST98nlo.LHgrid'
     elseif((lhainput .ge. 29060) .and. (lhainput .le. 29065)) then
        lhaset = 29060
        lhaname=lhapath(1:lhapathlen)//'/MRST98dis.LHgrid'
     elseif((lhainput .ge. 29070) .and. (lhainput .le. 29071)) then
        lhaset = 29070
        lhaname=lhapath(1:lhapathlen)//'/MRST98ht.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! Fermi Family
  elseif((lhainput .ge. 30000) .and. (lhainput .le. 39999)) then
     if((lhainput .ge. 30100) .and. (lhainput .le. 30200)) then
        lhaset = 30100
        lhaname=lhapath(1:lhapathlen)//'/Fermi2002_100.LHpdf'
     elseif((lhainput .ge. 31000) .and. (lhainput .le. 32000)) then
        lhaset = 31000
        lhaname=lhapath(1:lhapathlen)//'/Fermi2002_1000.LHpdf'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! Alekhin Family
  elseif((lhainput .ge. 40000) .and. (lhainput .le. 49999)) then
     if((lhainput .ge. 40100) .and. (lhainput .le. 40200)) then
        lhaset = 40100
        lhaname=lhapath(1:lhapathlen)//'/Alekhin_100.LHpdf'
     elseif((lhainput .ge. 41000) .and. (lhainput .le. 41999)) then
        lhaset = 41000
        lhaname=lhapath(1:lhapathlen)//'/Alekhin_1000.LHpdf'
     elseif((lhainput .ge. 40350) .and. (lhainput .le. 40367)) then
        lhaset = 40350
        lhaname=lhapath(1:lhapathlen)//'/a02m_lo.LHgrid'
        xmin = 1.0d-7
        q2min = 0.8d0
        q2max = 2.0d08
     elseif((lhainput .ge. 40450) .and. (lhainput .le. 40467)) then
        lhaset = 40450
        lhaname=lhapath(1:lhapathlen)//'/a02m_nlo.LHgrid'
        xmin = 1.0d-7
        q2min = 0.8d0
        q2max = 2.0d08
     elseif((lhainput .ge. 40550) .and. (lhainput .le. 40567)) then
        lhaset = 40550 
        lhaname=lhapath(1:lhapathlen)//'/a02m_nnlo.LHgrid'
        xmin = 1.0d-7
        q2min = 0.8d0
        q2max = 2.0d08
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! Botje Family
  elseif((lhainput .ge. 50000) .and. (lhainput .le. 59999)) then
     if((lhainput .ge. 50100) .and. (lhainput .le. 50200)) then
        lhaset = 50100
        lhaname=lhapath(1:lhapathlen)//'/Botje_100.LHpdf'
     elseif((lhainput .ge. 51000) .and. (lhainput .le. 51999)) then
        lhaset = 51000
        lhaname=lhapath(1:lhapathlen)//'/Botje_1000.LHpdf'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! ZEUS Family
  elseif((lhainput .ge. 60000) .and. (lhainput .le. 69999)) then
     q2min = 0.3d0
     q2max = 2.0d05
     if((lhainput .ge. 60000) .and. (lhainput .le. 60022)) then
        lhaset = 60000
        lhaname=lhapath(1:lhapathlen)//'/ZEUS2002_TR.LHpdf'
     elseif((lhainput .ge. 60100) .and. (lhainput .le. 60122)) then
        lhaset = 60100
        lhaname=lhapath(1:lhapathlen)//'/ZEUS2002_ZM.LHpdf'
     elseif((lhainput .ge. 60200) .and. (lhainput .le. 60222)) then
        lhaset = 60200
        lhaname=lhapath(1:lhapathlen)//'/ZEUS2002_FF.LHpdf'
     elseif((lhainput .ge. 60300) .and. (lhainput .le. 60322)) then
        lhaset = 60300
        lhaname=lhapath(1:lhapathlen)//'/ZEUS2005_ZJ.LHpdf'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! H1 Family
  elseif((lhainput .ge. 70000) .and. (lhainput .le. 79999)) then
     q2min = 1.5d0
     q2max = 3.5d04
     xmin = 5.7d-5
     if((lhainput .ge. 70050) .and. (lhainput .le. 70050)) then
        lhaset = 70050
        lhaname=lhapath(1:lhapathlen)//'/H12000ms.LHgrid'
     elseif((lhainput .ge. 70051) .and. (lhainput .le. 70070)) then
        lhaset = 70050
        lhaname=lhapath(1:lhapathlen)//'/H12000msE.LHgrid'
     elseif((lhainput .ge. 70150) .and. (lhainput .le. 70150)) then
        lhaset = 70150
        lhaname=lhapath(1:lhapathlen)//'/H12000dis.LHgrid'
     elseif((lhainput .ge. 70151) .and. (lhainput .le. 70170)) then
        lhaset = 70150
        lhaname=lhapath(1:lhapathlen)//'/H12000disE.LHgrid'
     elseif((lhainput .ge. 70250) .and. (lhainput .le. 70250)) then
        lhaset = 70250
        lhaname=lhapath(1:lhapathlen)//'/H12000lo.LHgrid'
     elseif((lhainput .ge. 70251) .and. (lhainput .le. 70270)) then
        lhaset = 70250
        lhaname=lhapath(1:lhapathlen)//'/H12000loE.LHgrid'
        ! Temporarily removed on returning to original H12000 files
        ! elseif((lhainput .ge. 70350) .and. (lhainput .le. 70350)) then
        ! lhaset = 70350
        ! lhaname=lhapath(1:lhapathlen)//'/H12000lo2.LHgrid'
        ! elseif((lhainput .ge. 70351) .and. (lhainput .le. 70370)) then
        ! lhaset = 70350
        ! lhaname=lhapath(1:lhapathlen)//'/H12000lo2E.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! GRV Family
  elseif((lhainput .ge. 80000) .and. (lhainput .le. 89999)) then
     q2min = 0.8d0
     q2max = 2.0d06
     xmin = 1.0d-9
     if((lhainput .ge. 80050) .and. (lhainput .le. 80051)) then
        lhaset = 80050
        lhaname=lhapath(1:lhapathlen)//'/GRV98nlo.LHgrid'
     elseif((lhainput .ge. 80060) .and. (lhainput .le. 80060)) then
        lhaset = 80060
        lhaname=lhapath(1:lhapathlen)//'/GRV98lo.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! 
     ! Pions
     ! 
     ! OW-PI Family
  elseif((lhainput .ge. 210) .and. (lhainput .le. 212)) then
     q2min = 4.0d0
     q2max = 2.0d03
     xmin = 5.0d-03
     xmax = 0.9998d0
     if((lhainput .ge. 210) .and. (lhainput .le. 212)) then
        lhaset = 210
        lhaname=lhapath(1:lhapathlen)//'/OWPI.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! SMRS-PI Family
  elseif((lhainput .ge. 230) .and. (lhainput .le. 233)) then
     q2min = 5.0d0
     q2max = 1.31d06
     xmin = 1.0d-05
     xmax = 0.9998d0
     if((lhainput .ge. 230) .and. (lhainput .le. 233)) then
        lhaset = 230
        lhaname=lhapath(1:lhapathlen)//'/SMRSPI.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! GRV-PI Family
  elseif((lhainput .ge. 250) .and. (lhainput .le. 252)) then
     q2max = 1.00d06
     xmin = 1.0d-05
     xmax = 0.9998d0
     if((lhainput .ge. 250) .and. (lhainput .le. 251)) then
        q2min = 3.0d-1
        lhaset = 250
        lhaname=lhapath(1:lhapathlen)//'/GRVPI1.LHgrid'
     elseif((lhainput .ge. 252) .and. (lhainput .le. 252)) then
        q2min = 2.5d-1
        lhaset = 252
        lhaname=lhapath(1:lhapathlen)//'/GRVPI0.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! ABFKW-PI Family
  elseif((lhainput .ge. 260) .and. (lhainput .le. 263)) then
     q2min = 2.0d0
     q2max = 1.00d08
     xmin = 1.0d-03
     xmax = 0.9998d0
     if((lhainput .ge. 260) .and. (lhainput .le. 263)) then
        lhaset = 260
        lhaname=lhapath(1:lhapathlen)//'/ABFKWPI.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! 
     ! photons
     ! 
     ! DO-G Family
  elseif((lhainput .ge. 310) .and. (lhainput .le. 312)) then
     q2min = 1.0d01
     q2max = 1.00d04
     xmin = 1.0d-05
     xmax = 0.9d0
     if((lhainput .ge. 310) .and. (lhainput .le. 311)) then
        lhaset = 310
        lhaname=lhapath(1:lhapathlen)//'/DOG0.LHgrid'
     elseif((lhainput .ge. 312) .and. (lhainput .le. 312)) then
        lhaset = 312
        lhaname=lhapath(1:lhapathlen)//'/DOG1.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! DG-G Family
  elseif((lhainput .ge. 320) .and. (lhainput .le. 324)) then
     xmin = 1.0d-05
     xmax = 0.9998d0
     lhaset = 320
     if((lhainput .ge. 320) .and. (lhainput .le. 321)) then
        q2min = 1.0d0
        q2max = 1.0d04
        ! lhaset = 320
        lhaname=lhapath(1:lhapathlen)//'/DGG.LHgrid'
     elseif((lhainput .ge. 322) .and. (lhainput .le. 322)) then
        q2min = 1.0d0
        q2max = 5.0d01
        ! lhaset = 322
        lhaname=lhapath(1:lhapathlen)//'/DGG.LHgrid'
     elseif((lhainput .ge. 323) .and. (lhainput .le. 323)) then
        q2min = 2.0d1
        q2max = 5.0d02
        ! lhaset = 323
        lhaname=lhapath(1:lhapathlen)//'/DGG.LHgrid'
     elseif((lhainput .ge. 324) .and. (lhainput .le. 324)) then
        q2min = 2.0d2
        q2max = 1.0d04
        ! lhaset = 324
        lhaname=lhapath(1:lhapathlen)//'/DGG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! LAC/GAL-G Family
  elseif((lhainput .ge. 330) .and. (lhainput .le. 334)) then
     q2min = 4.0d00
     q2max = 1.0d05
     xmin = 1.0d-04
     xmax = 0.9998d0
     lhaset = 330
     if((lhainput .ge. 330) .and. (lhainput .le. 332)) then
        ! lhaset = 330
        lhaname=lhapath(1:lhapathlen)//'/LACG.LHgrid'
     elseif((lhainput .ge. 333) .and. (lhainput .le. 333)) then
        q2min = 1.0d00
        ! lhaset = 333
        lhaname=lhapath(1:lhapathlen)//'/LACG.LHgrid'
     elseif((lhainput .ge. 334) .and. (lhainput .le. 334)) then
        q2min = 4.0d00
        ! lhaset = 334
        lhaname=lhapath(1:lhapathlen)//'/LACG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! GSG/GSG96-G Family
  elseif((lhainput .ge. 340) .and. (lhainput .le. 345)) then
     q2min = 5.3d00
     q2max = 1.0d08
     xmin = 5.0d-04
     xmax = 0.9998d0
     if((lhainput .ge. 340) .and. (lhainput .le. 341)) then
        lhaset = 340
        lhaname=lhapath(1:lhapathlen)//'/GSG1.LHgrid'
     elseif((lhainput .ge. 342) .and. (lhainput .le. 343)) then
        lhaset = 342
        lhaname=lhapath(1:lhapathlen)//'/GSG0.LHgrid'
     elseif((lhainput .ge. 344) .and. (lhainput .le. 344)) then
        lhaset = 344
        lhaname=lhapath(1:lhapathlen)//'/GSG961.LHgrid'
     elseif((lhainput .ge. 345) .and. (lhainput .le. 345)) then
        lhaset = 345
        lhaname=lhapath(1:lhapathlen)//'/GSG960.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! GRV-G Family
  elseif((lhainput .ge. 350) .and. (lhainput .le. 354)) then
     q2min = 3.0d-1
     q2max = 1.0d06
     xmin = 1.0d-05
     xmax = 0.9998d0
     if((lhainput .ge. 350) .and. (lhainput .le. 352)) then
        lhaset = 350
        lhaname=lhapath(1:lhapathlen)//'/GRVG1.LHgrid'
     elseif((lhainput .ge. 353) .and. (lhainput .le. 353)) then
        q2min = 2.5d-1
        lhaset = 353
        lhaname=lhapath(1:lhapathlen)//'/GRVG0.LHgrid'
     elseif((lhainput .ge. 354) .and. (lhainput .le. 354)) then
        q2min = 6.0d-1
        q2max = 5.0d04
        lhaset = 354
        lhaname=lhapath(1:lhapathlen)//'/GRVG0.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! ACFGP-G Family
  elseif((lhainput .ge. 360) .and. (lhainput .le. 363)) then
     q2min = 2.0d00
     q2max = 5.5d05
     xmin = 1.37d-03
     xmax = 0.9998d0
     if((lhainput .ge. 360) .and. (lhainput .le. 363)) then
        lhaset = 360
        lhaname=lhapath(1:lhapathlen)//'/ACFGPG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! WHIT-G Family
  elseif((lhainput .ge. 380) .and. (lhainput .le. 386)) then
     q2min = 4.0d00
     q2max = 2.5d03
     xmin = 1.0d-03
     xmax = 0.9998d0
     if((lhainput .ge. 380) .and. (lhainput .le. 386)) then
        lhaset = 380
        lhaname=lhapath(1:lhapathlen)//'/WHITG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! SAS-G Family
  elseif ((lhainput .ge. 390) .and. (lhainput .le. 398)) then
     q2max = 5.0d04
     xmin = 1.0d-05
     xmax = 0.9998d0
     lhaset = 390
     if ((lhainput .ge. 390) .and. (lhainput .le. 392)) then
        q2min = 3.6d-1
        ! lhaset = 390
        lhaname=lhapath(1:lhapathlen)//'/SASG.LHgrid'
     elseif((lhainput .ge. 393) .and. (lhainput .le. 394)) then
        q2min = 4.0d00
        !           lhaset = 393
        lhaname=lhapath(1:lhapathlen)//'/SASG.LHgrid'
     elseif((lhainput .ge. 395) .and. (lhainput .le. 396)) then
        q2min = 3.6d-1
        !           lhaset = 395
        lhaname=lhapath(1:lhapathlen)//'/SASG.LHgrid'
     elseif((lhainput .ge. 397) .and. (lhainput .le. 398)) then
        q2min = 4.0d00
        !           lhaset = 397
        lhaname=lhapath(1:lhapathlen)//'/SASG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! Unknown Family ?! Giving up
  else
     write(lhaprint,5150)  lhaset
     stop
  endif

  lhamemb=lhainput-lhaset
  ! Now work out if we have already called this set/member
  iset = 0
  do j=1,nsets
     if (lhaname.eq.lhanames(j).and. &
          lhamemb.eq.lhamembers(j)) then
        iset = j
     endif
  enddo
  if (iset.eq.0) then
     nsets=nsets+1
     if (nsets.gt.nmxset) then
        if (LHASILENT.ne.1) then
           print *, "WARNING: too many sets initialised"
           print *,"overwriting from set 1 again"
        endif
        nsets = 1
        ! stop
     endif
     iset=nsets
     lhanames(iset)=lhaname
     lhanumbers(iset)=lhainput
     lhamembers(iset)=lhamemb
     xxmin(iset)=xmin
     xxmax(iset)=xmax
     qq2min(iset)=q2min
     qq2max(iset)=q2max
     call initpdfsetm(iset,lhaname)
     call numberpdfm(iset,lhaallmem)
     if(lhasilent .ne. 1) then
        write(lhaprint,5151)
        write(lhaprint,5152) lhaname
        write(lhaprint,5153) lhaallmem
        write(lhaprint,5154)
     endif
     if ((lhamemb.lt.0) .or. (lhamemb.gt.lhaallmem)) then
        write(lhaprint,5155)  lhamemb
        write(lhaprint,5156)  lhaallmem
        stop
     endif

     ! print *,'calling initpdf',lhamemb 
     ! print *,'calling initpdfm ',iset,lhaname,lhamemb
     ! print *,'LHAGLUE .... initializing set,member ',iset,lhamemb
     call initpdfm(iset,lhamemb)
  endif
  !  the rest is done every time pdfset is called
  !print *,'setting nset to:',iset
  call setnset(iset)
  call setnmem(iset,lhamemb)
  xmin = xxmin(iset)
  xmax = xxmax(iset)
  q2min=qq2min(iset)
  q2max=qq2max(iset)
  call GetLam4M(iset,LHAMEMB,qcdlha4)
  call GetLam5M(iset,LHAMEMB,qcdlha5)

  QMZ = 91.1876D0
  alphasLHA = alphasPDFM(iset,QMZ)
  if(lhasilent .ne. 1) write(lhaprint,5158) alphasLHA

  if(lhaparm(17).EQ.'LHAPDF') then
     nptypepdfl = 1      ! Proton PDFs
     nflpdfl = 4
     qcdl4 = qcdlha4
     qcdl5 = qcdlha5
     if (LHASILENT .NE. 1) write(lhaprint,5159) qcdl4, qcdl5
  else
     nptypepdfl = 1      ! Proton PDFs
     nflpdfl = 4
     alambda = 0.192d0
     qcdlha4 = alambda
     qcdlha5 = alambda
     if (parm(1).EQ.'NPTYPE') then        !  PYTHIA
        qcdl4 = alambda
        qcdl5 = alambda
     endif
  endif

  if(parm(1).EQ.'NPTYPE') then  !  herwig
     axmin  = xmin
     axmax  = xmax
     aq2min = q2min
     aq2max = q2max
  endif
  ! Formats for initialization information.
5150 format(1X,'WRONG LHAPDF set number =',I12,' given! STOP EXE!')
5151 format(1X,'==============================================')
5152 format(1X,'PDFset name ',A80)
5153 format(1X,'with ',I10,' members')
5154 format(1X,'====  initialized. ===========================')
5155 format(1X,'LHAPDF problem => YOU asked for member = ',I10)
5156 format(1X,'Valid range is: 0 - ',I10,' Execution stopped.')
  !5157 format(1X,'Number of flavors for PDF is:',I4)
5158 format(1X,'Strong coupling at Mz for PDF is:',F9.5)
5159 format(1X,'Will use for PYTHIA QCDL4, QCDL5:',2F9.5)

  return
end subroutine pdfset


!********************************************************************
! -- STRUCTA
! -- copy of PDFLIB to use the eks98 nuclear correction factors

subroutine structa(x,q,a,upv,dnv,usea,dsea,str,chm,bot,top,glu)
  implicit double precision (a-h,o-z)
  character*20 lparm
  call getlhaparm(15,lparm)
  if(lparm.eq.'EPS08') then
     call eps08(x,q,a,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
  else
     call eks98(x,q,a,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
  endif
  call structm(x,q,upv,dnv,usea,dsea,str,chm,bot,top,glu)
  upv = ruv*upv
  dnv = rdv*dnv
  usea = ru*usea
  dsea = rd*dsea
  str = rs*str
  chm = rc*chm
  bot = rb*bot
  top = rt*top
  glu = rg*glu
  return
end subroutine structa


!*********************************************************************
! STRUCTM
! Gives parton distributions according to the LHAPDF interface.
! Two evolution codes used:
!   EVLCTEQ for CTEQ PDF sets
!   QCDNUM  for Other PDF sets
! 
! Author: Dimitri Bourilkov  bourilkov@mailaps.org
! 
! v4.0  21-Mar-2005  Photon/pion/new p PDFs, updated for LHAPDF v4
! v3.0  23-Jan-2004
! 
! interface to LHAPDF library
subroutine structm(dx,dq,upv,dnv,usea,dsea,str,chm,bot,top,glu)
  ! double precision and integer declarations.
  implicit double precision(a-h, o-z)
  implicit integer(i-n)
  ! Common blocks
  include 'commonlhapdf.inc'
  include 'commonlhasets.inc'
  include 'commonlhacontrol.inc'
  include 'commonlhaglsta.inc'
  ! interface to lhapdflib.
  double precision qcdlha4, qcdlha5
  integer nfllha
  common/lhapdfr/qcdlha4, qcdlha5, nfllha
  save /lhapdfr/
  integer lhaextrp
  common/lhapdfe/lhaextrp
  save /lhapdfe/
  ! interface to pdflib.
  common/w50513/xmin,xmax,q2min,q2max
  save /w50513/
  double precision xmin,xmax,q2min,q2max
  ! local variables
  double precision upv,dnv,usea,dsea,str,chm,bot,top,glu
  double precision dx,dq,x,q,f(-6:6),photon,gluino

  x = dx
  q = dq
  q2 = q**2
  ! statistics
  if(lhaparm(16).ne.'NOSTAT') then
     totnum = totnum+1.d0
     if(x .lt. xmin) xminnum = xminnum+1.d0
     if(x .gt. xmax) xmaxnum = xmaxnum+1.d0
     if(q2 .lt. q2min) q2minnum = q2minnum+1.d0
     if(q2 .gt. q2max) q2maxnum = q2maxnum+1.d0
  endif

  ! range of validity e.g. 10^-6 < x < 1, q2min < q^2 extended by
  ! freezing x*f(x,q2) at borders.
  if(lhaextrp .ne. 1) then    ! safe mode == "freeze"
     xin=max(xmin,min(xmax,x))
     q=sqrt(max(0d0,q2min,min(q2max,q2)))
  else                        ! adventurous mode == own risk !
     xin=x
  endif

  call getnset(iset)
  !print *,'calling evolvepdfm:',iset

  ! fix to allow STRUCTM to work for photon PDFs (Herwig does this)
  ! set P2 = 0.0d0 and IP2 = 0
  if(lhanumbers(iset).ge.300.and.lhanumbers(iset).le.399) then  
     p2 = 0.0d0
     ip2 = 0
     call evolvepdfpm(iset,xin,q,p2,ip2,f)
  else if (lhanumbers(iset).ge.20460.and.lhanumbers(iset).le.20462) then
     call evolvepdfphotonm(iset,xin,q,f,photon)
  else if (lhanumbers(iset).ge.20670.and.lhanumbers(iset).le.20677) then
     call evolvepdfgluinom(iset,xin,q,f,gluino)
  else
     call evolvepdfm(iset,xin,q,f)
  endif
  glu = f(0)
  dsea = f(-1)
  dnv = f(1) - dsea
  usea = f(-2)
  upv = f(2) - usea
  str = f(3)
  chm = f(4)
  bot = f(5)
  top = f(6)

  return
end subroutine structm


!*********************************************************************
! STRUCTP
! Gives parton distributions according to the LHAPDF interface.
! Used for photons.
! 
! v4.0  21-Mar-2005  Photon/pion/new p PDFs, updated for LHAPDF v4
! 
! Interface to LHAPDF library
subroutine structp(dx,dq2,p2,ip2,upv,dnv,usea,dsea,str,chm,bot,top,glu)
  ! Double precision and integer declarations.
  implicit double precision(a-h, o-z)
  implicit integer(i-n)
  ! Common blocks
  include 'parmsetup.inc'
  include 'commonlhapdf.inc'
  include 'commonlhacontrol.inc'
  include 'commonlhaglsta.inc'
  ! Interface to LHAPDFLIB.
  double precision qcdlha4, qcdlha5
  integer nfllha
  common/lhapdfr/qcdlha4, qcdlha5, nfllha
  save /lhapdfr/
  integer lhaextrp
  common/lhapdfe/lhaextrp
  save /lhapdfe/
  ! Interface to PDFLIB.
  common/w50513/xmin,xmax,q2min,q2max
  save /w50513/
  double precision xmin,xmax,q2min,q2max
  ! Local variables
  double precision upv,dnv,usea,dsea,str,chm,bot,top,glu
  double precision dx,dq2,q2,x,q,f(-6:6)

  x = dx
  q2 = dq2
  ! Statistics
  if(lhaparm(16).ne.'NOSTAT') then
     totnup = totnup+1.d0
     if(x .lt. xmin) xminnup = xminnup+1.d0
     if(x .gt. xmax) xmaxnup = xmaxnup+1.d0
     if(q2 .lt. q2min) q2minnup = q2minnup+1.d0
     if(q2 .gt. q2max) q2maxnup = q2maxnup+1.d0
  endif

  ! Range of validity e.g. 10^-6 < x < 1, Q2MIN < Q^2 extended by
  ! freezing x*f(x,Q2) at borders.
  q = dsqrt(q2)
  if(lhaextrp .ne. 1) then    ! safe mode == "freeze"
     xin=max(xmin,min(xmax,x))
     q=sqrt(max(0d0,q2min,min(q2max,q2)))
  else                        ! adventurous mode == OWN RISK !
     xin=x
  endif
  call getnset(iset)
  call evolvepdfpm(iset,xin,q,p2,ip2,f)
  glu = f(0)
  dsea = f(-1)
  dnv = f(1) - dsea
  usea = f(-2)
  upv = f(2) - usea
  str = f(3)
  chm = f(4)
  bot = f(5)
  top = f(6)
  return
end subroutine structp


!*********************************************************************
! PDFSTA
! For statistics ON structure functions (under/over-flows)
! 
! Author: Dimitri Bourilkov  bourilkov@mailaps.org
! 
! 
! first introduced in v4.0  28-Apr-2005 
! 
subroutine pdfsta
  ! Double precision and integer declarations.
  implicit double precision(a-h, o-z)
  implicit integer(i-n)
  ! Common blocks
  include 'commonlhaglsta.inc'
  ! Interface to LHAPDFLIB.

  print *
  print *,'===== PDFSTA statistics for PDF under/over-flows ===='
  print *
  print *,'====== STRUCTM statistics for nucleon/pion PDFs ====='
  print *
  print *,'  total # of calls ',TOTNUM
  if(totnum .gt. 0.d0) then
     percbelow = 100.d0*xminnum/totnum
     percabove = 100.d0*xmaxnum/totnum
     print *,'  X  below PDF min ',xminnum,' or ',percbelow, ' %'
     print *,'  X  above PDF max ',xmaxnum,' or ',percabove, ' %'
     percbelow = 100.d0*q2minnum/totnum
     percabove = 100.d0*q2maxnum/totnum
     print *,'  Q2 below PDF min ',q2minnum,' or ',percbelow, ' %'
     print *,'  Q2 above PDF max ',q2maxnum,' or ',percabove, ' %'
  endif
  print *
  print *,'========= STRUCTP statistics for photon PDFs ========'
  print *
  print *,'  total # of calls ',totnup
  if(totnup .gt. 0.d0) then
     percbelow = 100.d0*xminnup/totnup
     percabove = 100.d0*xmaxnup/totnup
     print *,'  X  below PDF min ',xminnup,' or ',percbelow, ' %'
     print *,'  X  above PDF max ',xmaxnup,' or ',percabove, ' %'
     percbelow = 100.d0*q2minnup/totnup
     percabove = 100.d0*q2maxnup/totnup
     print *,'  Q2 below PDF min ',q2minnup,' or ',percbelow, ' %'
     print *,'  Q2 above PDF max ',q2maxnup,' or ',percabove, ' %'
  endif
  print *
  return
end subroutine pdfsta


subroutine pftopdg(dx,dscale,dxpdf)
  !include "pdf/expdp.inc"
  double precision dx,dscale,dupv,ddnv,dusea,ddsea,dstr,dchm,dbot,dtop,dgl,dxpdf(-6:6)
  ! Call STRUCTM in PDFLIB to get flavour content
  call structm(dx,dscale,dupv,ddnv,dusea,ddsea,dstr,dchm,dbot,dtop,dgl)
  ! Convert flavour convention of PDFLIB to PDG convention
  dxpdf(0) = dgl
  dxpdf(1) = ddnv + ddsea
  dxpdf(2) = dupv + dusea
  dxpdf(3) = dstr
  dxpdf(4) = dchm
  dxpdf(5) = dbot
  dxpdf(6) = dtop
  dxpdf(-1) = ddsea
  dxpdf(-2) = dusea
  dxpdf(-3) = dstr
  dxpdf(-4) = dchm
  dxpdf(-5) = dbot
  dxpdf(-6) = dtop
  return
end subroutine pftopdg
