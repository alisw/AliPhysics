*
* $Id$
*
* $Log$
* Revision 1.7  1998/09/25 09:33:41  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.6  1997/09/02 15:16:13  mclareni
* WINNT corrections
*
* Revision 1.5  1997/03/14 17:21:20  mclareni
* WNT mods
*
* Revision 1.4  1997/01/15 17:41:16  cernlib
* disable section calling rzstrip
*
* Revision 1.3  1996/10/17 09:35:04  cernlib
* make cfstat an integer function
*
* Revision 1.2  1996/10/16 13:02:53  cernlib
* Use CFSTAT with CFIO instead of STATF (for RFIO)
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZOPEN(LUNIN,CHDIR,CFNAME,CHOPTT,LRECL,ISTAT)
*
************************************************************************
*
*      Open a ZEBRA/RZ file.
*
* Input:
*      LUNIT    Logical unit number
*      CFNAME   File name
*      CHOPT    Character variable specifying the option
*               ' '  default, open file in readonly mode
*               'L'  create file with relative organization (VAX only)
*               'N'  open a new file
*               'S'  open file in shared readonly mode
*               'U'  open file in update mode
*               'SU' open file in shared update mode
*               '1'  open file read/write assume single user
*               'V'  open new RZ file on VSAM file
*               'W'  return in CHDIR directory name include
*                    logical unit number
*               'Y'  suppress LRECL consistency check
*               'P'  Preserve case of file name (Unix)
*               'C'  Use C I/O (Unix, VMS)
*               'X'  Exchange mode file
*
*      *LRECL*  Record length, if zero determine LRECL from input file
*
* Output:
*      CHDIR      Character string containing decoded logical unit number
*      *LRECL*    Determined record length
*      ISTAT      Status return code
*      IQUEST(10) LRECL
*      IQUEST(11) LUNPTR - C I/O pointer
*      IQUEST(12) Exchange mode flag - set if IOPTX.ne.0
*                 or if exchange mode bit is set in file
*
* 
* Called by <USER>
*
*  Author  : R.Brun,J.Shiers
*  Written : 03.05.86
*  Last mod: See below
*
* Changes        Date    Comments
* G.Folger     96/10/16  Use CFSTAT with CFIO instead of STATF (for RFIO)
* V.Fine       96/07/30  Disable using C I/O to check Fortran I/O under Windows NT
* J.Shiers     95/06/20  Warning for auto-recl determination only for >8192
* J.Shiers     95/05/24  Use STATF instead of INQUIRE in case of IOPTC
* J.Shiers     94/09/21  Increase chopt, correct LUN string
* F.Rademakers 94/08/29  Added protection in case CFOPEN fails (no read perm).
*                        Return correct LUN string in case of option W and
*                        C I/O. Print correct error message in case user
*                        opens file with wrong LRECL.
* J.Shiers     94/08/18  Increase buffer from 8192 to 8192+512
*                        to permit record length determination upto 8192
* R.Brun       94/06/10  Introduce file striping for PIAF
*                        New routine RZSTRIP called
* J.Shiers     94/05/18  Add QMDOS flag for parameter NWORD
* V.Fine       94/02/07  DEC flag to use DEC Fortran and Windows/NT
* A.Lomov      93/09/14  No longer force exchange mode for LINUX
* J.Shiers     93/06/30  Set IQUEST(12) to IMODEX
* J.Shiers     93/02/15  Set IMODEH to 0
* J.Shiers     92/10/19  Alpha compatible record length determination
* J.Shiers     92/07/13  QMDOS, QCFIO flags
* J.Shiers     92/07/07  Translate filename to uppercase before VMQFIL
* J.Shiers     92/03/03  Incorporate MVS mods as follows:
* J.Shiers     92/02/26  Use VMQFIL on VM systems to determine record length
* J.Shiers     92/02/21  Mods for DESY (length of prefix)
* J.Shiers     91/11/27  Add ACTION=READ/READWRITE in VM open statements
* J.Shiers     91/11/11  Add C I/O support
* M.Marquina   91/06/26  Integrate MVS version
#if defined(CERNLIB_IBMMVS)
* IBMMVS changes:
*    -   IUNIT auf 33XX fuer GSI
*    -   VSAM Option fuer MVS   (OPT V) (VSAM Relative record dataset)
*        To use this option, a dataclass RZFILE must be defined via
*        SMS with appropriate parameters for default size and extends.
*        at GSI we use also a dataclass RZDATA for larger files. These
*        files contain a .RZDATA in their name.
*        At GSI we use the following definitions for SMS :
*
* DATACLAS    STORCLAS    MGMTCLAS   STORGRP  RECORG SPACE
*
* RZFILE(def) STANDARD    USERDATA   NORMALx  VSAM RR 128,128
* RZDATA      RAWDATA     RAWDATA    RAWDATA  VSAM RR 512,512
*
*
#endif
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzclun.inc"
#include "zebra/quest.inc"
#include "zebra/rzbuff.inc"
#if defined(CERNLIB_QMUIX)
#include "zebra/rzcstr.inc"
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_QMAPO)||defined(CERNLIB_QMAPO9)||defined(CERNLIB_QMUIX)||defined(CERNLIB_WINNT)
      integer cfstat,statf,info(12)
#endif
      CHARACTER*(*) CFNAME,CHDIR,CHOPTT
      CHARACTER*9   SPACES
      CHARACTER*8   STAT
      CHARACTER*36  CHOPT
      CHARACTER*255 CHFILE
      LOGICAL       IEXIST
      CHARACTER*4   CHOPE
#if defined(CERNLIB_IBMVM)
      CHARACTER*13  CHTIME
      CHARACTER*1   RECFM
#endif
#if defined(CERNLIB_QMVAX)
      CHARACTER*10  CHORG
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCVX)||defined(CERNLIB_QMIRT)||defined(CERNLIB_QMIBM)||defined(CERNLIB_QMSGI)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMSUN)||defined(CERNLIB_QMALT)||defined(CERNLIB_QMHPX)||defined(CERNLIB_QMIBX)
      PARAMETER (NWORD = 8704)
#endif
#if (defined(CERNLIB_QMAPO)||defined(CERNLIB_QMVAX)||defined(CERNLIB_QMNXT)||defined(CERNLIB_QMLNX))&&(!defined(CERNLIB_QMAPO9))
      PARAMETER (NWORD = 8704)
#endif
#if defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
      PARAMETER (NWORD = 8704)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
      PARAMETER (IBYTES=8)
#endif
#if defined(CERNLIB_QMCRY)
      DIMENSION IBUFF(8704)
#endif
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMSGI)||defined(CERNLIB_QMVMI)||(defined(CERNLIB_QFDEC))
      PARAMETER (IBYTES=1)
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QMVAX))&&(!defined(CERNLIB_QMSGI))&&(!defined(CERNLIB_QMVMI))&&(!defined(CERNLIB_QFDEC))
      PARAMETER (IBYTES=4)
#endif
#if defined(CERNLIB_QMIBM)
      CHARACTER*9 CHACT
#endif
#if defined(CERNLIB_QMIBMFVS)
      CHARACTER*5 CHREC
      CHARACTER*3 CHLUN
#endif
#if defined(CERNLIB_IBMMVS)
      CHARACTER*20   PREFIX
      CHARACTER*4    IUNIT, ITRK, IRECFM
      CHARACTER*8 cRecl
      CHARACTER*256  cTSO
      CHARACTER*20   cDataClas
* dummy record for VSAM
      INTEGER RECORD(2048) /2048*0/
      INTEGER        mvsams, kdffil
      LOGICAL*4      EXS, OPN
      INTEGER*4      ISPACE(3), IDCB(2)
*
      DATA           ISPACE / 30, 0, 0 /
      DATA           IDCB   / 4096, 0 /
      DATA           ITRK   / 'TRK' /, IRECFM / 'F' /
#endif
#if (defined(CERNLIB_IBMMVS))&&(!defined(CERNLIB_NEWLIB))&&(!defined(CERNLIB_GSI))
      DATA           IUNIT  / 'HSM' /
#endif
#if (defined(CERNLIB_IBMMVS))&&(!defined(CERNLIB_NEWLIB))&&(defined(CERNLIB_GSI))
      DATA           IUNIT  / '33XX' /
#endif
#if (defined(CERNLIB_IBMMVS))&&(defined(CERNLIB_NEWLIB))
      DATA           IUNIT  / 'FAST' /
#endif
*
*-----------------------------------------------------------------------
*
      CHOPT=CHOPTT
      CALL CLTOU(CHOPT)
 
      IOPT1=INDEX(CHOPT,'1')
      IOPTC=INDEX(CHOPT,'C')
      IOPTL=INDEX(CHOPT,'L')
      IOPTN=INDEX(CHOPT,'N')
      IOPTS=INDEX(CHOPT,'S')
      IOPTP=INDEX(CHOPT,'P')
      IOPTU=INDEX(CHOPT,'U')
      IOPTV=INDEX(CHOPT,'V')
      IOPTW=INDEX(CHOPT,'W')
      IOPTX=INDEX(CHOPT,'X')
      IOPTY=INDEX(CHOPT,'Y')
 
      LUNIT=LUNIN
      IQUEST(10) = 0
      IQUEST(11) = 0
      IQUEST(12) = 0
      IMODEC     = IOPTC
      IMODEX     = IOPTX
 
      LRECI      = LRECL
      LRECL2     = 0
#if (defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))
C     RZfile with Exchange mode for NonPPC-LINUX
*     IMODEX     = 1
#endif
      IMODEH     = 0
#if !defined(CERNLIB_QCFIO)
*SELF,IF=-QMCRY,IF=-QMVAX,IF=-QMCV64,IF=-QMAPO,IF=-QMAPO9,IF=-QMUIX,IF=-QMDOS.
      IF(IOPTC.NE.0) THEN
         WRITE(IQPRNT,*) 'RZOPEN. option C ignored - valid only ',
     +      'for MSDOS, Unix and VMS systems'
         IOPTC = 0
      ENDIF
#endif
#if !defined(CERNLIB_QMIBMFVS)
      CHFILE=CFNAME
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(!defined(CERNLIB_IBMMVS))
      CHFILE='/'//CFNAME
*
*     Search all disks if filemode not specified
*
      LCHF = LENOCC(CHFILE)
      CALL CTRANS('.',' ',CHFILE,1,LCHF)
      IF(INDEX(CHFILE(1:LCHF),' ').EQ.INDEXB(CHFILE(1:LCHF),' '))
     +   CHFILE = CHFILE(1:LCHF) // ' *'
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))
*--   Construct MVS file name
*--   Don't add prefix if the first character of file name is a dot
      CALL KPREFI (PREFIX, NCHPRE)
      IF ( CFNAME(1:1) .EQ. '.' ) THEN
         NCH    = LEN(CFNAME)
         CHFILE = '/'//CFNAME(2:NCH)
      ELSE
         CHFILE = '/'//PREFIX(1:NCHPRE)//CFNAME
      ENDIF
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCVX)||defined(CERNLIB_QMAPO)||defined(CERNLIB_QMAPO9)||defined(CERNLIB_QMUIX)||defined(CERNLIB_QMDOS)||defined(CERNLIB_QMLNX)||defined(CERNLIB_WINNT)
      IF(IOPTP.EQ.0)CALL CUTOL(CHFILE)
#endif
#if defined(CERNLIB_QMVAX)
      IF(IOPTL.NE.0) THEN
         CHORG = 'RELATIVE'
      ELSE
         CHORG = 'SEQUENTIAL'
      ENDIF
#endif
      IPASS=0
      kstrip=0
   10 CONTINUE
      IF(IOPTN.NE.0)THEN
#if !defined(CERNLIB_QMVAX)
         STAT='UNKNOWN'
#endif
#if defined(CERNLIB_QMVAX)
         STAT='NEW'
#endif
         IF(LRECI.LE.0) THEN
            WRITE(IQPRNT,10000) LRECI
10000 FORMAT(' RZOPEN. - invalid record length: ',I6)
            ISTAT = 1
            GOTO 70
         ELSEIF(LRECI.GE.8191) THEN
            WRITE(IQPRNT,10100) LRECI
10100 FORMAT(' RZOPEN. record length:',I6,
     +       ' > maximum safe value (8191 words).')
            IF(LRECI.GT.8192) WRITE(IQPRNT,10200)
10200 FORMAT(' RZOPEN. Automatic record length determination will not',
     +       ' work with this file.')
            WRITE(IQPRNT,10300)
10300 FORMAT(' RZOPEN. You may have problems transferring your',
     +       ' file to other systems ',/,
     +       '         or writing it to tape.')
         ENDIF
      ELSE
#if defined(CERNLIB_QMAPO)
         IF(IOPT1.EQ.0.AND.IOPTU.EQ.0)THEN
            STAT='READONLY'
         ELSE
            STAT='OLD'
         ENDIF
#endif
#if !defined(CERNLIB_QMAPO)
         STAT='OLD'
#endif
*
*     File should already exist. Issue inquire even if IOPTC
*
         LENF = LENOCC(CHFILE)
 
#ifndef CERNLIB_QFMSOFT
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_QMAPO)||defined(CERNLIB_QMAPO9)||defined(CERNLIB_QMUIX)||defined(CERNLIB_WINNT)
         IF(IOPTC.EQ.0) THEN
#endif
#endif
            INQUIRE(FILE=CHFILE,EXIST=IEXIST)
            ISTATF = 0

#ifndef CERNLIB_QFMSOFT
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_QMAPO)||defined(CERNLIB_QMAPO9)||defined(CERNLIB_QMUIX)||defined(CERNLIB_WINNT)
         ELSE
            IEXIST = CFSTAT(CHFILE(1:LENF),INFO).EQ.0
            ISTATF = 1
         ENDIF
#endif
#endif
 
         IF(.NOT.IEXIST) THEN
            WRITE(IQPRNT,*) 'RZOPEN. Error - input file ',
#if defined(CERNLIB_IBM)
     +         CHFILE(2:LENF),' does not exist'
#endif
#if !defined(CERNLIB_IBM)
     +         CHFILE(1:LENF),' does not exist'
#endif
            ISTAT = 2
            GOTO 70
         ENDIF
#if defined(CERNLIB_QMVAX)
         IF(IOPTC.EQ.0) INQUIRE(FILE=CHFILE,ORGANIZATION=CHORG)
#endif
         IF(LRECL.EQ.0) THEN
*
*     LRECL=0 was specified, try to determine correct record length
*     from the file itself.
*
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
            ICOUNT = NWORD
            IF(IOPTC.EQ.0) THEN
#endif
#if defined(CERNLIB_QMAPO9)
               OPEN(UNIT=LUNIT,FILE=CHFILE,FORM='UNFORMATTED',
     +              STATUS='READONLY')
#endif
#if defined(CERNLIB_QMVAX)
               INQUIRE(FILE=CHFILE,RECL=LRECL)
               LRECL = LRECL / 4
               GOTO 40
#endif
#if defined(CERNLIB_QMAPO9)
               INQUIRE(UNIT=LUNIT,RECL=LRECL)
               CLOSE(UNIT=LUNIT)
               GOTO 40
#endif
#if (defined(CERNLIB_QMIBM))&&(defined(CERNLIB_IBMVM))
*
*     Use KERNLIB routine VMQFIL to obtain record length
*     Previous method (below) does not work for VMSTAGEd files
*
               CALL CLTOU(CHFILE(1:LCHF))
               CALL VMQFIL(CHFILE(2:),RECFM,LRECL,NRECS,NBLOCKS,CHTIME,
     +                     ISTAT,IRC)
               LRECL = LRECL/4
               GOTO 40
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(!defined(CERNLIB_IBMVM))
*
*      Open file for READ to determine record length. JDS
*
               OPEN(UNIT=LUNIT,FILE=CHFILE,FORM='UNFORMATTED',
     +              ACTION='READ',
     +              STATUS=STAT)
               READ(UNIT=LUNIT,NUM=LRECL) ITEST
               LRECL = LRECL/4
*      Close...
               CLOSE(LUNIT)
               GOTO 40
#endif
*
*     Record length is stored in file but in record # 2
*     (rec # 1 is used for locks)
*
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
               IF(IMODEX.NE.0) ICOUNT = ICOUNT / 2
#endif
#if defined(CERNLIB_QMCRY)
               OPEN(LUNIT,FILE=CHFILE,FORM='UNFORMATTED',STATUS='OLD',
     +              RECL=IBYTES*ICOUNT,ACCESS='DIRECT')
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
            ENDIF
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
   20       NREAD  = ICOUNT
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
            IF(IOPTC.EQ.0) THEN
#endif
#if (defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL))&&(!defined(CERNLIB_QMAPO))
               OPEN(LUNIT,FILE=CHFILE,FORM='UNFORMATTED',STATUS='OLD',
     +              RECL=IBYTES*NREAD,ACCESS='DIRECT',IOSTAT=ISTAT)
               IF(ISTAT.NE.0)GOTO 60
#endif
#if defined(CERNLIB_QMAPO)
               OPEN(LUNIT,FILE=CHFILE,FORM='UNFORMATTED',
     +              STATUS='READONLY',RECL=IBYTES*NREAD,ACCESS='DIRECT'
     +              ,IOSTAT=ISTAT)
               IF(ISTAT.NE.0)GOTO 60
#endif
#if defined(CERNLIB_QMCRY)
               IF(IOPTX.EQ.0) THEN
                  CALL READ(LUNIT,ITEST,NREAD,IOS,NUS)
               ELSE
                  CALL READ(LUNIT,IBUFF,NREAD,IOS,NUS)
               ENDIF
               IF(IOS.EQ.2) THEN
#endif
#if defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
               READ(LUNIT,REC=1,IOSTAT=IOS) (ITEST(JW),JW=1,NREAD)
               IF(IOS.NE.0) THEN
                  CLOSE(LUNIT)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
                  ICOUNT = ICOUNT * .75
                  GOTO 20
               ENDIF
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
            ELSE
               CALL CFOPEN(LUNPTR,0,NREAD,'r',0,CHFILE,IOS)
               IF (IOS .NE. 0) THEN
                  ISTAT = -1
                  GOTO 70
               ENDIF
               NWTAK = NREAD
#endif
#if defined(CERNLIB_QMCRY)
               IF(IOPTX.NE.0) THEN
                  CALL CFGET(LUNPTR,0,NREAD,NWTAK,IBUFF,IOS)
               ELSE
                  CALL CFGET(LUNPTR,0,NREAD,NWTAK,ITEST,IOS)
               ENDIF
#endif
#if defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
               CALL CFGET(LUNPTR,0,NREAD,NWTAK,ITEST,IOS)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
               IF(IOS.NE.0) THEN
                  CALL CFCLOS(LUNPTR,0)
                  ICOUNT = ICOUNT * .75
                  GOTO 20
               ENDIF
            ENDIF
#endif
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
            IF(IOPTX.NE.0) CALL VXINVB(ITEST(1),NREAD)
#endif
#if defined(CERNLIB_QMCRY)
            IF(IOPTX.NE.0) THEN
               CALL UNPACK(IBUFF(1),32,ITEST(1),NREAD)
            ENDIF
#endif
#if defined(CERNLIB_QMCV64)
            IF(IOPTX.NE.0) THEN
               CALL UNPAK32(ITEST(1),ITEST(1),100)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
*
*     Work out record length
*
            DO 30  J=1, NWORD
*
*     ITEST(J+25) is the pointer in the 2nd record to the
*     file descriptor block, which by definition must be
*     in the same record
*
               IF(ITEST(J+25).GT.0.AND.ITEST(J+25).LE.J) THEN
*
*     Possible record length
*
                  LRC = ITEST(J+ITEST(J+25)+1)
                  IF(LRC.EQ.J) THEN
                     LE  = ITEST(J+30)
                     LD  = ITEST(J+24)
                     NRD = ITEST(J+LD)
*     Does directory size match (record length)*(number of records)?
                     IF(NRD*LRC.NE.LE) GOTO 30
                     LRECL = J
                     IF(IOPTC.EQ.0) THEN
                        CLOSE(LUNIT)
                     ELSE
                        CALL CFCLOS(LUNPTR,0)
                     ENDIF
                     GOTO 40
                  ENDIF
               ENDIF
   30       CONTINUE
            IF(IOPTC.EQ.0) THEN
               CLOSE(LUNIT)
            ELSE
               CALL CFCLOS(LUNPTR,0)
            ENDIF
 
#endif
#if (defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL))&&(!defined(CERNLIB_FQXISN))
            IF(IOPTX.EQ.0.AND.IPASS.EQ.0) THEN
               WRITE(IQPRNT,10400)
10400 FORMAT(' RZOPEN. problems determining record length - ',
     +       ' trying EXCHANGE mode.')
               IOPTX  = 1
               IMODEX = 1
               IPASS  = 1
               GOTO 10
            ENDIF
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
            WRITE(IQPRNT,*) ' RZOPEN. Error in the input file'
            ISTAT = 3
            GOTO 70
#endif
         ENDIF
      ENDIF
   40 CONTINUE
*
#if defined(CERNLIB_QMCV64)||defined(CERNLIB_QMCRY)||defined(CERNLIB_RZFRECL)
      IF(IOPTC.EQ.0) THEN
      NBYTES = IBYTES
#endif
#if defined(CERNLIB_QMCV64)||defined(CERNLIB_QMCRY)
      IF(IMODEX.NE.0) NBYTES = NBYTES / 2
#endif
#if (!defined(CERNLIB_QMIBMFVS))&&(!defined(CERNLIB_QMVAX))
      OPEN(UNIT=LUNIT,FILE=CHFILE,FORM='UNFORMATTED',
     +     RECL=NBYTES*LRECL,ACCESS='DIRECT',STATUS=STAT,IOSTAT=ISTAT)
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(!defined(CERNLIB_IBMMVS))
*
*     CHOPT: SU  = shared,update
*            S   = shared,readonly
*            U   = update
*            ' ' = readonly
*            1   = single user read/write
*
      CHACT = 'READ'
      IF(IOPTN.NE.0.OR.IOPTU.NE.0.OR.IOPT1.NE.0) CHACT = 'READWRITE'
      IF(IOPTN.EQ.0)THEN
         INQUIRE(FILE=CHFILE,EXIST=IEXIST)
         IF(.NOT.IEXIST) GOTO 60
      ENDIF
*
      CALL FILEINF(ISTAT,'MAXREC',2)
      IF(ISTAT.NE.0)GOTO 60
      IF(IOPTN.NE.0)THEN
         OPEN(UNIT=LUNIT,FILE=CHFILE,FORM='UNFORMATTED',RECL=4*LRECL,
     +        ACTION=CHACT,
     +        ACCESS='DIRECT',STATUS=STAT,IOSTAT=ISTAT)
         IF(ISTAT.NE.0)GOTO 60
         CLOSE(LUNIT)
      ENDIF
      CALL FILEINF(ISTAT,'MAXREC',16777215)
      IF(ISTAT.NE.0)GOTO 60
         OPEN(UNIT=LUNIT,FILE=CHFILE,FORM='UNFORMATTED',RECL=4*LRECL,
     +        ACTION=CHACT,
     +        ACCESS='DIRECT',STATUS=STAT,IOSTAT=ISTAT)
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))
*
*     CHOPT: SU  = shared,update
*            S   = shared,readonly
*            U   = update
*            ' ' = readonly
*            1   = single user read/write
*
      CHACT = 'READ'
      IF(IOPTN.NE.0.OR.IOPTU.NE.0.OR.IOPT1.NE.0) CHACT = 'READWRITE'
      IRECL = LRECL * 4
      IF(IRECL.LE.0) IRECL=IDCB(1)
 
      IF(IOPTN.EQ.0)THEN
 
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))&&(!defined(CERNLIB_GSI))
*
*--      open file in read/write mode, read mode only if file
*--           name start with a dot (file belonging to an other userid
*--           from an other group - racf protection)
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))&&(!defined(CERNLIB_NEWLIB))
         IF ( CFNAME(1:1) .EQ. '.'         .AND.
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))&&(defined(CERNLIB_NEWLIB))
         IF ( CFNAME(1:1) .EQ. '.')        THEN
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))&&(!defined(CERNLIB_NEWLIB))
     +        CFNAME(2:5) .NE. PREFIX(1:4))THEN
*SELF,IF=QMIBMFVS,IF=IBMMVS,IF=NEWLIB.
*    +        CFNAME(2:7) .NE. PREFIX(1:6))THEN
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))
            CHACT = 'READ'
         ELSE
            CHACT = 'READWRITE'
         ENDIF
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))
 
         OPEN ( UNIT=LUNIT, STATUS='OLD', FILE=CHFILE, ACCESS='DIRECT'
     +,         FORM='UNFORMATTED', IOSTAT=ISTAT
     +,         RECL=IRECL, ACTION=CHACT)
         IF (ISTAT.NE.0) THEN
            WRITE(IQPRNT,*) 'RZOPEN: OPEN Error file ',CHFILE,
     +                   ' IOSTAT= ', ISTAT
            GOTO 60
         ENDIF
      ENDIF
*
      IF(IOPTN.NE.0)THEN
* detect, whether file really does not exist if opt N given:
         INQUIRE(FILE=CHFILE,EXIST=IEXIST)
         IF(IEXIST) THEN
            WRITE(IQPRNT,*) 'RZOPEN. Warning,input file ',CHFILE
     +       (2:LENOCC(CHFILE)),' already exists !'
            STAT='UNKNOWN'
         ELSE
            STAT='NEW'
*--      define file parameters
            IF( INDEX(CHFILE,'(') .NE. 0 ) THEN
               IF(IOPTV.EQ.0) THEN
                  ISPACE(3) = 28
               ELSE
                  WRITE(IQPRNT,*)
     +             'RZOPEN: No PDS allowed for VSAM files !'
                  ISTAT = -1
                  GOTO 60
               ENDIF
            ENDIF
            CALL FILEINF ( ISTAT, 'DEVICE', IUNIT, ITRK, ISPACE(1)
     +,                 'SECOND', ISPACE(2), 'DIR', ISPACE(3)
     +,                 'RECFM', IRECFM, 'LRECL', IRECL
     +,                 'BLKSIZE', IDCB(2) )
            IF (ISTAT.NE.0) THEN
               WRITE(IQPRNT,*) 'RZOPEN: FILEINF Error, ISTAT = ',ISTAT
               GOTO 60
            ENDIF
*--      If VSAM, define Cluster
            IF(IOPTV.NE.0) THEN
               WRITE(cRecl,'(I6)') iRecl
               cDataClas = 'DATACLAS(RZFILE)'
               IF ( INDEX( CHFILE, '.RZDATA' ) .GT. 0 )
     +       cDataClas = 'DATACLAS(RZDATA)'
               cTSO = 'DEFINE CLUSTER(NAME('''//
     +           CHFILE(2:lenocc(CHFILE)) //
     +          ''') '//cDataClas// ' recordsize(' // cRecl // cRecl //
     +          ') shareoptions(2) NUMBERED UNIQUE NONSPANNED recovery)'
               iStat = mvsams( cTSO )
               IF ( iStat .NE. 0 ) THEN
                  WRITE(IQPRNT,*)
     +             'RZOPEN: Define Cluster Error, ISTAT = ', ISTAT
                  GOTO 60
               ENDIF
            ENDIF
         ENDIF
*--      OPEN FILE
         OPEN ( UNIT=LUNIT, STATUS=STAT
     +,         FILE=CHFILE, ACCESS='DIRECT'
     +,         FORM='UNFORMATTED', IOSTAT=ISTAT
     +,         RECL=IRECL, ACTION='READWRITE' )
      ENDIF
#endif
#if defined(CERNLIB_QMVAX)
*
*     CHOPT: SU  = shared,update
*            S   = shared,readonly
*            U   = update
*            ' ' = readonly
*            1   = single user read/write
*
      IF(IOPTS.NE.0.AND.IOPTU.NE.0)THEN
         OPEN(UNIT=LUNIT,FILE=CFNAME,FORM='UNFORMATTED',RECL=LRECL,
     +        ORGANIZATION=CHORG,
     +        ACCESS='DIRECT',SHARED,STATUS=STAT,IOSTAT=ISTAT)
      ELSEIF(IOPTS.NE.0.AND.IOPTU.EQ.0)THEN
         OPEN(UNIT=LUNIT,FILE=CFNAME,FORM='UNFORMATTED',RECL=LRECL,
     +        ORGANIZATION=CHORG,
     +        ACCESS='DIRECT',SHARED,READONLY,STATUS=STAT,IOSTAT=ISTAT)
      ELSEIF(IOPTS.EQ.0.AND.IOPTU.NE.0)THEN
         OPEN(UNIT=LUNIT,FILE=CFNAME,FORM='UNFORMATTED',RECL=LRECL,
     +        ORGANIZATION=CHORG,
     +        ACCESS='DIRECT',SHARED,STATUS=STAT,IOSTAT=ISTAT)
      ELSEIF(IOPT1.NE.0)THEN
         OPEN(UNIT=LUNIT,FILE=CFNAME,FORM='UNFORMATTED',RECL=LRECL,
     +        ORGANIZATION=CHORG,
     +        ACCESS='DIRECT',STATUS=STAT,IOSTAT=ISTAT)
      ELSEIF(IOPTS.EQ.0.AND.IOPTU.EQ.0.AND.IOPTN.EQ.0)THEN
         OPEN(UNIT=LUNIT,FILE=CFNAME,FORM='UNFORMATTED',RECL=LRECL,
     +        ORGANIZATION=CHORG,
     +        ACCESS='DIRECT',READONLY,STATUS=STAT,IOSTAT=ISTAT)
      ELSEIF(IOPTN.NE.0)THEN
         OPEN(UNIT=LUNIT,FILE=CFNAME,FORM='UNFORMATTED',RECL=LRECL,
     +        ORGANIZATION=CHORG,
     +        ACCESS='DIRECT',SHARED,STATUS=STAT,IOSTAT=ISTAT)
      ENDIF
#endif
#if (defined(CERNLIB_QMIBMFVS))&&(defined(CERNLIB_IBMMVS))
*
* --- for new datasets, do a dummy initial load
*
      IF ((IOPTV.NE.0).AND.(STAT.EQ.'NEW')) THEN
         CALL RZIODO(LUNIT, irecl/4, 1, record, 2 ) ! dummy write
         CALL RZIODO(LUNIT, irecl/4, 1, record, 1 ) ! dummy read
         WRITE(IQPRNT,*)  'RZOPEN: VSAM FILE ', CHFILE(2:LENOCC(CHFILE))
     +,   ' succesfully created'
      ENDIF
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
      ELSE
*     CHOPT: SU  = shared,update
*            S   = shared,readonly
*            U   = update
*            ' ' = readonly
*            1   = single user read/write
         CHOPE = 'r'
         IF(IOPTU.NE.0.OR.IOPT1.NE.0) CHOPE = 'r+'
         IF(IOPTN.NE.0) CHOPE = 'w+'
         JRECL = LRECL
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
         IF(IMODEX.NE.0) JRECL = LRECL / 2
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
         CALL CFOPEN(LUNPTR,0,JRECL,CHOPE,0,CHFILE,ISTAT)
         LUNIT = 1000 + LUNPTR
      ENDIF
#endif
      IF(ISTAT.NE.0)GOTO 60
      IF(IOPTY.NE.0)GOTO 50
      if(kstrip.ne.0)go to 50
*
*           Check consistency of LRECL
*
      IF(IOPTN.EQ.0.AND.IPASS.EQ.0.AND.ISTAT.EQ.0)THEN
         IMODEX=IOPTX
         IZRECL=LRECL
         CALL RZIODO(LUNIT,50,2,ITEST,1)
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
         CALL VXINVB(ITEST(9),1)
         IF(JBIT(ITEST(9),12).NE.0)THEN
            IMODEX=1
            CALL RZIODO(LUNIT,50,2,ITEST,1)
         ELSE
            CALL VXINVB(ITEST(9),1)
         ENDIF
#endif
         LB=ITEST(25)
*
*     Protection against bad files
*
         IF(LB.GT.8187) THEN
            WRITE(IQPRNT,10500) CHFILE(1:LENOCC(CHFILE))
10500       FORMAT(' RZOPEN: cannot determine record length.',
     +             ' File ',A,' probably not in RZ format')
            LRECP=-1
            ISTAT=2
#if defined(CERNLIB_IBM)
            CLOSE(LUNIT)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
            IF(IOPTC.EQ.0) THEN
               CLOSE(LUNIT)
            ELSE
               CALL CFCLOS(LUNIT-1000,0)
            ENDIF
#endif
            GOTO 70
         ENDIF
         IF(LB.GT.48) CALL RZIODO(LUNIT,LB+6,2,ITEST,1)
         LRECP=ITEST(LB+1)
         IQUEST(1)=0
         IF(LRECP.NE.LRECL)THEN
            LRECL2=LRECL
            LRECL=0
#if defined(CERNLIB_IBM)
            CLOSE(LUNIT)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)||defined(CERNLIB_RZFRECL)
            IF(IOPTC.EQ.0) THEN
               CLOSE(LUNIT)
            ELSE
               CALL CFCLOS(LUNIT-1000,0)
            ENDIF
#endif
 
#if defined(CERNLIB_IBMMVS)
            IF(IPASS.EQ.0.AND.LRECI.EQ.0) THEN
#endif
#if !defined(CERNLIB_IBMMVS)
            IF(IPASS.EQ.0) THEN
#endif
               IPASS=1
               GOTO 10
            ELSE
               WRITE(IQPRNT,*) 'Cannot determine record length'
               ISTAT = 1
               GOTO 70
            ENDIF
         ENDIF
      ENDIF
*
      IF (IPASS.NE.0 .AND. LRECL2.NE.0) THEN
         WRITE(IQPRNT,10600) LRECL2,LRECL
10600    FORMAT(' RZOPEN:  LRECL inconsistant - ',
     +          ' file was opened with LRECL = ',I6,
     +          ' should be LRECL = ',I6)
      ENDIF
*
*          If option 'W' build CHDIR
*
   50 IF(IOPTW.NE.0)THEN
         IF (IOPTC .EQ. 0) THEN
            LUN = LUNIT
         ELSE
            LUN = LUNIT - 1000
         ENDIF
         IF(LUN.LT.10)WRITE(CHDIR,10700)LUN
         IF(LUN.GE.10)WRITE(CHDIR,10800)LUN
10700    FORMAT('LUN',I1,'    ')
10800    FORMAT('LUN',I2,'   ')
      ENDIF
*
   60 CONTINUE
      IQUEST(10) = LRECL
      IQUEST(11) = LUNIT
      IQUEST(12) = IMODEX
#if defined(CERNLIB_QMUIX)
*-*  Fill structure if file is striped
      if(kstrip.ne.0.and.istat.eq.0)then
         if(ioptc.eq.0)then
            lun=lunit
         else
            lun=lunit-1000
         endif
         nstrip(lun)=nst
         nrstrip(lun)=nrs
         istrip(lun)=1
         islast=lun
         rznames(lun)=rzsfile(maxstrip)
      endif
#endif
   70 CONTINUE
      END
