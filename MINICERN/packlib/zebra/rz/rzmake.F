*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:01  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZMAKE(LUNIN,CHDIR,NWKEY,CHFORM,CHTAG,NRECP,CHOPT)
*
************************************************************************
*
*           Routine to create a new RZ file
*           To use an already existing file CALL RZFILE
* Input:
*   LUNP    Logical unit number associated with  the RZ file.   A FORTRAN
*           OPEN statement must precede the call to RZFILE.
*           Starting address of the memory area which will contain the RZ
*           information ('M' option)
*   CHDIR   Character variable specifying  the name of the  top directory
*           to be associated with unit LUN.
*   NWKEY   Number of words associated to a key (maximum 5)
*   CHFORM  Character variable describing each element  of the key vector
*           'B' Bit string but not zero
*           'H' Hollerith (4 characters)
*           'I' Integer (nonzero)
*           Ex: CHFORM='IIH' for NWKEY=3 and the 2 first keys are integer
*               and the third one is Hollerith
*   CHTAG   Character array defined as CHARACTER*8 CHTAG(NWKEY).
*           Each  element of  the  array allows  the  description of  the
*           corresponding element in the key vector with a tag of up to 8
*           characters.
*   NRECP   Number of physical records for primary allocation
*   CHOPT   Character variable specifying the selected options.
*           medium
*             default
*                   Disk
*             'M'   Memory
*                   In this  case the user  must have allocated  at least
*                   NRECP*LUNP words of memory starting at address LUN.
*           mode
*             default
*                   Native mode
*             'X'   Exchange mode
*           other
*             'F'   Format NRECP records (unless 'M')
*             'C'   C I/O (unless 'M')
*                   LRECL (words) taken from IQUEST(10)
*             'N'   New format for Cycle information (default is old)
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 14.09.93 No longer force exchange mode for LINUX
*          : 09.03.94 S.Banerjee (Change in cycle structure)
*          : 30.01.95 J.Shiers. Permit nrecp>65000 for new format
*
************************************************************************
*
#include "zebra/zunit.inc"
#include "zebra/zstate.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzcycle.inc"
#if defined(CERNLIB_QMVAX)
      CHARACTER*16 CHORG
#endif
      CHARACTER    CHOPT*(*),CHDIR*(*),CHFORM*(*)
      CHARACTER*16 CHTOP
      CHARACTER*(*)  CHTAG(*)
      DIMENSION    IOPTV(5),IHDIR(2)
      EQUIVALENCE (IOPTM,IOPTV(1)), (IOPTX,IOPTV(2))
     +,           (IOPTF,IOPTV(3)), (IOPTC,IOPTV(4))
     +,           (IOPTN,IOPTV(5))
*
*-----------------------------------------------------------------------
*
      IQUEST(1)=0
      LOGLV = MIN(NQLOGD,4)
      LOGLV = MAX(LOGLV,-3)
      LUNP  = LUNIN
*
      CALL UOPTC(CHOPT,'MXFCN',IOPTV)
#if !defined(CERNLIB_QCFIO)
      IF(IOPTC.NE.0) THEN
         WRITE(IQPRNT,*) 'RZMAKE. option C ignored - valid only ',
     +      'for MSDOS, Unix and VMS systems'
         IOPTC = 0
      ENDIF
#endif
      IMODEX=IOPTX
      IMODEC=IOPTC
      IF(IOPTC.NE.0) LUNP = IQUEST(11)
#if defined(CERNLIB_FQXISN)
      IMODEX=1
#endif
#if defined(CERNLIB_QMLNX)
C     RZfile maked always with Exchange mode for LINUX
*     IMODEX=1
#endif
      IRELAT=0
*
*                Check NWKEY and NRECP
*
      IF(NWKEY.LE.0.OR.NWKEY.GT.KNMAX)THEN
         IF(LOGLV.GE.-2) WRITE(IQLOG,9010)
 9010    FORMAT(' RZMAKE. NWKEY input value is invalid')
         IQUEST(1) =1
         IQUEST(11)=NWKEY
         GO TO 99
      ENDIF
      IF(NRECP.LT.2.OR.(NRECP.GT.65000.AND.IOPTN.EQ.0))THEN
         IF(LOGLV.GE.-2) WRITE(IQLOG,9011)
 9011    FORMAT(' RZMAKE. NRECP input value is invalid')
         IQUEST(1) =1
         IQUEST(11)=NRECP
         GO TO 99
      ENDIF
*
*          Save existing material (if any)
*
      CALL RZSAVE
*
*            Find record length (as specified in the OPEN statement)
*
*          A, Memory option. LUN contains the buffer address
*                            and the value of LUNP is the block length
*
      IF(IOPTM.NE.0)THEN
         LRECP=LUNP
         IF(LRECP.LT.100.OR.LRECP.GT.10000)LRECP=1024
      ELSE
*
*          B, Standard option DISK. Use information as specified
*             in the Fortran OPEN statement
*
          IF(IOPTC.EQ.0) THEN
#if (!defined(CERNLIB_QMVAX))&&(!defined(CERNLIB_QMVDS))
            INQUIRE(UNIT=LUNP,RECL=LRECB)
#endif
#if defined(CERNLIB_QMVDS)
            LRECB=4096
#endif
#if defined(CERNLIB_QMVAX)
            INQUIRE(UNIT=LUNP,RECL=LRECB,ORGANIZATION=CHORG)
            IF(CHORG.EQ.'RELATIVE')IRELAT=1
#endif
*
#if defined(CERNLIB_RZBYTES)
            LRECP=LRECB/4
#endif
#if !defined(CERNLIB_RZBYTES)
            LRECP=LRECB
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
            IF(IOPTX.EQ.0) THEN
               LRECP=LRECB/8
            ELSE
               LRECP=LRECB/4
            ENDIF
#endif
         ELSE
*
*     Take LRECL from IQUEST(10) in case of C I/O option
*
            LRECP = IQUEST(10)
         ENDIF
      ENDIF
*
      LUN = LUNP
      IZRECL = LRECP
      IF(LUN.LE.0.AND.IOPTM.EQ.0)THEN
         IF(LOGLV.GE.-2) WRITE(IQLOG,9012)
 9012    FORMAT(' RZMAKE. LUN input value is invalid')
         IQUEST(1) =1
         IQUEST(11)=LUN
         GO TO 99
      ENDIF
      IF(LRECP.LT.50)THEN
         IF(LOGLV.GE.-2) WRITE(IQLOG,9013)
 9013    FORMAT(' RZMAKE. LRECP input value less than 50')
         IQUEST(1) =1
         IQUEST(11)=LRECP
         GO TO 99
      ENDIF
      IF(LOGLV.GE.0) WRITE(IQLOG,9014) LUNP,LRECP,CHOPT
 9014 FORMAT(' RZMAKE. Unit ',I6,' Initializing with LREC=',I6,
     +', OPT= ',A)
      CALL MZSDIV (0,-7)
*
*           Check if LUN not already defined
*
      LRZ=LQRS
  10  IF(LRZ.NE.0)THEN
         IF(IQ(KQSP+LRZ-5).EQ.LUN)THEN
            IF(LOGLV.GE.-2) WRITE(IQLOG,9015)
 9015       FORMAT(' RZMAKE. Logical unit number already in use')
            IQUEST(1) =1
            IQUEST(11)=LUN
            GO TO 99
         ELSE
            LRZ=LQ(KQSP+LRZ)
            GO TO 10
         ENDIF
      ENDIF
*
*            First call to RZMAKE, create link area
*
      IF(LQRS.EQ.0)THEN
         CALL MZLINK(JQPDVS,'RZCL',LTOP,LTOP,LFROM)
         CALL MZBOOK(JQPDVS,LRZ0,LQRS,1,'RZ0 ',2,2,36,2,0)
         IQ(KQSP+LRZ0-5)=0
         ISAVE = 1
         NHPWD = 0
         CALL VBLANK(IHPWD,2)
      ENDIF
      NCHD  = LEN(CHDIR)
      IF(NCHD.GT.16)NCHD=16
      CHTOP = CHDIR(1:NCHD)
*
*            Create control bank
*
      IDTIME=0
      CALL RZDATE(IDTIME,IDATE,ITIME,2)
      KTAGS = KKDES+(NWKEY-1)/10+1
      NREC  = NRECP
      LREC  = LRECP
      NWREC = (NREC-1)/32 +1
      NW    = 50+NWREC
      NRD   = (NW-1)/LREC +1
      NWL   = NRD*LREC
      LD    = KTAGS+2*NWKEY
      LB    = LD+NRD+1
      LS    = LB+3+NWREC
      LK    = LS
      LF    = LS
*
      CALL MZBOOK (JQPDVS,LTOP,LQRS,1,'RZ  ',10,9,NWL,2,0)
*
*            Disk or memory
*
      IF(IOPTM.EQ.0)THEN
         IQ(KQSP+LTOP-5) = LUN
*
*            C I/O?
         IF(IOPTC.NE.0) CALL SBIT1(IQ(KQSP+LTOP),5)
      ELSE
         NMEM=IQ(KQSP+LRZ0)+1
         IQ(KQSP+LRZ0)=NMEM
         IQ(KQSP+LTOP-5)=-NMEM
         IF(2*NMEM.GT.IQ(KQSP+LRZ0-1))THEN
            CALL MZPUSH(JQPDVS,LRZ0,0,10,'I')
         ENDIF
         IQ(KQSP+LRZ0+2*NMEM-1)=LOCF(LUNP)-LOCF(IQ(1))+1
         IQ(KQSP+LRZ0+2*NMEM  )=LRECP
         LUN=-NMEM
      ENDIF
*
*            Pre-format file
*
      IF((IOPTF.NE.0).AND.(IOPTM.EQ.0))THEN
         DO 100 I=2,NRECP
  100    CALL RZIODO(LUN,LREC,I,IQ(KQSP+LTOP+1),2)
         IF(IQUEST(1).NE.0)THEN
            IF(LOGLV.GE.-1) WRITE(IQLOG,1000) I-1
 1000       FORMAT(' RZMAKE. Could only pre-format',I6,' records')
            IQUEST(1)=0
         ENDIF
      ENDIF
*
*            Write empty record for locks
*
      CALL RZIODO(LUN,LREC,1,IQ(KQSP+LTOP+1),2)
      IF(IQUEST(1).NE.0) GO TO 99
*
*            Build top-directory parameters
*
      CALL SBIT1(IQ(KQSP+LTOP),2)
      CALL VBLANK(IQ(KQSP+LTOP+1),4)
      CALL UCTOH(CHDIR,IQ(KQSP+LTOP+1),4,NCHD)
      CALL ZHTOI(IQ(KQSP+LTOP+1),IQ(KQSP+LTOP+1),4)
#if defined(CERNLIB_QMVAX)
*
*             Set ORGANIZATION type
*
      IF(IRELAT.NE.0)CALL SBIT1(IQ(KQSP+LTOP),4)
#endif
*
      NHPWD = 0
      CALL VBLANK(IHPWD,2)
      CALL UCOPY(IHPWD,IQ(KQSP+LTOP+KPW1),2)
      IQ(KQSP+LTOP+KPW1+2) = NCHD
      IF(IMODEX.GT.0)THEN
         CALL SBIT1(IQ(KQSP+LTOP+KPW1+2),12)
      ENDIF
      IQ(KQSP+LTOP+KDATEC) = IDTIME
      IQ(KQSP+LTOP+KDATEM) = IDTIME
      IQ(KQSP+LTOP+KQUOTA) = NREC
      IQ(KQSP+LTOP+KRUSED) = NRD
      IQ(KQSP+LTOP+KWUSED) = NWL
      IF (IOPTN.NE.0) THEN
        WRITE(IQLOG,7001)
 7001 FORMAT(' RZMAKE. new RZ format selected.',/,
     +       '         This file will not be readable with versions',
     +       ' of RZ prior to release 94B')
        IQ(KQSP+LTOP+KRZVER) = 1
      ELSE
        IQ(KQSP+LTOP+KRZVER) = 0
      ENDIF
      IQ(KQSP+LTOP+KIP1)   = 2
      IQ(KQSP+LTOP+KNFREE) = NWL-LF
      IQ(KQSP+LTOP+KLD)    = LD
      IQ(KQSP+LTOP+KLB)    = LB
      IQ(KQSP+LTOP+KLS)    = LS
      IQ(KQSP+LTOP+KLK)    = LK
      IQ(KQSP+LTOP+KLF)    = LF
      IQ(KQSP+LTOP+KLC)    = NWL+1
      IQ(KQSP+LTOP+KLE)    = NWL
      IQ(KQSP+LTOP+KNWKEY) = NWKEY
      IQ(KQSP+LTOP+LD)     = NRD
      IQ(KQSP+LTOP+LB)     = NWREC
      IQ(KQSP+LTOP+LB+1)   = LREC
      IQ(KQSP+LTOP+LB+2)   = IDTIME
*
      NCHF=LEN(CHFORM)
      NCH =LEN(CHTAG(1))
      IF(NCH.GT.8)NCH=8
      DO 20 I=1,NWKEY
         IF(NCH.LT.8)CALL VBLANK(IHDIR,2)
         CALL UCTOH(CHTAG(I),IHDIR,4,NCH)
         CALL UCOPY(IHDIR,IQ(KQSP+LTOP+KTAGS+2*(I-1)),2)
         IFORM=2
         IF(I.LE.NCHF)THEN
            IF(CHFORM(I:I).EQ.'B')IFORM=1
            IF(CHFORM(I:I).EQ.'H')IFORM=3
            IF(CHFORM(I:I).EQ.'A')IFORM=4
         ENDIF
         IKDES=(I-1)/10
         IKBIT1=3*I-30*IKDES-2
         CALL SBYT(IFORM,IQ(KQSP+LTOP+KKDES+IKDES),IKBIT1,3)
  20  CONTINUE
      CALL ZHTOI(IQ(KQSP+LTOP+KTAGS),IQ(KQSP+LTOP+KTAGS),2*NWKEY)
      DO 30 I=1,NRD
         IQ(KQSP+LTOP+LD+I)=I+1
         CALL SBIT1(IQ(KQSP+LTOP+LB+3),I+1)
  30  CONTINUE
*
*            Store default LOG level
*
      LOGL = LOGLV + 3
      CALL SBYT(LOGL,IQ(KQSP+LTOP),15,3)
      CALL RZVCYC(LTOP)
*
*            Allocate free records
*
      CALL MZBOOK(JQPDVS,LFREE,LTOP,-2,'RZFR',0,0,3,2,0)
      IQ(KQSP+LFREE-5)=LUN
      IQ(KQSP+LFREE+1)=1
      IQ(KQSP+LFREE+2)=NRD+2
      IQ(KQSP+LFREE+3)=NREC
*
*            Allocate space for used records
*
      CALL MZBOOK(JQPDVS,LUSED,LTOP,-3,'RZUS',0,0,21,2,0)
*
      IQ(KQSP+LUSED-5)=LUN
      LRIN  = 0
      LPURG = 0
      LROUT = 0
      LCDIR = LTOP
      NLCDIR= 1
      NLNDIR= 1
      NLPAT = 1
      CHCDIR(1)=CHTOP
      CHNDIR(1)=CHTOP
      IQUEST(1)=0
*
  99  RETURN
      END
