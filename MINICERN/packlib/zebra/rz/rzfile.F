*
* $Id$
*
* $Log$
* Revision 1.5  1998/09/25 09:33:35  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.4  1997/09/02 15:16:11  mclareni
* WINNT corrections
*
* Revision 1.3  1997/03/14 17:21:19  mclareni
* WNT mods
*
* Revision 1.2  1996/04/24 17:26:50  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZFILE(LUNIN,CHDIR,CHOPT)
*
************************************************************************
*
*           Routine to access an already existing RZ file
*           To create a new RZ file routine RZMAKE should be used
* Input:
*   LUNP    Logical unit number associated with  the RZ file.   A FORTRAN
*           OPEN statement must precede the call to RZFILE.
*           Starting address of the memory area which will contain the RZ
*           information ('M' option)
*   CHDIR   Character variable specifying  the name of the  top directory
*           to be associated with unit LUN.
*   CHOPT   Character variable specifying the selected options.
*           medium
*          default
*                   Disk
*             'M'   Memory mode
*                   In this case space has already been allocated starting
*                   at address LUNP.
*           other
*             'U'   UPDATE mode  default is READ mode only
*             'S'   SHARED mode  default is exclusive access
*             '1'   UPDATE mode  with only one user (no LOCKs required)
*             'L'   Show all locked directories and lock-ids
*             'D'   Reset LOCKing word in first record
*             'C'   C I/O
*             'H'   Hook user routine to RZIODO
*             'X'   Exchange mode file
*
*             'B'   Rebuild bit map of file occupancy from file itself
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 07.04.86
*  Last mod: 22.09.94 JDS - include Z=RZCYCLE and call to RZVCYC
*            06.07.95 JDS - return RZ file version in IQUEST(13)
*
************************************************************************
*
#include "zebra/zunit.inc"
#include "zebra/zstate.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
#include "zebra/rzbuff.inc"
#if defined(CERNLIB_QMVAX)
      CHARACTER*16 CHORG
#endif
      CHARACTER    CHOPT*(*),CHDIR*(*)
      CHARACTER*16 CHTOP
      DIMENSION    IOPTV(10)
      EQUIVALENCE (IOPTM,IOPTV(1)), (IOPTU,IOPTV(2))
      EQUIVALENCE (IOPTS,IOPTV(3)), (IOPTL,IOPTV(4))
      EQUIVALENCE (IOPT1,IOPTV(5)), (IOPTD,IOPTV(6))
      EQUIVALENCE (IOPTC,IOPTV(7)), (IOPTX,IOPTV(8))
      EQUIVALENCE (IOPTB,IOPTV(9)), (IOPTH,IOPTV(10))
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      LOGLV = MIN(NQLOGD,4)
      LOGLV = MAX(LOGLV,-3)
      LUNSA = LUN
      LUNP  = LUNIN
*
*          Save existing material (if any)
*
      CALL RZSAVE
*
      CALL UOPTC (CHOPT,'MUSL1DCXBH',IOPTV)
#if !defined(CERNLIB_QCFIO)
*SELF,IF=-QMCRY,IF=-QMVAX,IF=-QMCV64,IF=-QMAPO,IF=-QMAPO9,IF=-QMUIX,IF=-QMDOS.
      IF(IOPTC.NE.0) THEN
         WRITE(IQPRNT,*) 'RZFILE. option C ignored - valid only ',
     +      'for MSDOS, Unix and VMS systems'
         IOPTC = 0
       ENDIF
#endif
      IRELAT=0
      IMODEC=IOPTC
      IMODEH=IOPTH
      IMODEX=IOPTX
#if (defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))
C     RZfile with Exchange mode for NonPPC-LINUX
*     IMODEX=1
#endif
*
*     Take LRECL and LUNP from IQUEST(10-11) in case of C I/O
*
      IF(IOPTC.NE.0) THEN
         LRECP = IQUEST(10)
         LUNP  = IQUEST(11)
      ENDIF
*
*     Take LRECL and address of user routine from IQUEST(10-11)
*     in case of user I/O routine
*
      IF(IOPTH.NE.0) THEN
         LRECP  = IQUEST(10)
         LUN    = IQUEST(11)
         LUSER  = LUNIN
      ENDIF
*
*            Find record length (as specified in the OPEN statement)
*
*          A, Memory option. LUN contains the buffer address
*                            and the value of LUNP is the block length
*
      IF(IOPTM.NE.0)THEN
         LRECP=1024
         LUN=-99
      ELSEIF(IOPTH.EQ.0) THEN
*
*          B, Standard option DISK. Use information as specified
*             in the Fortran OPEN statement
*
#if defined(CERNLIB_QMVAX)
         IF(IOPTC.EQ.0) THEN
            INQUIRE(UNIT=LUNP,ORGANIZATION=CHORG)
            IF(CHORG.EQ.'RELATIVE')IRELAT=1
         ENDIF
#endif
*
         IZRECL=LRECP
         CALL RZIODO(LUNP,50,2,ITEST,1)
*
*      If option X not specified, determine mode (eXchange, native)
*      from file
*
         IF(IOPTX.EQ.0) THEN
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
            CALL VXINVB(ITEST(9),1)
#endif
            IF(JBIT(ITEST(9),12).NE.0)THEN
               IMODEX=1
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
               CALL RZIODO(LUNP,50,2,ITEST,1)
            ELSE
               CALL VXINVB(ITEST(9),1)
#endif
            ENDIF
         ENDIF
 
         IF(IQUEST(1).NE.0)GO TO 30
         LB=ITEST(KLB)
         IF(LB.GT.48)CALL RZIODO(LUNP,LB+6,2,ITEST,1)
         IF(LB.GT.100)THEN
            IF(LOGLV.GE.-1) WRITE(IQLOG,10000)
10000       FORMAT(' RZFILE. WARNING!! Top directory is big')
         ENDIF
         LRECP=ITEST(LB+1)
#if defined(CERNLIB_QMALT)||defined(CERNLIB_QMAPO)||defined(CERNLIB_QMDOS)||defined(CERNLIB_WINNT)
         IF(IOPTC.EQ.0.AND.IOPTH.EQ.0) THEN
            INQUIRE(UNIT=LUNP,RECL=LRECL)
* DEC Fortran takes "longword" units
#endif
#if (defined(CERNLIB_QFDEC))&&(defined(CERNLIB_QMDOS))
            IF(LRECP.NE.LRECL)THEN
#endif
#if (defined(CERNLIB_QMALT)||defined(CERNLIB_QMAPO)||defined(CERNLIB_QMDOS)||defined(CERNLIB_WINNT))&&(!defined(CERNLIB_QFDEC))
            IF(LRECP.NE.LRECL/4)THEN
#endif
#if defined(CERNLIB_QMALT)||defined(CERNLIB_QMAPO)||defined(CERNLIB_QMDOS)||defined(CERNLIB_WINNT) 
               IQUEST(1)=1
               IF(LOGLV.GE.-2) WRITE(IQLOG,10100)LUNP,LRECP,LRECL/4
10100          FORMAT(' RZFILE. Unit ',I6,'RECL on file ',I5,
     +                ' incompatible with RECL in OPEN =',I5)
               GO TO 30
            ENDIF
         ENDIF
#endif
         LUN=LUNP
         IQUEST(1)=0
      ENDIF
*
      IF(LOGLV.GE.0) WRITE(IQLOG,10200) LUN,LRECP,CHOPT
10200 FORMAT(' RZFILE. UNIT ',I6,' Initializing with LREC=',I6,
     +', OPT= ',A)
      CALL MZSDIV (0,-7)
*
*           Check if LUN not already defined
*
      LRZ=LQRS
   10 IF(LRZ.NE.0)THEN
         IF(IQ(KQSP+LRZ-5).EQ.LUN)THEN
            IQUEST(1)=1
            IF(LOGLV.GE.-2) WRITE(IQLOG,10300)
10300       FORMAT(' RZFILE. Unit is already in use')
            LUN=LUNSA
            GO TO 30
         ELSE
            LRZ=LQ(KQSP+LRZ)
            GO TO 10
         ENDIF
      ENDIF
*
*            First call to RZFILE, create link area
*
      IF(LQRS.EQ.0)THEN
         CALL MZLINK(JQPDVS,'RZCL',LTOP,LTOP,LFROM)
         CALL MZBOOK (JQPDVS,LRZ0,LQRS,1,'RZ0 ',2,2,36,2,0)
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
      CALL MZBOOK(JQPDVS,LTOP,LQRS,1,'RZ  ',10,9,LRECP,2,0)
*
*            Disk or memory
*
      IF(IOPTM.EQ.0)THEN
         IQ(KQSP+LTOP-5) = LUN
*
*            C I/O?
*
         IF(IOPTC.NE.0) CALL SBIT1(IQ(KQSP+LTOP),5)
*
*            user I/O?
*
         IF(IOPTH.NE.0) THEN
            CALL SBIT1(IQ(KQSP+LTOP),6)
            CALL SBYT(LUSER,IQ(KQSP+LTOP),7,7)
         ENDIF
      ELSE
         NMEM=IQ(KQSP+LRZ0)+1
         IQ(KQSP+LRZ0)=NMEM
         IQ(KQSP+LTOP-5)=-NMEM
         IF(2*NMEM.GT.IQ(KQSP+LRZ0-1))THEN
            CALL MZPUSH(JQPDVS,LRZ0,0,10,' ')
         ENDIF
         IQ(KQSP+LRZ0+2*NMEM-1)=LOCF(LUNP)-LOCF(IQ(1))+1
         IQ(KQSP+LRZ0+2*NMEM  )=LRECP
         LUN=-NMEM
      ENDIF
*
*            Read 1st record of directory
*
      CALL RZIODO(LUN,LRECP,2,IQ(KQSP+LTOP+1),1)
      IF(IQUEST(1).NE.0)GO TO 30
      LD   = IQ(KQSP+LTOP+KLD)
      LB   = IQ(KQSP+LTOP+KLB)
      LREC = IQ(KQSP+LTOP+LB+1)
      NRD  = IQ(KQSP+LTOP+LD)
#if defined(CERNLIB_FQXISN)
*
*     Set exchange mode bit
*
      CALL SBIT1(IQ(KQSP+LTOP+KPW1+2),12)
#endif
      IMODEX=JBIT(IQ(KQSP+LTOP+KPW1+2),12)
*
*            Increase size of control bank if required
*            and read all records for top directory
*
      NPUSH=NRD*LREC-LRECP
      IF(NPUSH.NE.0)CALL MZPUSH(JQPDVS,LTOP,0,NPUSH,'I')
      DO 20 I=2,NRD
         CALL RZIODO(LUN,LREC,IQ(KQSP+LTOP+LD+I),
     +               IQ(KQSP+LTOP+(I-1)*LREC+1),1)
         IF(IQUEST(1).NE.0)GO TO 30
   20 CONTINUE
      CALL VBLANK(IQ(KQSP+LTOP+1),4)
      CALL UCTOH(CHDIR,IQ(KQSP+LTOP+1),4,NCHD)
      CALL ZHTOI(IQ(KQSP+LTOP+1),IQ(KQSP+LTOP+1),4)
      CALL SBYT(NCHD,IQ(KQSP+LTOP+KPW1+2),1,5)
      CALL UCOPY(IQ(KQSP+LTOP+KPW1),IHPWD,2)
      NHPWD=JBYT(IQ(KQSP+LTOP+KPW1+2),6,5)
      IQ(KQSP+LTOP+KIRIN)=0
      IQ(KQSP+LTOP+KIROUT)=0
#if defined(CERNLIB_NOTNEW)
*
*     Check that the file is not in the NEW format
*
      IF (IQ(KQSP+LTOP+KRZVER).NE.0) THEN
         CALL ZFATAM
     +      (' RZFILE. file cannot be processed by this version of RZ')
      ENDIF
#endif
 
      LFREE = 0
      LUSED = 0
      LRIN  = 0
      LPURG = 0
      LROUT = 0
      LCDIR = LTOP
      NLCDIR= 1
      NLNDIR= 1
      NLPAT = 1
      CHCDIR(1)=CHTOP
      CHNDIR(1)=CHTOP
*
*             Reset LOCKing word in record 1
*
      IF(IOPTD.NE.0)THEN
         CALL RZDLOK
      ENDIF
*
*             Show locks
*
      IF(IOPTL.NE.0)THEN
         CALL RZLLOK
      ENDIF
#if defined(CERNLIB_QMVAX)
*
*             Set ORGANIZATION type
*
      IF(IRELAT.NE.0)THEN
         UNLOCK(UNIT=LUN)
         CALL SBIT1(IQ(KQSP+LTOP),4)
      ENDIF
#endif
*
*            Store default LOG level
*
      LOGL = LOGLV + 3
      CALL SBYT(LOGL,IQ(KQSP+LTOP),15,3)
*
*     RZ version
*
      CALL RZVCYC(LTOP)
      IQUEST(13) = IQ(KQSP+LTOP+KRZVER)
*
*     Rebuild bit map?
*
      IF(IOPTB.NE.0) CALL RZVERI('//'//CHTOP(1:NCHD),'B')
*
*             UPDATE mode only
*
      CALL SBIT1(IQ(KQSP+LTOP),1)
      IF(IOPTU.NE.0.OR.IOPT1.NE.0)THEN
*
*            Allocate free records
*
         CALL SBIT0(IQ(KQSP+LTOP),1)
         CALL MZBOOK(JQPDVS,LFREE,LTOP,-2,'RZFR',0,0,21,2,0)
         IQ(KQSP+LFREE-5)=LUN
*
*        IF(IOPTU.EQ.0.AND.IOPT1.EQ.0)THEN
         IF(IOPTS.EQ.0)THEN
            CALL SBIT1(IQ(KQSP+LTOP),3)
            CALL RZLOCK('RZFILE')
            IF(IQUEST(1).NE.0)THEN
               CALL SBIT1(IQ(KQSP+LTOP),1)
               IQ1=IQUEST(1)
               CALL MZDROP(JQPDVS,LFREE,' ')
               LFREE=0
               IQUEST(1)=2+IQ1
               GO TO 30
            ENDIF
         ELSE
            CALL SBIT0(IQ(KQSP+LTOP),3)
         ENDIF
*
*            Allocate space for used records
*
         CALL MZBOOK(JQPDVS,LUSED,LTOP,-3,'RZUS',0,0,21,2,0)
         IQ(KQSP+LUSED-5)=LUN
      ENDIF
      IQUEST(7)=IQ(KQSP+LCDIR+KNKEYS)
      IQUEST(8)=IQ(KQSP+LCDIR+KNWKEY)
*
   30 RETURN
      END
