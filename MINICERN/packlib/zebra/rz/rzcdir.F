*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:40  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZCDIR(CHPATH,CHOPT)
*
************************************************************************
*
*      To set,read or print the Current Working Directory
* Input:
*   *CHPATH* Character  variable  specifying  the  pathname  of  the  CWD
*            (default option).
*            Unless several RZ files are open at the same time,  the path
*            name can  be specified  either as a  path starting  with the
*            character  '/',   in  which case  an  absolute  pathname  is
*            intended for the given top directory.  When several RZ files
*            are open,   an absolute  pathname must  start with  a double
*            slash '//' and  the top directory.   When  the pathname does
*            not start with a '/', the pathname is prefixed with the path
*            of the CWD.
*            CHPATH = ' ' means the CWD (useful with the 'U' option)
*   CHOPT    Character variable specifying the option
*            'R'   Read the CWD pathname into CHPATH
*            'P'   Print the CWD
*            ' '   Set the CWD
*            'U'   The same  as the  default but  the time  stamp in  the
*                  directory in memory is checked  against the one on the
*                  file and if needed the  directory in memory is brought
*                  up  to date. This option  should be  used when the user
*                  expects that  directories can be  changed concurrently
*                  by another user and he wants to use the latest version
*            'Q'   quiet - do not print message if directory does not exist
*            'K'   Keep current directory in memory. Do not drop
* Output:
*   *CHPATH* Character variable  containing the complete pathname  of the
*            current working directory (with 'R' option only).
*
* Called by <USER>,RZINPA,RZLDIR
*
*  Author  : R.Brun DD/US/PD
*  Written : 02.04.86
*  Last mod: 09.06.93 JDS. Protection against no RZ files open
*
*  IQUEST(1) = 1 : RZIODO error
*              2 : unknown directory
*              3 : directory overwritten
*              4 : no control bank - RZFILE/RZMAKE not called
*              5 : no open files (RZEND called for all)
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
      DIMENSION    IOPTV(5)
      EQUIVALENCE (IOPTR,IOPTV(1)), (IOPTP,IOPTV(2)), (IOPTU,IOPTV(3))
      EQUIVALENCE (IOPTK,IOPTV(4)), (IOPTQ,IOPTV(5))
      CHARACTER*(*) CHPATH,CHOPT
      CHARACTER*1   COPTQ
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
 
      CALL UOPTC (CHOPT,'RPUKQ',IOPTV)
      IF(IOPTK.NE.0) IOPTU=0
      IF(IOPTR.NE.0) CHPATH = ' '
*
*     Check for RZ control bank (i.e. did we call RZFILE/MAKE yet)
*
      IF(LQRS.EQ.0) THEN
         IQUEST(1) = 4
         GOTO 999
      ENDIF
*
*     Check that we have at least one RZ file still open
*
      LRZ=LQRS
   10 IF(LRZ.EQ.0) GOTO 20
      IF(IQ(KQSP+LRZ-5).NE.0) GOTO 30
      LRZ=LQ(KQSP+LRZ)
      GO TO 10
   20 CONTINUE
      IQUEST(1) = 5
      GOTO 999
   30 CONTINUE
*
*             Read Working directory
*
      IF(IOPTR.NE.0)THEN
         CALL RZPAFF(CHCDIR,NLCDIR,CHPATH)
         GO TO 999
      ENDIF
*
*             Print Working directory
*
      IF(IOPTP.NE.0)THEN
         CALL RZPAFF(CHCDIR,NLCDIR,CHL)
         WRITE(IQPRNT,10000)CHL(1:LENOCC(CHL))
10000 FORMAT(' Current Working Directory = ',A)
         GO TO 999
      ENDIF
*
*             Set Working directory
*             Mark old working directory to be dropped
*             Write current directory if modified
*
      COPTQ = ' '
      IF(IOPTQ.NE.0) COPTQ = 'Q'
      IF(LCDIR.NE.0.AND.ISAVE.NE.0.AND.IOPTK.EQ.0)THEN
         LBANK=LCDIR
   40    IF(LBANK.NE.LTOP)THEN
            LUP=LQ(KQSP+LBANK+1)
            IF(IOPTU.EQ.0)THEN
               CALL SBIT1(IQ(KQSP+LBANK),IQDROP)
            ELSE
               CALL MZDROP(JQPDVS,LBANK,' ')
               IQ(KQSP+LTOP+KIRIN)=0
            ENDIF
            LBANK=LUP
            IF(LBANK.NE.0)GO TO 40
         ENDIF
      ENDIF
*
*             Read Top Dir In Case Of U Option
*
      IF(IOPTU.NE.0)CALL RZRTOP
*
      IF(ISAVE.NE.0)THEN
         CALL RZSAVE
      ENDIF
*
*             Set new directory
*
      CALL RZPATH(CHPATH)
      CALL RZFDIR('RZCDIR',LT,LDIR,COPTQ)
      IF(LDIR.NE.0)THEN
         NLCDIR= NLPAT
         LCDIR = LDIR
         LTOP  = LT
         DO 50 I=1,NLPAT
            CHCDIR(I)=CHPAT(I)
   50    CONTINUE
      ELSE
*     Already set by RZFDIR
*        IQUEST(1)=1
         IF(LCDIR.NE.0)CALL SBIT0(IQ(KQSP+LCDIR),IQDROP)
         GO TO 999
      ENDIF
 
      LFREE  = LQ(KQSP+LTOP-2)
      LUSED  = LQ(KQSP+LTOP-3)
      LPURG  = LQ(KQSP+LTOP-5)
      LROUT  = LQ(KQSP+LTOP-6)
      LRIN   = LQ(KQSP+LTOP-7)
      LB     = IQ(KQSP+LTOP+KLB)
      LREC   = IQ(KQSP+LTOP+LB+1)
      LUN    = IQ(KQSP+LTOP-5)
      IZRECL = IQ(KQSP+LTOP+LB+1)
      IMODEC = JBIT(IQ(KQSP+LTOP),5)
      IMODEH = JBIT(IQ(KQSP+LTOP),6)
#if defined(CERNLIB_FQXISN)
*
*     Set exchange mode bit
*
      CALL SBIT1(IQ(KQSP+LTOP+KPW1+2),12)
      IMODEX = 1
#endif
#if !defined(CERNLIB_FQXISN)
      IMODEX = JBIT(IQ(KQSP+LTOP+KPW1+2),12)
#endif
#if defined(CERNLIB_QMVAX)
      IRELAT = JBIT(IQ(KQSP+LTOP),4)
#endif
      IQUEST(7)=IQ(KQSP+LCDIR+KNKEYS)
      IQUEST(8)=IQ(KQSP+LCDIR+KNWKEY)
      IQUEST(9)=IQ(KQSP+LCDIR+KNSD)
      IQUEST(10)=IQ(KQSP+LCDIR+KQUOTA)
      IQUEST(11)=LCDIR
      IQUEST(12)=LTOP
      IQUEST(13)=IQ(KQSP+LCDIR+KLK)
      CALL RZDATE(IQ(KQSP+LCDIR+KDATEC),IDATEC,ITIMEC,1)
      CALL RZDATE(IQ(KQSP+LCDIR+KDATEM),IDATEM,ITIMEM,1)
      IQUEST(14)=IDATEC
      IQUEST(15)=ITIMEC
      IQUEST(16)=IDATEM
      IQUEST(17)=ITIMEM
      IQUEST(18)=IQ(KQSP+LCDIR+KRUSED)
      IQUEST(19)=IQ(KQSP+LCDIR+KMEGA)
      IQUEST(20)=IQ(KQSP+LCDIR+KWUSED)
      IQUEST(21)=IQ(KQSP+LCDIR+IQ(KQSP+LCDIR+KLD))
*
*             Check password
*
      IF(JBYT(IQ(KQSP+LCDIR+KPW1+2),6,5).NE.0)THEN
         IF(IQ(KQSP+LCDIR+KPW1).NE.IHPWD(1).OR.
     +      IQ(KQSP+LCDIR+KPW1+1).NE.IHPWD(2))THEN
            CALL SBIT1(IQ(KQSP+LCDIR),1)
         ELSE
            CALL SBIT0(IQ(KQSP+LCDIR),1)
         ENDIF
      ENDIF
      IF(JBIT(IQ(KQSP+LTOP),1).NE.0)CALL SBIT1(IQ(KQSP+LCDIR),1)
*
  999 END
