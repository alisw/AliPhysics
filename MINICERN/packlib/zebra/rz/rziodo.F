*
* $Id$
*
* $Log$
* Revision 1.6  1998/09/25 09:33:38  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.5  1998/03/13 16:51:18  mclareni
* Put the record length warning behind a CERNLIB_DEBUG flag to avoid spurious warnings
*
* Revision 1.4  1997/09/02 15:16:12  mclareni
* WINNT corrections
*
* Revision 1.3  1997/05/14 08:33:39  couet
* - Bug fixed by S.O'Neale. atlas problems with cernlib 97a, with rfio/cio
* the record was not correct in rziodo. Now rzfdir.F rest the correct one.
*
* Revision 1.2  1996/04/24 17:26:56  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZIODO(LUNRZ,JREC,IREC1,IBUF,IRW)
*
**********************************************************************
*
*        To read/write a block from disk or memory
*
*       LUNRZ Logical unit number for disk (LUN>0)
*       JREC  Record length
*       IREC  RECORD NUMBER
*       IBUF  ADDRESS OF BUFFER TO BE READ/WRITTEN ON IREC
*       IRW   =1 TO READ
*       IRW   =2 TO WRITE
*
* Called by RZCOP1,RZDELT,RZDFIR,RZFILE,RZFREE,RZLLOK,RZLOCK
*           RZMAKE,RZMDIR,RZOUT,RZREAD,RZSAVE,RZOPEN
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 10.06.94 Implement File striping for PIAF (R.Brun)
*                     New routine RZSTRIR called
*  Last mod: 26.10.93 IQUEST(1) = 101 in case of READ error,
*                                 102 in case of WRITE error
*
************************************************************************
*
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzbuff.inc"
#include "zebra/rzcount.inc"
#if defined(CERNLIB_QMUIX)
#include "zebra/rzcstr.inc"
#endif
      DIMENSION IBUF(JREC)
      PARAMETER (MEDIUM=0)
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbyt.inc"

*
*     I/O statistics
*
      RZXIO(IRW) = RZXIO(IRW) + JREC
 
      IREC=IREC1
      IF(LUNRZ.GT.0)THEN
         NERR=0
#if defined(CERNLIB_DEBUG)
C    Too many false alarms occur if we leave this in production code.
         IF ( IREC.GT.2 .AND. IZRECL.NE.JREC ) THEN
            WRITE(IQLOG,10010) LUNRZ,IREC,IZRECL,JREC
10010       FORMAT(' RZIODO. Potential Record Length Problem   LUN =',I6,
     +             ' IREC =',I6,' IZRECL =',I6,'  JREC =',I6)
         END IF
#endif
         IF(IMODEH.NE.0) THEN
*
            IQUEST(1) = JBYT(IQ(KQSP+LTOP),7,7)
            IQUEST(2) = JREC
            IQUEST(4) = IREC
            IOWAY     = IRW - 1
*
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
*
*     Byte swap
*
            IF(IRW.EQ.2.AND.IMODEX.NE.0) CALL VXINVB(IBUF,JREC)
#endif
            CALL JUMPST(LUNRZ)
            CALL JUMPX2(IBUF,IOWAY)
 
            IF(IQUEST(1).NE.0) IQUEST(1) = 100 + IRW
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
*
*     Byte swap
*
            IF(IMODEX.NE.0) CALL VXINVB(IBUF,JREC)
#endif
         ELSE
*
*     Read a record
*
   10       IF (IRW.EQ.1)THEN
*
*     Fortran I/O
*
               IF(IMODEC.EQ.0) THEN
#if defined(CERNLIB_QMUIX)
*-* Case of File striping
*
      if(nstrip(lunrz).gt.0)then
         call rzstrir(lunrz,irec)
      endif
*
#endif
#if defined(CERNLIB_QMVDS)
                  CALL READ_ASS(LUNRZ,IREC,ISTAT,IBUF)
                  IF(ISTAT.EQ.1)GO TO 20
#endif
#if (!defined(CERNLIB_QMVDS))&&(!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCV64))
                  READ (UNIT=LUNRZ,REC=IREC,ERR=20,IOSTAT=ISTAT)IBUF
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
*
*     Exchange mode on Cray and Convex - read JREC/2 64 bit
*     words (=JREC 32 bit words)
*
                  IF(IMODEX.NE.0) THEN
                     CALL RZIPHD(LUNRZ,JREC/2,IREC,ITEST,ISTAT)
                     IF(ISTAT.NE.0) GOTO 20
#endif
#if defined(CERNLIB_QMCV64)
*
*     Unpack 32 bit words to 64 bit words
*
                     CALL UPAK32(ITEST,IBUF,JREC)
#endif
#if defined(CERNLIB_QMCRY)
*
*     Unpack 32 bit words to 64 bit words
*
                     CALL UNPACK(ITEST,32,IBUF,JREC)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
*
*     Native mode
*
                  ELSE
                     READ(UNIT=LUNRZ,REC=IREC,ERR=20,IOSTAT=ISTAT)IBUF
                  ENDIF
#endif
*
*     C I/O
*
               ELSE
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
                  IF(IMODEX.NE.0) THEN
*
*      Exchange mode?
*
                     CALL CFSEEK(LUNRZ-1000,MEDIUM,IZRECL/2,IREC-1,
     +               ISTAT)
                     NWTAK = JREC / 2
                     CALL CFGET(LUNRZ-1000,MEDIUM,JREC,NWTAK,ITEST,
     +               ISTAT)
                     IF(ISTAT.NE.0) GOTO 20
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
*
*     Unpack
*
#endif
#if defined(CERNLIB_QMCRY)
                     CALL UNPACK(ITEST,32,IBUF,JREC)
#endif
#if defined(CERNLIB_QMCV64)
                     CALL UPAK32(ITEST,IBUF,JREC)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
                  ELSE
#endif
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMVDS)
*
*     C I/O not supported
*
                     ISTAT = 99
#endif
#if defined(CERNLIB_QMUIX)
*-* Case of File striping
*
      if(nstrip(lunrz-1000).gt.0)then
         call rzstrir(lunrz-1000,irec)
      endif
*
#endif
#if (!defined(CERNLIB_QMVDS))&&(!defined(CERNLIB_QMIBM))
                     CALL CFSEEK(LUNRZ-1000,MEDIUM,IZRECL,IREC-1,ISTAT)
                     NWTAK = JREC
                     CALL CFGET(LUNRZ-1000,MEDIUM,JREC,NWTAK,IBUF,
     +               ISTAT)
                     IF(ISTAT.NE.0) GOTO 20
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
                  ENDIF
#endif
               ENDIF
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
*
*     Byte swap if exchange mode
*
               IF(IMODEX.NE.0) CALL VXINVB(IBUF,JREC)
#endif
*
*     Write
*
            ELSE
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
*
*     Byte swap if exchange mode
*
               IF(IMODEX.NE.0) CALL VXINVB(IBUF,JREC)
#endif
               IF(IMODEC.EQ.0) THEN
#if defined(CERNLIB_QMVDS)
                  CALL WRITE_ASS(LUNRZ,IREC,ISTAT,IBUF)
                  IF(ISTAT.EQ.1)GO TO 20
#endif
#if (!defined(CERNLIB_QMVDS))&&(!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCV64))
                  WRITE(UNIT=LUNRZ,REC=IREC,ERR=20,IOSTAT=ISTAT)IBUF
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
                  IF(IMODEX.NE.0) THEN
*
*     Pack
*
#endif
#if defined(CERNLIB_QMCRY)
                     CALL PACK(IBUF,32,ITEST,JREC)
#endif
#if defined(CERNLIB_QMCV64)
                     CALL PACK32(IBUF,ITEST,JREC)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
                     CALL RZOPHD(LUNRZ,IZRECL/2,IREC,ITEST,ISTAT)
                     IF(ISTAT.NE.0) GOTO 20
                  ELSE
                     WRITE(UNIT=LUNRZ,REC=IREC,ERR=20,IOSTAT=ISTAT)
     +               IBUF
                  ENDIF
#endif
               ELSE
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
                  IF(IMODEX.NE.0) THEN
*
*     Pack
*
#endif
#if defined(CERNLIB_QMCRY)
                     CALL PACK(IBUF,32,ITEST,JREC)
#endif
#if defined(CERNLIB_QMCV64)
                     CALL PACK32(IBUF,ITEST,JREC)
#endif
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCV64)
                     JREC = JREC / 2
                     CALL CFSEEK(LUNRZ-1000,MEDIUM,IZRECL/2,IREC-1,
     +               ISTAT)
                     IF(ISTAT.NE.0) GOTO 20
                     CALL CFPUT(LUNRZ-1000,MEDIUM,JREC/2,ITEST,ISTAT)
                     IF(ISTAT.NE.0) GOTO 20
                  ELSE
                     CALL CFSEEK(LUNRZ-1000,MEDIUM,IZRECL,IREC-1,ISTAT)
                     IF(ISTAT.NE.0) GOTO 20
                     CALL CFPUT(LUNRZ-1000,MEDIUM,JREC,IBUF,ISTAT)
                     IF(ISTAT.NE.0) GOTO 20
                  ENDIF
#endif
#if defined(CERNLIB_QMIBM)||defined(CERNLIB_QMVDS)
*
*     C I/O not supported
*
                  ISTAT = 99
#endif
#if (!defined(CERNLIB_QMVDS))&&(!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QMIBM))
                  CALL CFSEEK(LUNRZ-1000,MEDIUM,IZRECL,IREC-1,ISTAT)
                  IF(ISTAT.NE.0) GOTO 20
                  CALL CFPUT(LUNRZ-1000,MEDIUM,JREC,IBUF,ISTAT)
                  IF(ISTAT.NE.0) GOTO 20
#endif
               ENDIF
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMVMI)||defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_WINNT)
*
*     Byte swap back
*
               IF(IMODEX.NE.0) CALL VXINVB(IBUF,JREC)
#endif
            ENDIF
            RETURN
   20       NERR=NERR+1
#if defined(CERNLIB_QMVAX)
            CALL LIB$WAIT(0.1)
#endif
#if defined(CERNLIB_QMCRY)
            IC = SLEEPF(1)
#endif
            IF(NERR.LT.100)GO TO 10
            IQUEST(1)=100+IRW
            WRITE(IQLOG,1000)IREC,LUNRZ,ISTAT
 1000       FORMAT(' RZIODO. Error at record =',I5,' LUN =',I6,
     +             ' IOSTAT =',I6)
         ENDIF
      ELSE
         KOF=IQ(KQSP+LRZ0-2*LUNRZ-1)+IQ(KQSP+LRZ0-2*LUNRZ)*(IREC-1)
         IF (IRW.EQ.1)THEN
            CALL UCOPY(IQ(KOF),IBUF,JREC)
         ELSE
            CALL UCOPY(IBUF,IQ(KOF),JREC)
         ENDIF
      ENDIF
*
      END
