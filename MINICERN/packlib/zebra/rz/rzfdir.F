*
* $Id$
*
* $Log$
* Revision 1.3  1997/05/14 08:33:37  couet
* - Bug fixed by S.O'Neale. atlas problems with cernlib 97a, with rfio/cio
* the record was not correct in rziodo. Now rzfdir.F rest the correct one.
*
* Revision 1.2  1996/04/24 17:26:48  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:22  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZFDIR(CHROUT,LT,LDIR,CHOPT)
*
************************************************************************
*
*        Check if Pathname stored in CHPAT is a valid directory
*        If YES then LT points to the TOP directory
*                    LDIR points to the directory
*
* Called by RZCDIR,RZCOPY
*
*  Author  : R.Brun DD/US/PD
*  Written : 03.04.86
*  Last mod: 13.05.97 S.O'Neale  Store IZRECL, bug fixes for messages
*          : 02.06.93 JDS. Bug fix in 'unknown directory' warning
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
#include "zebra/rzcycle.inc"
#include "zebra/rzclun.inc"
      CHARACTER*(*) CHROUT
      CHARACTER*(*) CHOPT
      DIMENSION IHDIR(4)
      LOGICAL RZSAME
      INTEGER FQUOTA
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      IOPTQ = INDEX(CHOPT,'Q')
      LT=0
      LDIR=0
      IF(LQRS.EQ.0) GOTO 110
      IF(NLPAT.LE.0)THEN
         CHL='??? '
         GOTO 90
      ENDIF
*
*            Find LT
*
      CALL VBLANK(IHDIR,4)
      CALL UCTOH(CHPAT(1),IHDIR,4,16)
      CALL ZHTOI(IHDIR,IHDIR,4)
      LRZ=LQRS
   10 IF(.NOT.RZSAME(IHDIR,IQ(KQSP+LRZ+1),4))THEN
         LRZ  = LQ(KQSP+LRZ)
         IF(LRZ.GT.0)GOTO 10
         GOTO 80
      ENDIF
      LTEMP  = LRZ
      LT     = LRZ
      LDIR   = LRZ
      CALL RZVCYC(LT)
      IF(NLPAT.LT.2)GOTO 110
      LBT    = IQ(KQSP+LRZ+KLB)
      LREF   = IQ(KQSP+LRZ+LBT+1)
      LUNF   = IQ(KQSP+LRZ-5)
      FQUOTA = IQ(KQSP+LRZ+KQUOTA)
      LOGLV  = JBYT(IQ(KQSP+LT),15,3)-3
      IZRECL = IQ(KQSP+LT+LBT+1)      ! SWON: Needed by CFSEEK in RZIODO
      IMODEX = JBIT(IQ(KQSP+LT+KPW1+2),12)
      IMODEC = JBIT(IQ(KQSP+LT),5)
      IMODEH = JBIT(IQ(KQSP+LT),6)
*
*          Search levels down
*
      DO 60 IL=2,NLPAT
         CALL VBLANK(IHDIR,4)
         CALL UCTOH(CHPAT(IL),IHDIR,4,16)
         CALL ZHTOI(IHDIR,IHDIR,4)
         CALL SBIT0(IQ(KQSP+LRZ),IQDROP)
         NSDIR=IQ(KQSP+LRZ+KNSD)
         LS   =IQ(KQSP+LRZ+KLS)
         IF(NSDIR.LE.0)GOTO 80
*
*              Check if element in list of subdirectories
*
         DO 50 I=1,NSDIR
            IH=LS+7*(I-1)
            IF(RZSAME(IHDIR,IQ(KQSP+LRZ+IH),4))THEN
               IF (KVSCYC.EQ.0) THEN
                  IRS  = JBYT(IQ(KQSP+LRZ+IH+5),1,18)
               ELSE
                  IRS  = IQ(KQSP+LRZ+IH+5)
               ENDIF
*
*    Record number of this subdirectory < 0 or > file quota
*
               IQUEST(20) = 0
               IF(IRS.LE.0.OR.IRS.GT.FQUOTA) GOTO 100
               LRN  = LQ(KQSP+LRZ-1)
   20          IF(LRN.EQ.0)THEN
                  CALL MZBOOK(JQPDVS,LDIR,LRZ,-1,'RZ  ',6,6,LREF,2,-1)
                  LRZ=LDIR
                  CALL RZIODO(LUNF,LREF,IRS,IQ(KQSP+LRZ+1),1)
                  IF(IQUEST(1).NE.0) GOTO 70
                  LDS=IQ(KQSP+LRZ+KLD)
                  IF(LDS.GT.IQ(KQSP+LRZ-1)) GOTO 100
                  IF(LDS.LE.0) GOTO 100
                  NRDS=IQ(KQSP+LRZ+LDS)
                  IF(NRDS.GT.1)THEN
                     CALL MZPUSH(JQPDVS,LRZ,0,LREF*(NRDS-1),' ')
                     LDIR=LRZ
*
*      Number of records, record numbers
*
                     IQUEST(20) = NRDS
                     IQUEST(21) = IRS
                     DO 30 IR=2,NRDS
                        IRS=IQ(KQSP+LRZ+LDS+IR)
                        JR = 20 + IR
                        IF(JR.LE.100) IQUEST(JR) = IRS
                        IF(IRS.LE.0.OR.IRS.GT.FQUOTA) GOTO 100
                        CALL RZIODO(LUNF,LREF,IRS,
     +                              IQ(KQSP+LRZ+(IR-1)*LREF+1),1)
                        IF(IQUEST(1).NE.0)GOTO 70
   30                CONTINUE
                  ENDIF
               ELSE
   40             IF(RZSAME(IHDIR,IQ(KQSP+LRN+1),4))THEN
                     LRZ = LRN
                     LDIR= LRN
                     GOTO 60
                  ELSE
                     LRN=LQ(KQSP+LRN)
                     GOTO 20
                  ENDIF
               ENDIF
               GOTO 60
            ENDIF
   50    CONTINUE
         GOTO 80
   60 CONTINUE
      CALL SBIT0(IQ(KQSP+LDIR),IQDROP)
      LT=LTEMP
#if defined(CERNLIB_QMVAX)
      IF(IRELAT.NE.0)UNLOCK(UNIT=LUNF)
#endif
      GOTO 110
*
*          Errors
*
   70 CONTINUE
*
*     RZIODO error
*
      LDIR = 0
      IQUEST(1) = 1
#if defined(CERNLIB_QMVAX)
      IF(IRELAT.NE.0)UNLOCK(UNIT=LUNF)
#endif
      GOTO 110
   80 CALL RZPAFF(CHPAT,NLPAT,CHL)
#if defined(CERNLIB_QMVAX)
      IF(IRELAT.NE.0)UNLOCK(UNIT=LUNF)
#endif
   90 LDIR=0
      IQUEST(1) = 2       ! SWON: Write a message if "Unknown directory"
*     IF(LOGLV.GE.-2.AND.
*    +   IQUEST(1).EQ.0.AND.IOPTQ.EQ.0)THEN  (retain original code )
      IF(LOGLV.GE.-2.AND.IOPTQ.EQ.0)THEN
         WRITE(IQLOG,10000)CHROUT,CHL(1:LENOCC(CHL))
10000    FORMAT(1X,A,'. Unknown directory ',A)
      ENDIF
      GOTO 110
*
*         Directory overwritten
*
  100 CALL RZPAFF(CHPAT,NLPAT,CHL)
      IQUEST(1) = 3
#if defined(CERNLIB_QMVAX)
      IF(IRELAT.NE.0)UNLOCK(UNIT=LUNF)
#endif
      LDIR=0
      IF(LOGLV.GE.-2)THEN    ! SWON: Write a message if RZ is in trouble
         WRITE(IQLOG,10100)CHROUT,CHL(1:LENOCC(CHL))
10100    FORMAT(1X,A,'. Directory overwritten ',A)
      ENDIF
*
  110 RETURN
      END
