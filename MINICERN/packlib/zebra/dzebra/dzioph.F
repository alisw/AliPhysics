*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:11  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:06  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZIOPH (IFIRST,ILAST)
      SAVE KFOTYP
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzioc.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/bkfoparq.inc"

      CHARACTER CHROUT*(*),CHSTAK*6,KFOTYP(0:11)*1
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
      CHARACTER CDWORD*36
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
      CHARACTER CDWORD*20
#endif
      DOUBLE PRECISION DWORD
      INTEGER                IWORD
      REAL                         RWORD
      EQUIVALENCE     (DWORD,IWORD,RWORD)
      PARAMETER (CHROUT = 'DZIOPH')
      DATA KFOTYP /'U','B','I','F','D','H','*','S','*','N','*','L'/

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif







      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      CALL ZPAGE(IQPRNT,MIN((ILAST-IFIRST)/10+2,5))
      CQLINE = ' --------  DATA part of bank  --------'
      IF (IFLOPT(MPOSQQ).EQ.0) CALL DZTEXT(0,CDUMMQ,1)

      WRITE (CQLINE,'(1X,I7,'' /'',120X)') IFIRST
      INLINE =11

      JDATA  = 0
      JFOCUR = 0
      NWPRNT = 0


   10 ITYPE  = MFO(JFOCUR+1)
      IF (ITYPE.EQ.7)                                      GO TO 40
      NWSEC = MFO(JFOCUR+2)
      IF   (NWSEC)           20, 30, 60


   20 NWSEC = ILAST - JDATA
                                                           GO TO 60


   30 JDATA = JDATA + 1
      IWORD = IQ(KQS+LS+JDATA)
      NWSEC = IWORD
      IF (JDATA.GE.IFIRST) THEN
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
          WRITE(CQLINE(INLINE:INLINE+19),'('' *'',A1,''('',I15,'')'')')
     +     KFOTYP(MIN(ITYPE,6)),NWSEC
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
          WRITE(CQLINE(INLINE:INLINE+11),'('' *'',A1,''('', I7,'')'')')
     +     KFOTYP(MIN(ITYPE,6)),NWSEC
#endif
          NWPRNT = NWPRNT + 1
      ENDIF
                                                           GO TO 50


   40 JDATA = JDATA + 1
      IWORD = IQ(KQS+LS+JDATA)
      ITYPE = MOD (IWORD,16)
      NWSEC = IWORD/16
      IF (JDATA.GE.IFIRST) THEN
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
         WRITE(CQLINE(INLINE:INLINE+19),'('' *S('',A1,'','',I13,'')'')')
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
         WRITE(CQLINE(INLINE:INLINE+11),'('' *S('',A1,'','', I5,'')'')')
#endif
     +     KFOTYP(MIN(ITYPE,6)),NWSEC
          NWPRNT = NWPRNT + 1
      ENDIF

   50 IF (JDATA.GE.IFIRST) THEN
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
          IF (MOD(NWPRNT,6).EQ.0)  THEN
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
          IF (MOD(NWPRNT,10).EQ.0)  THEN
#endif
              IF (IFLOPT(MPOSQQ).NE.0)                    GO TO 999
              CALL DZTEXT(0,CDUMMQ,1)
              WRITE (CQLINE,'(1X,I7,'' /'',120X)') NWPRNT+IFIRST
              INLINE =11
          ELSE
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
              INLINE = INLINE + 20
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
              INLINE = INLINE + 12
#endif
          ENDIF
      ENDIF

      IF (ITYPE.GE.8) THEN
          ITYPE      = 0
          NWSEC      = ILAST - JDATA
          IQUEST(11) = -ITYPE
      ELSEIF (NWSEC.EQ.0) THEN
          ITYPE      = 0
                                                           GO TO 900
      ELSEIF (NWSEC.LT.0) THEN
          ITYPE      = 0
          NWSEC      = ILAST - JDATA
          IQUEST(11) = -ITYPE-32
      ENDIF


   60 IDBLE  = 0
      DO 100 I=MAX(JDATA+1,IFIRST),MIN(JDATA+NWSEC,ILAST)
          IWORD = IQ(KQS+LS+I)
          IF (IDBLE.EQ.1)           THEN
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
              WRITE(CQLINE(INLINE:INLINE+19),'(A,'')'')') CDWORD(18:)
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
              WRITE(CQLINE(INLINE:INLINE+11),'(A,'')'')') CDWORD(10:)
#endif
              IDBLE  = 0
          ELSEIF (ITYPE.EQ.IFOHOQ)  THEN
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
              WRITE(CQLINE(INLINE:INLINE+19),'(7X,''"'',A12)') IWORD
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
              WRITE(CQLINE(INLINE:INLINE+11),'(7X,''"'',A4)')  IWORD
#endif
          ELSEIF (ITYPE.EQ.IFOINQ)  THEN
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
              WRITE(CQLINE(INLINE:INLINE+19),'(I20)')          IWORD
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
              WRITE(CQLINE(INLINE:INLINE+11),'(I12)')          IWORD
#endif
          ELSEIF (ITYPE.EQ.IFOFLQ)  THEN
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
              WRITE(CQLINE(INLINE:INLINE+19),'(G20.14)')       RWORD
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
              WRITE(CQLINE(INLINE:INLINE+11),'(G12.4)')        RWORD
#endif
          ELSEIF (ITYPE.EQ.IFODOQ)  THEN
              CALL UCOPY(Q(KQS+LS+I),DWORD,2)
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
              WRITE(CDWORD,'(D36.29)') DWORD
              WRITE(CQLINE(INLINE:INLINE+19),'('' D('',A)') CDWORD(:17)
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
              WRITE(CDWORD,'(D20.13)') DWORD
              WRITE(CQLINE(INLINE:INLINE+11),'('' D('',A)') CDWORD(:9)
#endif
              IDBLE  = 1
          ELSE
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
              WRITE(CQLINE(INLINE:INLINE+19),'(''  Z '',Z16)') IWORD
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
              WRITE(CQLINE(INLINE:INLINE+11),'(''  Z '',Z8)')  IWORD
#endif
          ENDIF
          NWPRNT = NWPRNT + 1

#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
          IF (MOD(NWPRNT,6).EQ.0)  THEN
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
          IF (MOD(NWPRNT,10).EQ.0)  THEN
#endif
              IF (IFLOPT(MPOSQQ).NE.0)                    GO TO 999
              CALL DZTEXT(0,CDUMMQ,1)
              WRITE (CQLINE,'(1X,I7,'' /'',120X)') NWPRNT+IFIRST
              INLINE =11
          ELSE
#if defined(CERNLIB_QMCRY)||defined(CERNLIB_QMCDC)
              INLINE = INLINE + 20
#endif
#if (!defined(CERNLIB_QMCRY))&&(!defined(CERNLIB_QMCDC))
              INLINE = INLINE + 12
#endif
          ENDIF

  100 CONTINUE


      JDATA = JDATA + NWSEC

      IF (JDATA.GE.ILAST)                                  GO TO 900

      IF (JDATA.LT.ILAST)   THEN
          JFOCUR = JFOCUR + 2
          IF (JFOCUR.LT.JFOEND)                            GO TO 10
          JFOCUR = JFOREP
                                                           GO TO 10
      ENDIF


  900 IF (INLINE.NE.11.AND.IFLOPT(MPOSQQ).EQ.0) CALL DZTEXT(0,CDUMMQ,1)

  999 CQSTAK(MCQSIQ:) = CHSTAK
      END
