*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:12:06  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:17  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZSTOR (IXSTOR,CHNAME,CHOPT
     +,                  IFENCE,LV,LLR,LLD,LIMIT,LAST)

C-    Initialize new Zebra store region, user called

#include "zebra/zbcd.inc"
#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcwk.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXSTOR(9),IFENCE(9)
      DIMENSION    LV(9),LLR(9),LLD(9),LIMIT(9),LAST(9)
      DIMENSION    MMSYSL(5), NAMELA(2), NAMESY(2)
      CHARACTER    *(*) CHNAME,CHOPT
#if defined(CERNLIB_QMVDS)
      SAVE         MMSYSL, NAMELA, NAMESY, NAMWSP, NAMEDV
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZST, 4HOR   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZSTOR /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZSTOR  ')
#endif
#if defined(CERNLIB_QHOLL)
      DATA  MMSYSL / 4HSYSL,0,0,101,2/
      DATA  NAMELA / 4Hsyst, 4Hem   /
      DATA  NAMESY / 4Hsyst, 4Hem   /
      DATA  NAMWSP / 4Hqwsp /
      DATA  NAMEDV / 4HQDIV /
#endif
#if !defined(CERNLIB_QHOLL)
      DATA  MMSYSL / 0,0,0,101,2/
#endif

#include "zebra/q_sbit1.inc"
#include "zebra/q_shiftl.inc"
#include "zebra/q_locf.inc"


C--                Clear Zebra tables on first entry

      IF (NQSTOR.NE.-1)            GO TO 13
      CALL VZERO (NQOFFT,32)
      LQATAB = LOCF (IQTABV(1)) - 1
      LQASTO = LOCF (LQ(1)) - 1
      LQBTIS = LQATAB - LQASTO
      LQWKTB = LOCF(IQWKTB(1)) - LQASTO
      LQWKFZ = LOCF(IQWKFZ(1)) - LQASTO
      NQTSYS = LOCF(IQDN2(20)) - LQATAB
      NQWKTB = NQWKTT

C-      KQFT=342 relies on LQFSTA(1) to be LQSTA(1+342) in /MZCC/
      KQFT = 342
#if defined(CERNLIB_QPRINT)
      IF (NQLOGD.GE.-1)
     +WRITE (IQLOG,9011) LQATAB,LQATAB
 9011 FORMAT (1X/' MZSTOR.  ZEBRA table base TAB(0) in /MZCC/ at adr'
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))
     F,I12,1X,O11,' OCT')
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))&&(!defined(CERNLIB_B64))
     F,I12,1X,Z11,' HEX')
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))&&(defined(CERNLIB_B64))
     F,I12,1X,Z16,' HEX')
#endif
   13 CONTINUE
#if !defined(CERNLIB_QHOLL)
      CALL UCTOH ('SYSL',    MMSYSL, 4,4)
      CALL UCTOH ('system  ',NAMELA, 4,8)
      CALL UCTOH ('system  ',NAMESY, 4,8)
      CALL UCTOH ('qwsp'    ,NAMWSP, 4,4)
      CALL UCTOH ('QDIV'    ,NAMEDV, 4,4)
#endif

#include "zebra/qtrace.inc"

      CALL UOPTC (CHOPT,'Q:',IQUEST)
      LOGQ   = IQUEST(1)
      IFLSPL = IQUEST(2)

      JQSTOR = NQSTOR + 1
      CALL VZERO (KQT,27)

C--                Calculate store off-set

      LQSTOR = LOCF(LV(1)) - 1
      KQS    = LQSTOR - LQASTO

      NFEND  = (LQSTOR+1) - LOCF(IFENCE(1))
      NQFEND = NFEND

C--                Printing name of store

      NQSNAM(1) = IQBLAN
      NQSNAM(2) = IQBLAN
      N = MIN (8, LEN(CHNAME))
      IF (N.NE.0)  CALL UCTOH (CHNAME,NQSNAM,4,N)

C--                Set log level

      NQLOGL = NQLOGD
      IF (LOGQ.NE.0)  NQLOGL=-2

C--                Permanent links et al.

      NQSTRU = LOCF(LLR(1)) - (LQSTOR+1)
      NQREF  = LOCF(LLD(1)) - (LQSTOR+1)
      NQLINK = NQREF
      LQ2END = LOCF(LIMIT(1)) - LQSTOR
      NDATAT = LOCF(LAST(1))  - LQSTOR

C--                Calculate table off-set

      NDATA = NDATAT
      LOCT  = LQATAB
      IF (JQSTOR.NE.0)  THEN
          NDATA = NDATA  - NQTSYS
          NQSNAM(6) = NDATA
          LOCT  = LQSTOR + NDATA
          KQT   = LOCT   - LQATAB
          NDATA = NDATA - 4
          CALL VFILL (LQ(KQS+NDATA),10,IQNIL)
        ENDIF

#if defined(CERNLIB_QPRINT)
      IF (NQLOGL.GE.-1)
     +WRITE (IQLOG,9021) JQSTOR,NQSNAM(1),NQSNAM(2)
     +,                  LQSTOR,LOCT,LQSTOR,LOCT,KQS,KQT,KQS,KQT
     +,                  NQSTRU,NQREF,LQ2END,NDATAT,NFEND
 9021 FORMAT (1X/' MZSTOR.  Initialize Store',I3,'  in ',2A4,
     F/10X,'with Store/Table at absolute adrs',2I12
#endif
#if (defined(CERNLIB_QPRINT))&&(!defined(CERNLIB_HEX))
     F/40X,'OCT',2(1X,O11)/40X,'OCT',2(1X,O11)
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))&&(!defined(CERNLIB_B64))
     F/40X,'HEX',2(1X,Z11)/40X,'HEX',2(1X,Z11)
#endif
#if (defined(CERNLIB_QPRINT))&&(defined(CERNLIB_HEX))&&(defined(CERNLIB_B64))
     F/30X,'HEX',2(1X,Z16)/30X,'HEX',2(1X,Z16)
#endif
#if defined(CERNLIB_QPRINT)
     F/30X,'relative adrs',2I12
     F/10X,'with',I6,' Str. in',I6,' Links in',I7,' Low words in'
     F,I8,' words.'
     F/10X,'This store has a fence of',I5,' words.')
#endif

C--                Set minimum sizes

      NSYS   =  400
      NQMINR =   40
      NWF    = 2000
      IF (JQSTOR.EQ.0)  NQMINR=164

C--                Check parameters valid

      IF (NQSTRU.LT.0)               GO TO 91
      IF (NQREF .LT.NQSTRU)          GO TO 91
      IF (NDATAT.LT.NQLINK+NWF)      GO TO 91
      IF (LQ2END.LT.NQLINK+NQMINR)   GO TO 91
      IF (NFEND .LT.1)               GO TO 92
      IF (NFEND .GE.1001)            GO TO 92
      IF (IFLSPL.EQ.1)  THEN
          IF (JQSTOR.EQ.0)           GO TO 96
          GO TO 39
        ENDIF

#if (!defined(CERNLIB_QSINGLST))&&(defined(CERNLIB_QDEBUG))
C--                Check overlapping stores

      IF (JQSTOR.EQ.0)             GO TO 41
      KSA = KQS - NQFEND
      KSE = KQS + NDATAT

      DO 36  JSTO=1,JQSTOR
      JT  = NQOFFT(JSTO)
      JS  = NQOFFS(JSTO)
      JSA = JS  - IQTABV(JT+2)
      JSE = JS  + LQSTA(JT+21)
      JTA = JT  + LQBTIS
      JTE = JTA + NQTSYS

      IF (KSE.GT.JTA .AND. KSA.LT.JTE)    GO TO 94
      IF (KSE.GT.JSA .AND. KSA.LT.JSE)    GO TO 95
   36 CONTINUE
#endif
   39 IF (JQSTOR.GE.16)            GO TO 93

C----              Initialize divisions 1 + 2 + system

   41 NQOFFT(JQSTOR+1) = KQT
      NQOFFS(JQSTOR+1) = KQS
      NQALLO(JQSTOR+1) = IFLSPL
      CALL VZERO (IQTABV(KQT+1),NQTSYS)
      CALL VBLANK (IQDN1(KQT+1), 40)
      NQSTOR = NQSTOR + 1

      LQ(KQS+NDATA-1) = IQNIL
      LQ(KQS+NDATA)   = IQNIL

      NDATA = NDATA - 2
      LQSTA(KQT+21) = NDATA

      JQDVLL = 2
      JQDVSY = 20
      LQSTA(KQT+20)  = NDATA
      LQEND(KQT+20)  = NDATA
      NQDMAX(KQT+20) = NDATA
      IQMODE(KQT+20) = 1
      IQKIND(KQT+20) = ISHFTL (1, 23)
      IQRNO(KQT+20)  = 9437183
      IQDN1(KQT+20)  = NAMESY(1)
      IQDN2(KQT+20)  = NAMESY(2)

      LQSTA(KQT+2)  = NDATA - NSYS
      LQEND(KQT+2)  = LQSTA(KQT+2)
      NQDMAX(KQT+2) = NDATA
      IQMODE(KQT+2) = 1
      IQKIND(KQT+2) = MSBIT1 (2, 21)
      IQRCU(KQT+2)  = 3
      IQRTO(KQT+2)  = ISHFTL (3,20)
      IQRNO(KQT+2)  = 9437183
      IQDN1(KQT+2)  = NAMEDV
      IQDN2(KQT+2)  = IQNUM(3)

      LQSTA(KQT+1)  = NQLINK + 1
      LQEND(KQT+1)  = LQSTA(KQT+1)
      NQDMAX(KQT+1) = NDATA
      IQKIND(KQT+1) = MSBIT1 (1, 21)
      IQRCU(KQT+1)  = 3
      IQRTO(KQT+1)  = ISHFTL (3,20)
      IQRNO(KQT+1)  = 9437183
      IQDN1(KQT+1)  = NAMEDV
      IQDN2(KQT+1)  = IQNUM(2)

      CALL UCOPY (IQCUR,IQTABV(KQT+1),16)
      CALL VFILL (IFENCE,NFEND,IQNIL)
      IF (NQLINK.NE.0)  CALL VZERO (LV,NQLINK)

C--                Return IXSTOR

      IF (JQSTOR.EQ.0)  THEN
          IF (IXSTOR(1).EQ.0)      GO TO 71
        ENDIF
      IDN = ISHFTL (JQSTOR,26)
      IXSTOR(1) = IDN

C----              Create system link table bank

   71 JQDIVI = JQDVSY
      CALL MZLIFT (-7,LSYS,0,2,MMSYSL,0)
      LQSYSS(KQT+1) = LSYS

      NALL   = LOCF(IQTDUM(1)) - LOCF(LQSYSS(1))
      NSTR   = LOCF(LQSYSR(1)) - LOCF(LQSYSS(1))

      LOCAR  = LOCF (LQSYSS(KQT+1)) - LQSTOR
      LOCARE = LOCAR + NALL

C--                Working space

      IQ(KQS+LSYS+1) = 11
      IQ(KQS+LSYS+2) = 1
      IQ(KQS+LSYS+3) = 1 + NQLINK
      IQ(KQS+LSYS+4) = NQSTRU
      IQ(KQS+LSYS+5) = NAMWSP
      IQ(KQS+LSYS+6) = IQBLAN

C--                System link area

      IQ(KQS+LSYS+7) = LOCAR
      IQ(KQS+LSYS+8) = LOCARE
      IQ(KQS+LSYS+9) = NSTR
      IQ(KQS+LSYS+10)= NAMELA(1)
      IQ(KQS+LSYS+11)= NAMELA(2)

C--                Range of possible values for an origin-link

      IQTABV(KQT+13) = MIN (1, LOCAR)
      IQTABV(KQT+14) = MAX (LQSTA(KQT+21), LOCARE)

#include "zebra/qtrace99.inc"
      RETURN

C------            Error conditions

#if (!defined(CERNLIB_QSINGLST))&&(defined(CERNLIB_QDEBUG))
   95 NQCASE = 1
   94 NQCASE = NQCASE - 2
      NQFATA = 3
      IQUEST(20) = JSTO - 1
      IQUEST(21) = NQPNAM(JT+1)
      IQUEST(22) = NQPNAM(JT+2)
#endif

   96 NQCASE = NQCASE + 3
   93 NQCASE = NQCASE + 1
   92 NQCASE = NQCASE + 1
   91 NQCASE = NQCASE + 1
      NQFATA = NQFATA + 9
      IQUEST(11) = NQSNAM(1)
      IQUEST(12) = NQSNAM(2)
      IQUEST(13) = NFEND
      IQUEST(14) = NQSTRU
      IQUEST(15) = NQLINK
      IQUEST(16) = LQ2END
      IQUEST(17) = NDATAT
      IQUEST(18) = NQMINR
      IQUEST(19) = NWF
#include "zebra/qtofatal.inc"
      END
*      ==================================================
#include "zebra/qcardl.inc"
