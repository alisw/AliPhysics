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
      SUBROUTINE DZIOPD (IFIRST,ILAST)
      SAVE KFOTYP
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzioc.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/bkfoparq.inc"

      PARAMETER (NENLNQ = 5 )
      PARAMETER (NLNGRQ = 10)
      PARAMETER (NLISTQ = 6 )
      PARAMETER (NLINEQ = 24)

      CHARACTER CHROUT*(*),CHSTAK*6,KFOTYP(0:11)*1,CDWORD*31
      DOUBLE PRECISION DWORD
      INTEGER                IWORD
      REAL                         RWORD
      EQUIVALENCE     (DWORD,IWORD,RWORD)
      PARAMETER (CHROUT = 'DZIOPD')
      DATA KFOTYP /'U','B','I','F','D','H','*','S','*','N','*','L'/

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif







      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      INLINE = NLISTQ
      ILINE  = 1
      NLINES = MIN((ILAST-IFIRST)/NENLNQ+1,NLNGRQ)
      DO 5 I=1,NLINES
    5 CQMAP(I) = ' '
      CALL ZPAGE(IQPRNT,NLINES+2)
      CQLINE      = ' DATA part of bank'
      CQLINE(60:) = '--------------------'
      IF(IFLOPT(MPOSQQ).EQ.0) CALL DZTEXT(0,CDUMMQ,1)

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
          WRITE(CQMAP(ILINE)(INLINE+1:INLINE+NLINEQ),
     +     '(I6,'' *'',A1,''('',I13,'')'')')
     +     JDATA,KFOTYP(MIN(ITYPE,6)),NWSEC
          NWPRNT = NWPRNT + 1
      ENDIF
                                                           GO TO 50


   40 JDATA = JDATA + 1
      IWORD = IQ(KQS+LS+JDATA)
      ITYPE = MOD (IWORD,16)
      NWSEC = IWORD/16
      IF (JDATA.GE.IFIRST) THEN
          WRITE(CQMAP(ILINE)(INLINE+1:INLINE+NLINEQ),
     +     '(I6,'' *S('',A1,'','',I11,'')'')')
     +     JDATA,KFOTYP(MIN(ITYPE,6)),NWSEC
          NWPRNT = NWPRNT + 1
      ENDIF

   50 IF (JDATA.GE.IFIRST) THEN
          IF (       (INLINE.EQ.NLISTQ+(NENLNQ-1)*NLINEQ)
     X         .AND. (ILINE.EQ.NLINES)                    )  THEN
              IF (IFLOPT(MPOSQQ).NE.0)                     GO TO 999
              IF (NWPRNT.GT.NLNGRQ*NENLNQ) CALL DZTEXT(1,CDUMMQ,1)
              CALL DZTEXT(0,CDUMMQ,NLINES)
              INLINE = NLISTQ
              ILINE  = 1
              NLINES = MIN((ILAST-JDATA-1)/NENLNQ+1,NLNGRQ)
              DO 55 I=1,NLINES
   55         CQMAP(I) = ' '
          ELSEIF (ILINE.EQ.NLINES) THEN
              ILINE  = 1
              INLINE = INLINE + NLINEQ
          ELSE
              ILINE  = ILINE + 1
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
              WRITE(CQMAP(ILINE)(INLINE+1:INLINE+NLINEQ),
     +                      '(I6,1X,A,'')'')') NWPRNT+IFIRST,CDWORD(16:)
              IDBLE  = 0
          ELSEIF (ITYPE.EQ.IFOHOQ)  THEN
              WRITE(CQMAP(ILINE)(INLINE+1:INLINE+NLINEQ),
     +                      '(I6,5X,''"'',A12)') NWPRNT+IFIRST,IWORD
          ELSEIF (ITYPE.EQ.IFOINQ)  THEN
              WRITE(CQMAP(ILINE)(INLINE+1:INLINE+NLINEQ),
     +                               '(I6,I18)') NWPRNT+IFIRST,IWORD
          ELSEIF (ITYPE.EQ.IFOFLQ)  THEN
              WRITE(CQMAP(ILINE)(INLINE+1:INLINE+NLINEQ),
     +                             '(I6,G18.8)') NWPRNT+IFIRST,RWORD
          ELSEIF (ITYPE.EQ.IFODOQ)  THEN
              CALL UCOPY(Q(KQS+LS+I),DWORD,2)
              WRITE(CDWORD,'(D31.24)') DWORD
              WRITE(CQMAP(ILINE)(INLINE+1:INLINE+NLINEQ),
     +                '(I6,1X,''D('',A)') NWPRNT+IFIRST,CDWORD(:15)
              IDBLE  = 1
          ELSE
              WRITE(CQMAP(ILINE)(INLINE+1:INLINE+NLINEQ),
     +                        '(I6,'' Z'',Z16)') NWPRNT+IFIRST,IWORD
          ENDIF
          NWPRNT = NWPRNT + 1

          IF (       (INLINE.EQ.NLISTQ+(NENLNQ-1)*NLINEQ)
     X         .AND. (ILINE.EQ.NLINES)                    )  THEN
              IF (IFLOPT(MPOSQQ).NE.0)                     GO TO 999
              IF (NWPRNT.GT.NLNGRQ*NENLNQ) CALL DZTEXT(1,CDUMMQ,1)
              CALL DZTEXT(0,CDUMMQ,NLINES)
              INLINE = NLISTQ
              ILINE  = 1
              NLINES = MIN((ILAST-NWPRNT-1)/NENLNQ+1,NLNGRQ)
              DO 65 J=1,NLINES
   65         CQMAP(J) = ' '
          ELSEIF (ILINE.EQ.NLINES) THEN
              ILINE  = 1
              INLINE = INLINE + NLINEQ
          ELSE
              ILINE  = ILINE + 1
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


  900 IF (IFLOPT(MPOSQQ).NE.0)                             GO TO 999
          IF (INLINE.NE.NLISTQ.OR.ILINE.NE.1) THEN
              IF (NWPRNT.GT.NLNGRQ*NENLNQ) CALL DZTEXT(1,CDUMMQ,1)
              CALL DZTEXT(0,CDUMMQ,NLINES)
          ENDIF

  999 CQSTAK(MCQSIQ:) = CHSTAK
      END
