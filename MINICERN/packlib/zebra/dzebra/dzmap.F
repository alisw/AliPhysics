*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:12  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*------------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZMAP
#include "zebra/bankparq.inc"
#include "zebra/divparq.inc"
#include "zebra/storparq.inc"
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/mzcn.inc"
#include "zebra/zbcd.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"


      PARAMETER ( NLMAPQ = 7 )
      PARAMETER ( ISIDEQ = 6 )
      PARAMETER ( ILINKQ = 15)
      PARAMETER ( NLINKQ = 14)
      PARAMETER ( IMAD1Q = 1 , IMAD2Q= 8)
      PARAMETER ( IMTG1Q = 9 , IMTG2Q= 9)
      PARAMETER ( IMID1Q = 10, IMID2Q= 13)

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZMAP' )

#include "zebra/q_jbit.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif


      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      IDOPT  = IFLOPT(MPOSDQ)


      IF (LN.GT.0)      THEN
          CALL MZCHLN(NCHEKQ,LN)
          IF (IQFOUL.NE.0)                                 GO TO 998
      ELSE
          CALL MZCHLS(NCHEKQ,LS)
          IF (IQFOUL.NE.0)                                 GO TO 998
      ENDIF

      IF (IQND.LT.0) THEN
          IF (IFLOPT(MPOSHQ).NE.0) THEN

              WRITE(CQLINE,
#if !defined(CERNLIB_OCTMAP)
     W         '(1X,''(*HO*'',1X,I8,''('',Z8,
#endif
#if defined(CERNLIB_OCTMAP)
     W         '(1X,''(*HO*'',1X,I8,''('',O8,
#endif
     W          '') -- HOLE of '',I8,'' words'')') IQLN,(IQLN+LQSTOR)
*              MAP addresses are in BYTES
#if !defined(CERNLIB_WORDMAP)
     W                                                               *4
#endif
     W          ,-IQND

              CALL DZTEXT(0,CDUMMQ,1)
          ENDIF
          LX = IQNX
                                                           GO TO 999
      ELSE
          LS = IQLS
          NL = IQNL
          NS = IQNS
          ND = IQND
          LX = IQNX
      ENDIF

      IF (IFLOPT(MPOSKQ).EQ.0)         THEN
          JDROP  = JBIT(IQ(KQS+IQLS),IQDROP)
      ELSE
          JDROP  = 0
      ENDIF
      MARKD  = 0
      IF (JDROP.EQ.0)                  THEN
          IF(IFLOPT(MPOSDQ).NE.0)          THEN
              MARKD = JRSBYT(0,IQ(KQS+LS),IQMARK,1)
              IFLOPT(MPOSDQ) = 0
          ENDIF
          IF(IFLOPT(MPOSCQ).NE.0)          THEN
              MARKD = JRSBYT(0,IQ(KQS+LS),IQCRIT,1) + MARKD
          ENDIF
          IF (MARKD+IFLOPT(MPOSFQ).NE.0)       THEN
              IFLOPT(MPOSDQ) = 1
              CALL DZSHPR(LS,0,0,0,0)
                                                           GO TO 999
          ENDIF
          IF (MARKD+IFLOPT(MPOSEQ).NE.0)       THEN
              IFLOPT(MPOSDQ) = 1
              CALL DZSHPR(LS,0,0,0,-1)
                                                           GO TO 999
          ENDIF
      ENDIF

      CALL DZBKHD
      IF (IQUEST(1).NE.0)                                  GO TO 998


      IF (NL.EQ.0)                                         GO TO 999


      CQLINE = '      . LINKS'
      LAST   = LS - NL
      L      = LAST

      DO 43 J=1,NL
          IF (LQ(L+KQS).NE.0)                              GO TO 44
   43 L = L+1
   44 NP= LS - L

      IF (NP.EQ.0)                                         GO TO 999

      IF (NP.GT.NLMAPQ)                 THEN
          CQLINE(ISIDEQ+1:ISIDEQ+1) = '+'
          NP    = NLMAPQ
      ENDIF

      IF (JDROP.NE.0) CQLINE(ISIDEQ:ISIDEQ+1) = '**'

      L = LS + KQS
      DO 50  J=1,NP
          I       = (J-1)*NLINKQ + ILINKQ
          LINK    = LQ (L-J)
          WRITE(CQLINE(I+IMAD1Q:I+IMAD2Q),'(I8)') LINK
          IF (LINK.EQ.LNULL)                               GO TO 50
          CALL MZCHLS(NCHEKQ,LINK)
          WRITE(CQLINE(I+IMID1Q:I+IMID2Q),'(A4)') IQID

          IF (IQFOUL.EQ.0)                 THEN
              IF (JBIT(IQ(KQS+LINK),IQDROP).NE.0)   THEN
                  CQLINE(I+IMTG1Q:I+IMTG2Q) = '('
                  IF (JDROP.EQ.0) CQLINE(ISIDEQ+1:ISIDEQ+1) = 'F'
                  IF (IQND.LT.0)  CQLINE(I+IMID1Q:I+IMID2Q) = '*HO*'
              ELSE
                  IF (JDROP.NE.0.AND.J.LE.NS.AND.J.GT.1)        THEN
                      CQLINE(ISIDEQ+1:ISIDEQ+1) = 'F'
                  ENDIF
              ENDIF
          ELSEIF (IQFOUL.GT.0)             THEN
              IF (J.LE.NS)                     THEN
                  CQLINE(I+IMID1Q:I+IMID2Q) = '****'
                  CQLINE(ISIDEQ+1:ISIDEQ+1) = 'F'
              ELSE
                  CQLINE(I+IMID1Q:I+IMID2Q) = '-'
              ENDIF
          ELSE
              CQLINE(I+IMID1Q:I+IMID2Q) = '****'
              CQLINE(ISIDEQ+1:ISIDEQ+1) = 'F'
          ENDIF

   50 CONTINUE


      CALL DZTEXT(0,CDUMMQ,1)

                                                           GO TO 999

  998 IQUEST(1) = 1

  999 IFLOPT(MPOSDQ) = IDOPT
      CQSTAK(MCQSIQ:) = CHSTAK
      END
