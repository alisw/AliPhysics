*
* $Id$
*
* $Log$
* Revision 1.2  1997/03/14 17:19:51  mclareni
* WNT mods
*
* Revision 1.1.1.1.2.1  1997/01/21 11:33:41  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*   convert ASCII -> Display Code
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMCDC)
      NTHR = NWDO
  522 NTH  = MIN (NTHR,24)
      NCH  = 4*NTH
      CALL UPKCH (MS(JMS+1),IQUEST,NCH,IPKVU)
      DO 529  JCH=1,NCH
      J = MIN (MAX(IQUEST(JCH)-JCHMIN,0), JCHMAX)
  529 IQUEST(JCH) = MTAB(J+1)
      CALL PKCHAR (IQUEST,MT(JMT+1),NCH,IPKVP)
      JMT  = JMT + NTH
      JMS  = JMS + NTH
      NTHR = NTHR - NTH
      IF (NTHR.NE.0)               GO TO 522
      GO TO 801
#endif
*   left adjust 4 ASCII ch., blank fill
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMCRY)
      DO 519  JL=1,NWDO
      MT(JMT+1) = SHIFTL (MS(JMS+1), 32) .OR. IBL32
      JMT = JMT + 1
  519 JMS = JMS + 1
      GO TO 801
#endif
*   left adjust 4 ASCII ch., blank fill
*
* fzicvfh.inc
*
#if (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      DO 519  JL=1,NWDO
      MT(JMT+1) = ISHFT (MS(JMS+1), 32) .OR. IBL32
      JMT = JMT + 1
  519 JMS = JMS + 1
      GO TO 801
#endif
*   invert ASCII characters
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
#include "fzcvhpc.inc"
#endif
*   convert ASCII -> EBCDIC
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMIBM)
#include "fzcvhibm.inc"
#endif
*   convert ASCII -> EBCDIC
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMIBMD)
#include "fzcvhibm.inc"
#endif
*   invert ASCII characters
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMLNX)
#include "fzcvhlnx.inc"
#endif
*   invert ASCII characters
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMTMO)
#include "fzcvhtmo.inc"
#endif
*   invert ASCII characters
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMVAX)
#include "fzcvhvax.inc"
#endif
*   invert ASCII characters
*
* fzicvfh.inc
*
#if defined(CERNLIB_QMVMI)
#include "fzcvhvmi.inc"
#endif
