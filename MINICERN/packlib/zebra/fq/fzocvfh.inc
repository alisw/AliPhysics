*
* $Id$
*
* $Log$
* Revision 1.2  1997/03/14 17:19:52  mclareni
* WNT mods
*
* Revision 1.1.1.1.2.1  1997/01/21 11:33:42  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
*   convert Display Code -> ASCII
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMCDC)
C--       CDC Display code -> ASCII
      NTHR = NWDO
  522 NTH  = MIN (NTHR,24)
      NCH  = 4*NTH
      CALL UPKCH (MS(JMS+1),IQUEST,NCH,IPKVU)
      DO 529  JCH=1,NCH
      J = IQUEST(JCH)
  529 IQUEST(JCH) = MTAB(J+1)
      CALL PKCHAR (IQUEST,MT(JMT+1),NCH,IPKVP)
      JMT  = JMT + NTH
      JMS  = JMS + NTH
      NTHR = NTHR - NTH
      IF (NTHR.NE.0)               GO TO 522
      GO TO 801

*    --------------  sequences for input   ---------------------
#endif
*   right adjust 4 ASCII characters
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMCRY)
C--       Cray ASCII : keep only the first 4 characters
      DO 519  JL=1,NWDO
      MT(JMT+1) = SHIFTR (MS(JMS+1), 32)
      JMT = JMT + 1
  519 JMS = JMS + 1
      GO TO 801
#endif
*   right adjust 4 ASCII characters
*
* fzocvfh.inc
*
#if (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
C--       Convex ASCII : keep only the first 4 characters
      DO 519  JL=1,NWDO
      MT(JMT+1) = ISHFT (MS(JMS+1), -32)
      JMT = JMT + 1
  519 JMS = JMS + 1
      GO TO 801
#endif
*   invert ASCII characters
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
#include "fzcvhpc.inc"
*    --------------  sequences for input   ---------------------
#endif
*    --------------  conversions for input   -------------------
*   convert EBCDIC -> ASCII
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMIBM)
#include "fzcvhibm.inc"
*    --------------  conversions for input   -------------------
#endif
*   convert EBCDIC -> ASCII
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMIBMD)
#include "fzcvhibm.inc"
*    --------------  Added by MRO             ------------------
#endif
*   invert ASCII characters
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMLNX)
#include "fzcvhlnx.inc"
*    --------------  sequences for input   ---------------------
#endif
*   invert ASCII characters
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMTMO)
#include "fzcvhtmo.inc"
*    --------------  sequences for input   ---------------------
#endif
*    --------------  conversions for input   -------------------
*   invert ASCII characters
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMVAX)
#include "fzcvhvax.inc"
*    --------------  sequences for input   ---------------------
#endif
*   invert ASCII characters
*
* fzocvfh.inc
*
#if defined(CERNLIB_QMVMI)
#include "fzcvhvmi.inc"
*    --------------  sequences for input   ---------------------
#endif
