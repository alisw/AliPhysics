*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:12  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if defined(CERNLIB_FZDACC)
      SUBROUTINE FZUDAT (LUNP,IFROM)

C-    Update the forward reference DAT record

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzcx.inc"
#include "fzhci.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9)

*      Declaratives, DIMENSION etc.
#include "fzudatd1.inc"
* Ignoring t=pass

#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZUD, 4HAT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZUDAT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZUDAT  ')
#endif

*      Declaratives, DATA
#include "fzudatd2.inc"
* Ignoring t=pass


#include "zebra/qtrace.inc"

      LUNNX   = LUNP(1)
      CALL FZLOC (LUNNX,2)

      IF (IFIFOX.NE.2)             GO TO 999
      IF (IACMOX.GE.3)             GO TO 999
#if defined(CERNLIB_QDEBPRI)
      IF (LOGLVX.GE.2)  WRITE (IQLOG,9002) LUNX
 9002 FORMAT (' FZUDAT-  called for LUN=',I4)
#endif

      LBPARX = LQFX + INCBPX
      NWMIN  = IQ(KQSP+LBPARX+1)
      INCBUF = IQ(KQSP+LBPARX+2)
      L4STOX = KQSP+8 + LBPARX + INCBUF

C--       read the first record of the file
      LIN    = L4STOX
#if defined(CERNLIB_FQNEEDPK)
      IF (IUPAKX.EQ.0)  LIN = LIN + IQ(KQSP+LBPARX+3)

#endif
#if defined(CERNLIB_FZLIBC)
      IF (IACMOX.EQ.2)  THEN
          CALL CFTELL (IADOPX, MEDIUX, NWMIN, IPOSCR, ISW)
          IF (ISW.NE.0)            GO TO 999
          CALL CFREW (IADOPX, MEDIUX)
          NWR = NWMIN
          CALL CFGET (IADOPX, MEDIUX, NWMIN, NWR, LQ(LIN), ISW)
          CALL CFSEEK (IADOPX, MEDIUX, NWMIN, IPOSCR, ISTAT)
          GO TO 37
        ENDIF

#endif
#if defined(CERNLIB_FZFORTRAN)
#if defined(CERNLIB_QMIBMD)||defined(CERNLIB_QMIRTD)
      CALL ADREDKD (LUNX, LQ(LIN),NWMIN,1,ISW)
#else
      READ (LUNX,REC=1,IOSTAT=ISW)  (LQ(LIN+J-1),J=1,NWMIN)
#endif
#endif
   37 IF (ISW.NE.0)                GO TO 999

#if defined(CERNLIB_FQNEEDCV)
      IF (IUPAKX.NE.0)             GO TO 41
#endif
#if defined(CERNLIB_FQNEEDCV)
#include "fzudat31.inc"
#endif

C--       check the very first logical record is DAT
   41 CALL FZICHH (1, LQ(L4STOX),-1)
      IF (IQUEST(1).NE.0)          GO TO 999
      IF (NTLRI.NE.8)              GO TO 999
      L = L4STOX + 9
      IF (LQ(L).NE.2)              GO TO 999
      IF (LQ(L+3).NE.2)            GO TO 999

C--       insert the forward reference
      LQ(L+12) = IQ(KQSP+LQFX+34)
      LQ(L+13) = IQ(KQSP+LQFX+35)

C--       over-write the updated record
      LWR = L4STOX
#if defined(CERNLIB_FQNEEDCV)
      IF (IUPAKX.NE.0)             GO TO 47
#endif
#if defined(CERNLIB_FQNEEDCV)
#include "fzudat44.inc"
#endif
#if defined(CERNLIB_FQNEEDCV)
   47 CONTINUE

#endif
#if defined(CERNLIB_FZLIBC)
      IF (IACMOX.EQ.2)  THEN
          CALL CFREW (IADOPX, MEDIUX)
          CALL CFPUT (IADOPX, MEDIUX, NWMIN, LQ(LWR), ISW)
          CALL CFSEEK (IADOPX, MEDIUX, NWMIN, IPOSCR, ISTAT)
          GO TO 87
        ENDIF

#endif
#if defined(CERNLIB_FZFORTRAN)
#if defined(CERNLIB_QMIBMD)||defined(CERNLIB_QMIRTD)
      CALL ADWRTKD (LUNX, LQ(LWR),NWMIN,1,ISW)
      WRITE (LUNX,REC=1,IOSTAT=ISW)  (LQ(LWR+J-1),J=1,NWMIN)
#endif
#endif
   87 IF (ISW.NE.0)                GO TO 999
#if defined(CERNLIB_QDEBPRI)
      IF (IFROM.EQ.0)  THEN
          IF (LOGLVX.GE.2)   WRITE (IQLOG,9087) LUNX
        ELSE
          IF (LOGLVX.GE.-1)  WRITE (IQLOG,9088) LUNX
        ENDIF
 9087 FORMAT (' FZODAT.  LUN=',I4,' DAT forward ref. record updated')
 9088 FORMAT (' FZIDAT.  LUN=',I4,' DAT forward ref. record updated')
#endif

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
