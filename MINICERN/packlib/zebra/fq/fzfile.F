*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:23  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZFILE (LUNP,LRECP,CHOPT)

C-    Initialize Sequential Zebra I/O unit, User called

#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/mzcwk.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzci.inc"
#include "zebra/fzcx.inc"
#include "zebra/fzcf.inc"
#include "fzficc.inc"
C--------------    End CDE                             --------------
      DIMENSION    LUNP(9),LRECP(9)
      DIMENSION    MMFZ(5)
      CHARACTER    CHOPT*(*)
      CHARACTER    HOLD*12
      CHARACTER    VIDQQ*(*)
#if defined(CERNLIB_PY_VS5)
*      PARAMETER (VIDQQ= '@(#)? >')
#include "zebra/qfhead.inc"
#endif
#if !defined(CERNLIB_PY_VS5)
      PARAMETER   (VIDQQ = '@(#)' //
     +  'ZEFQ     3.77 '
     + // '>')
#endif
#if defined(CERNLIB_QMVDS)
      SAVE         MMFZ
#endif
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZFI, 4HLE   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZFILE /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZFILE  ')
#endif
#if defined(CERNLIB_QHOLL)
      DATA  MMFZ   / 4HFZ  , 2, 1, -1, 2 /
#endif
#if !defined(CERNLIB_QHOLL)
      DATA  MMFZ   / 0,      2, 1, -1, 2 /
#endif

#include "zebra/q_or.inc"
#include "zebra/q_shiftl.inc"

#include "zebra/qtrace.inc"

      LUN   = LUNP(1)
      LREC  = LRECP(1)
      LUNPTR= IQUEST(1)

      CALL FZFICR (0, CHOPT)

      LOGLV = NQLOGD
      IF (IOPTQ.NE.0)        LOGLV = MIN (LOGLV, -2)
#if defined(CERNLIB_QPRINT)
      IF (LOGLV.GE.0)  WRITE (IQLOG,9001) LUN,CHOPT
 9001 FORMAT (1X/' FZFILE.  LUN=',I3,' initialize for OPT= ',A)
#endif

      IF (NQSTOR.LT.0)
     +        CALL ZFATAM ('FZFILE - no dynamic store initialized.')

      IF (LUN.LE.0)                GO TO 92
      NEOFOP = -1
      IDAFO  = 0
      IACMO  = 0
      IADOP  = 0
      INCBP  = 0
      INCBUF = 0
      NWLIFT = 0
      NWXBUF = 0
      JOFFSI = 0
      JOFFSO = 0

#if defined(CERNLIB_FQXISN)
      IUPAK  = 1
#endif
#if !defined(CERNLIB_FQXISN)
      IUPAK = MIN (1, IOPTN+IOPTU+IOPTA)
#endif

      LRECX  = 30 * (LREC/30)
      IF (LRECX.LE.0)  LRECX = 900
      NWRECX = LRECX
#if !defined(CERNLIB_B32)||defined(CERNLIB_QMIRTD)
#include "fzfilelg.inc"
#endif

C--------------------------------------------------------
C----              Memory mode
C--------------------------------------------------------

      IF (IOPTM.EQ.0)              GO TO 24
#if (defined(CERNLIB_FZMEMORY))&&(!defined(CERNLIB_FQXISN))
      IF (IOPTN.EQ.0)  IDAFO = 1
#endif
#if defined(CERNLIB_FZMEMORY)
      MEDIUM = 6
      IFIFO  = 3
      IOPTR  = 0
      JOFFSI = 0
      GO TO 59
#endif
#if !defined(CERNLIB_FZMEMORY)
      GO TO 93
#endif

C--------------------------------------------------------
C----              Alfa mode
C--------------------------------------------------------

   24 MEDIUM = IOPTT
      IF (IOPTA.EQ.0)              GO TO 31
#if (defined(CERNLIB_FZALFA))&&(!defined(CERNLIB_FQXISN))
      IDAFO  = 1
#endif
#if defined(CERNLIB_FZALFA)
      IFIFO  = 4
      LREC   = 900
      NWREC  = 0
      JOFFSI = 0

      IF (IOPTO.NE.0)  THEN
          NWXBUF = 128
          JOFFSO = NWXBUF
        ENDIF
      GO TO 67
#endif
#if !defined(CERNLIB_FZALFA)
      GO TO 93
#endif

C-------------------------------------------------
C--------          Mode not Memory, not Alfa, check L,K,C
C-------------------------------------------------

   31 IFIFO = IABS (IOPTX)
#if !defined(CERNLIB_FQXISN)
      IF (IOPTN.EQ.0)  IDAFO = IFIFO
#endif

      IF (IFIFO.EQ.0)              GO TO 41

C--                Access channel

      IF (IOPTC.EQ.0)              GO TO 34
#if defined(CERNLIB_FZCHANNEL)
      MEDIUM = 4 + IOPTT
      IACMO  = 3
      IOPTR  = 0
      GO TO 36
#endif
#if !defined(CERNLIB_FZCHANNEL)
      GO TO  93
#endif

C--                Access L or K

   34 IF (IOPTL+IOPTK.EQ.0)        GO TO 36
#if defined(CERNLIB_FZLIBC)
      IACMO  = 2
      IADOP  = LUNPTR
      MEDIUM = IOPTT + 2*IOPTK
      IF (IADOP.LT.0)              GO TO 94
#endif
#if !defined(CERNLIB_FZLIBC)
      GO TO 93
#endif

C-------------------------------------------------
C----              Mode Direct-access
C-------------------------------------------------

   36 IF (IOPTD.EQ.0)              GO TO 51
#if (defined(CERNLIB_FZDACC))&&(!defined(CERNLIB_FZDACCH))
      IF (IOPTC.NE.0)              GO TO 93
#endif
#if (defined(CERNLIB_FZDACC))&&(!defined(CERNLIB_FZDACCL))
      IF (IOPTL+IOPTK.NE.0)        GO TO 93
#endif
#if defined(CERNLIB_FZDACC)
      IFIFO  = 2
      GO TO 59
#endif
#if !defined(CERNLIB_FZDACC)
      GO TO 93
#endif

C-------------------------------------------------
C----              Mode Sequential
C-------------------------------------------------

C------            Native File Format

   41 NEOFOP = IOPTVF(1)
#if defined(CERNLIB_FZFFNAT)
      IUPAK = 1
      IF (LREC.LE.0)  LREC=2440

#include "fzfilen1.inc"
* Ignoring t=pass

      IF (IOPTO.NE.0)  LREC = MIN (LREC, 2499)
      LRECPR = LREC
      LREC   = LREC - 2
      NWLIFT = 170

#include "fzfilen2.inc"
* Ignoring t=pass

      GO TO 71
#endif
#if !defined(CERNLIB_FZFFNAT)
      GO TO 93
#endif

C----             Sequential Exchange File Format

   51 IF (IOPTC.NE.0)              GO TO 59
      NEOFOP = IOPTVF(1)
      IF (IOPTL+IOPTK.NE.0)        GO TO 59

#include "fzfilex7.inc"
* Ignoring t=pass

   59 LREC  = LRECX
      NWREC = NWRECX

   67 LRECPR = LREC
      INCBP  = 61
      IF (IOPTI.NE.0)  INCBP = 141
      INCBUF = 7

C-------------------------------------------------
C----              Ready new unit
C-------------------------------------------------

   71 LUNI = 0
      LUNX = 0

#include "fzfiblki.inc"

      CALL FZLOC (LUN,0)
      IF (LUNF.NE.0)               GO TO 91

C--       Construct text string of options used
C--       also :  IOPTVF(J) = IABS(IOPTVF(J))

      IOPTVF(1) = NEOFOP
      CALL FZFICR (7, HOLD)

C-------------------------------------------------
C-                 Lift and fill control bank
C-------------------------------------------------

#if !defined(CERNLIB_QHOLL)
      CALL UCTOH ('FZ  ', MMFZ, 4,4)
#endif
      IF (NWLIFT.EQ.0)  NWLIFT = INCBP + INCBUF + LREC + NWXBUF + 4
      MMFZ(4) = NWLIFT
      CALL MZLIFT (JQPDVS,LQFF,LQFS,1,MMFZ,0)

      NEOF = MAX (0, NEOFOP)

      MSTAT = IOR (MEDIUM,ISHFTL(IFIFO, 3))
      MSTAT = IOR (MSTAT, ISHFTL(IDAFO, 6))
      MSTAT = IOR (MSTAT, ISHFTL(IACMO, 7))
      MSTAT = IOR (MSTAT, ISHFTL(IOPTI, 10))
      MSTAT = IOR (MSTAT, ISHFTL(IOPTO, 11))
      MSTAT = IOR (MSTAT, ISHFTL(NEOF,  12))
      MSTAT = IOR (MSTAT, ISHFTL(IOPTS, 14))
      MSTAT = IOR (MSTAT, ISHFTL(IUPAK, 15))

      IQ(KQS+LQFF-5) = LUN
      IQ(KQS+LQFF)   = IOR (IQ(KQS+LQFF), MSTAT)
      IQ(KQS+LQFF+1) = IADOP
      IQ(KQS+LQFF+3) = INCBP
      IQ(KQS+LQFF+4) = LOGLV
      IQ(KQS+LQFF+5) = LREC

      IF (IFIFO.NE.0)  THEN
          LBPAR = LQFF + INCBP
          IQ(KQS+LBPAR)   = LREC
          IQ(KQS+LBPAR+1) = NWREC
          IQ(KQS+LBPAR+2) = INCBUF
          IQ(KQS+LBPAR+3) = JOFFSI
          IQ(KQS+LBPAR+4) = JOFFSO
        ENDIF

C--                Print options used

#if defined(CERNLIB_QPRINT)
      NHOLD = LNBLNK (HOLD)
      IF (IFLERR.EQ.0)  THEN
          IF (LOGLV.LT.-1)         GO TO 88
          WRITE (IQLOG,9086) LRECPR,HOLD(1:NHOLD)
        ELSE
          WRITE (IQLOG,9087) LRECPR,HOLD(1:NHOLD)
        ENDIF

 9086 FORMAT (' FZFILE.  Use LREC=',I5,', options= ',A)
 9087 FORMAT (' FZFILE.  !!! Incompatible options, use LREC=',I5,
     F ', options= ',A)
#endif

C----              REWIND

   88 CALL FZLOC (LUN,0)
      IF (IOPTR.NE.0)        CALL FZMACH (0)
      IQUEST(1) = 0
      GO TO 999

C-------------------------------------------------
C-                 Error handling
C-------------------------------------------------

   91 IQUEST(1) = 1
      IF (IOPTP.NE.0)              GO TO 999
      CALL ZFATAM ('FZFILE - File already open.')

   92 IQUEST(1) = 2
      IF (IOPTP.NE.0)              GO TO 999
      CALL ZFATAM ('FZFILE - LUN invalid.')

   93 IQUEST(1) = 3
      IF (IOPTP.NE.0)              GO TO 999
      CALL ZFATAM ('FZFILE - Requested format not available.')
      CALL ZFATAM (VIDQQ)

#if defined(CERNLIB_FZLIBC)
   94 IQUEST(1) = 4
      IF (IOPTP.NE.0)              GO TO 999
      CALL ZFATAM ('FZFILE - L mode and IQUEST(1) < 0 on call.')
#endif

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
