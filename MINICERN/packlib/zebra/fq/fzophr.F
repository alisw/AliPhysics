*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:48  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZOPHR

C-    Write one physical record
C-    service routine to FZOUT

#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/mzcn.inc"
#include "zebra/mzioc.inc"
#include "zebra/fzcx.inc"
C--------------    End CDE                             --------------

#include "fzstamp.inc"

*      Declaratives, DIMENSION etc.
#include "fzophrd1.inc"
* Ignoring t=pass

#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOP, 4HHR   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOPHR /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOPHR  ')
#endif

*      Declaratives, DATA
#include "fzophrd2.inc"


#include "zebra/qtrace.inc"

C------            Book-keeping

      JOP = IDX(2)
      IF (JOP.LT.-1)               GO TO 81
      IQ(KQSP+LQFX+22) = IQ(KQSP+LQFX+22) + 1
      IQ(KQSP+LQFX+33) = IQ(KQSP+LQFX+33) + 1

      JFAST  = IQ(KQSP+LBPARX-6)
      NFASTX = IQ(KQSP+LBPARX-5)
      IF (JFAST.NE.0)              GO TO 17

C--                Steering record

      LQ(L4STAX)   = MCCW1
      LQ(L4STAX+1) = MCCW2
      LQ(L4STAX+2) = MCCW3
      LQ(L4STAX+3) = MCCW4
      LQ(L4STAX+7) = NFASTX
      IQ(KQSP+LQFX+23)= IQ(KQSP+LQFX+23) + 1
#if defined(CERNLIB_QDEVZE)
      J = LQ(L4STAX+5)
      IF (LOGLVX.GE.3)  WRITE (IQLOG,9716) J,NFASTX
 9716 FORMAT (' FZOPHR-  write steering record #',I5,' NFASTX=',2I4)
#endif
      GO TO 19

C--                Fast record

   17 NFASTX = NFASTX - 1
      IQ(KQSP+LBPARX-5) = NFASTX
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVX.GE.3)  WRITE (IQLOG,9717)
 9717 FORMAT (' FZOPHR-  write     fast record')
#endif
   19 IQ(KQSP+LBPARX-6) = NFASTX

C-----------       Output of the record       ----------------------

      LWR  = L4STAX
      NWR  = IQ(KQSP+LBPARX+1)
      NWU  = NWR - NWFILX

#if defined(CERNLIB_FZALFA)
      IF (IFIFOX.EQ.4)             GO TO 78
#endif
#if defined(CERNLIB_FQNEEDCV)

C--                Pack the record, if needed

      IF (IUPAKX.NE.0)             GO TO 24
      NWPK = MAXREX - NWFILX

#endif
#if defined(CERNLIB_FQNEEDCV)
#include "fzophr22.inc"
#endif
#if defined(CERNLIB_FQNEEDCV)
   24 CONTINUE

#endif
#if defined(CERNLIB_FZDACC)
      IF (IFIFOX.EQ.2)             GO TO 51
#endif
#if defined(CERNLIB_FZMEMORY)
      IF (IFIFOX.EQ.3)             GO TO 71
#endif

C-----------------------------------------------------------
C----              Output Sequential
C-----------------------------------------------------------

      IQ(KQSP+LQFX+20) = IQ(KQSP+LQFX+20) + NWR

#if defined(CERNLIB_FZCHANNEL)
C--       Output sequential channel

      IF (IACMOX.EQ.3)  THEN
          CALL JUMPST (IADOPX)
          ICODE     = 1
          IQUEST(1) = LUNX
          IQUEST(2) = NWR
          IQUEST(3) = ISTENX
          IQUEST(4) = 0
          IQUEST(5) = MEDIUX - 4
          CALL JUMPX2 (LQ(LWR),ICODE)
          GO TO 999
        ENDIF

#endif
#if defined(CERNLIB_FZLIBC)
C--       Output sequential with calls to the C library

      IF (IACMOX.EQ.2)  THEN
          CALL CFPUT (IADOPX, MEDIUX, NWR, LQ(LWR), ISW)
          IF (ISW.NE.0)  THEN
              IQUEST(1) = 19
              IQUEST(2) = 21
              IQUEST(3) = ISW
              IQUEST(4) = LUNX
              IQUEST(5) = IADOPX
              CALL ZTELL (19,0)
            ENDIF
          GO TO 999
        ENDIF

#endif
#if defined(CERNLIB_FZFORTRAN)
C--       Output sequential with Fortran calls

#include "fzophr44.inc"
* Ignoring t=pass
#endif
#if defined(CERNLIB_FZFORTRAN)
#include "fzophr45.inc"

#endif
#if defined(CERNLIB_FZDACC)
C-----------------------------------------------------------
C----              Output Direct Access
C-----------------------------------------------------------

   51 JREC = IQ(KQSP+LQFX+33)
      IQ(KQSP+LQFX+20) = IQ(KQSP+LQFX+20) + NWR

#endif
#if defined(CERNLIB_FZDACCH)
C--       Output random channel

      IF (IACMOX.EQ.3)  THEN
          CALL JUMPST (IADOPX)
          ICODE     = 1
          IQUEST(1) = LUNX
          IQUEST(2) = NWR
          IQUEST(3) = ISTENX
          IQUEST(4) = JREC
          IQUEST(5) = MEDIUX - 4
          CALL JUMPX2 (LQ(LWR),ICODE)
          GO TO 999
        ENDIF

#endif
#if defined(CERNLIB_FZDACCL)
C--       Output random with calls to the C library

      IF (IACMOX.EQ.2)  THEN
          CALL CFPUT (IADOPX, MEDIUX, NWR, LQ(LWR), ISW)
          IF (ISW.NE.0)  THEN
              IQUEST(1) = 19
              IQUEST(2) = 22
              IQUEST(3) = ISW
              IQUEST(4) = LUNX
              IQUEST(5) = IADOPX
              CALL ZTELL (19,0)
            ENDIF
          GO TO 999
        ENDIF

#endif
#if defined(CERNLIB_FZDACCF)
#include "fzophr55.inc"

#endif
#if defined(CERNLIB_FZMEMORY)
C-----------------------------------------------------------
C----              Output in memory mode
C-----------------------------------------------------------

   71 LBUF = IQ(KQSP+LQFX+1)
      CALL UCOPY (LQ(LWR),LQ(LBUF),NWU)
      IQ(KQSP+LQFX+1)  = IQ(KQSP+LQFX+1)  + NWU
      IQ(KQSP+LQFX+20) = IQ(KQSP+LQFX+20) + NWU
      GO TO 999

#endif
#if defined(CERNLIB_FZALFA)
C-----------------------------------------------------------
C----              Output ALFA format
C-----------------------------------------------------------

   78 CALL FZOASC (LUNX,L4STAX,L4ENDX,JFAST,L4STOX,0)
      IQ(KQSP+LQFX+20) = IQ(KQSP+LQFX+20) + MAXREX
      GO TO 999

#endif
C-----------------------------------------------------------
C--                ENDFILE
C-----------------------------------------------------------

C-       NEOFM = 1  EoF 1 only       IDX(1) = -2  EoF
C-               2  EOF 2 only                -3  EoD
C-               3  EOF 1 + 2

   81 IF (IFIFOX.GE.2)             GO TO 999
      IF (IACMOX.GE.3)             GO TO 999
      NEOFU = IQUEST(11)
      NEOFM = IQUEST(12)
#if defined(CERNLIB_FZLIBC)
      IF (IACMOX.EQ.2)  THEN
          CALL CFWEOF (IADOPX, MEDIUX, IQUEST(11))
          GO TO 89
        ENDIF
#endif

#include "fzophre1.inc"
* Ignoring t=pass

   87 DO 88 J=1,NEOFU
      ENDFILE LUNX
   88 CONTINUE
   89 CONTINUE
#if defined(CERNLIB_QPRINT)
      IF (LOGLVX.GE.0)   WRITE (IQLOG,9089) LUNX,IQUEST(11)
 9089 FORMAT (' FZOPHR.  LUN=',I4,' Write',I2,' System EOF')
#endif

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
