*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:10:49  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE FZOREC

C-    Record controls, exchange file format
C-    service routine to FZOUT

C-    Controlling parameter : IDX(2)
C-
C-    IDX(2)   > 0   start new logical record,
C-                   IDX(1) = length
C-                   IDX(2) = type :  1 start/end of run
C-                                    2 new d/s of new event
C-                                    3 new d/s of same event
C-                                    9 write emergency stop
C-
C-             = 0   up-date buffer parameters from L4CURX
C-                   write PhR if buffer completed :
C-                   either : buffer full and LR continues
C-                       or : LR ended on fast block
C-
C-             = -1  flush the buffer
C-             = -2  End-of-File
C-             = -3  End-of-Data

C-    Note : for each LR there is a final call with IDX(2)=0
C-    to save the buffer parameters; at this moment the buffer
C-    is flushed if the record is fast or almost full

#include "zebra/zmach.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqf.inc"
#include "zebra/fzcx.inc"
C--------------    End CDE                             --------------
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HFZOT, 4HRX   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HFZOREC /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'FZOREC  ')
#endif

#include "zebra/q_or.inc"
#include "zebra/q_jbit.inc"
#include "zebra/q_sbit1.inc"

#include "zebra/qtrace.inc"


C-    buffer parameters :
C-                     -6 JFAST = 0 : current buffer is steering rec
C-                           .NE. 0 : current buffer is fast record
C-                     -5 NFASTX, # of fast records still to be written
C-                     -4
C-                     -3 N4RESX, # of words in LR still to be done
C-                     -2 N4DONX, # of words stored in buffer
C-                     -1
C-              LBPARX +0   maximum size of buffer
C-                     +1   size of PhR, local machine words
C-                     +2   INCBUF : step to buffer
C-                     +3   (off-set for reading)
C-                     +4   off-set for writing
C-                          ( =0 normal, =128 ALFA)
C-                     -1   free to allow packing
C-              L4STOX +0   start of buffer
C-              L4STAX +0   start of buffer accumulation

      LBPARX = LQFX + INCBPX
      JFAST  = IQ(KQSP+LBPARX-6)
      NFASTX = IQ(KQSP+LBPARX-5)
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVX.GE.3)  WRITE (IQLOG,9701) IDX(2),JFAST,NFASTX
 9701 FORMAT (' FZOREC-  entered with IDX(2)=',I3,
     F' JFAST,NFASTX=',2I4)
#endif

      IF   (IDX(2).EQ.0)           GO TO 41

      N4RESX = IQ(KQSP+LBPARX-3)
      N4DONX = IQ(KQSP+LBPARX-2)
      INCBUF = IQ(KQSP+LBPARX+2)
      JOFFSO = IQ(KQSP+LBPARX+4)
      L4STOX = KQSP+8 + LBPARX + INCBUF
      L4STAX = L4STOX + JOFFSO
      L4ENDX = L4STAX + MAXREX
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVX.GE.3)  WRITE (IQLOG,9703) N4RESX,N4DONX
 9703 FORMAT (12X,'N4RESX,N4DONX=',2I7)
#endif
      IF (IDX(2).LT.0)             GO TO 51
      IF (IDX(2).EQ.9)             GO TO 30

C-----------------------------------------------------------
C------            IDX(2) > 0 :  start new logical record
C-----------------------------------------------------------

      IQ(KQSP+LQFX+21) = IQ(KQSP+LQFX+21) + 1

#if defined(CERNLIB_FZMEMORY)
      IF (IFIFOX.EQ.3)             GO TO 23
#endif
      IF (N4RESX.NE.0)             GO TO 31
      IF (N4DONX.NE.0)             GO TO 24

C--     LR starts new physical record

   23 LQ(L4STAX+4) = MAXREX
      IF (MEDIUX.LT.4)  THEN
          LQ(L4STAX+5) = IQ(KQSP+LQFX+23)
        ELSE
          LQ(L4STAX+5) = 0
        ENDIF
      LQ(L4STAX+6) = 8

      IQ(KQSP+LBPARX-6)= 0
      IQ(KQSP+LBPARX-5)= 0
      N4DONX = 8
      GO TO 25

   24 IF (LQ(L4STAX+6).EQ.0)  LQ(L4STAX+6)=N4DONX
   25 IQ(KQSP+LQFX+31) = IQ(KQSP+LQFX+33) + 1
      IQ(KQSP+LQFX+32) = N4DONX
      N4DONX = N4DONX + 2
      L4CURX = L4STAX + N4DONX
      LQ(L4CURX-2) = IDX(1)
      LQ(L4CURX-1) = IDX(2)
      ISTENX = 1

      IF (IDX(2).EQ.1)  THEN
          LQ(L4STAX+4) = MSBIT1(LQ(L4STAX+4),31-JRUNX)
          LQ(L4STAX+5) = 0
        ENDIF

      IQ(KQSP+LBPARX-3) = IDX(1)
      IQ(KQSP+LBPARX-2) = N4DONX

      NWOVFL = N4DONX + IDX(1) - MAXREX
      IF (NWOVFL.LE.0)             GO TO 991

      NFASTX = (NWOVFL-1)/MAXREX + 1
      IF (JBIT(MSTATX,15).EQ.0)   THEN
          NWUNUS = NFASTX*MAXREX - NWOVFL
          IF (NWUNUS .GE. 12)  NFASTX=NFASTX-1
        ENDIF
      IQ(KQSP+LBPARX-5) = NFASTX
      GO TO 991

C----              EMERGENCY STOP

   30 IDX(2) = 0
   31 IF (JFAST.NE.0)              GO TO 36
      IF (NFASTX.NE.0)             GO TO 33
      IF (N4RESX+N4DONX.GT.MAXREX) GO TO 33
      NWTOLR = LQ(L4STAX+6)
      IF (NWTOLR.EQ.0)             GO TO 36
      LQ(L4STAX+NWTOLR) = MAXREX - NWTOLR

   33 NWCL = MAXREX - N4DONX
      IF (NWCL.GT.0)  CALL VZERO (LQ(L4STAX+N4DONX),NWCL)
      CALL FZOPHR

   36 LQ(L4STAX+4) = MAXREX
      IF (MEDIUX.LT.4)  THEN
          LQ(L4STAX+5) = IQ(KQSP+LQFX+23)
        ELSE
          LQ(L4STAX+5) = 0
        ENDIF
      LQ(L4STAX+6) = 8
      LQ(L4STAX+4) = MSBIT1 (LQ(L4STAX+4),32)

      IQ(KQSP+LBPARX-6)= 0
      IQ(KQSP+LBPARX-5)= 0
      IQ(KQSP+LBPARX-3)= 0
      IQ(KQSP+LBPARX-2)= 0

      NWCL = MAXREX - 10
      LQ(L4STAX+8) = NWCL
      LQ(L4STAX+9) = 5
      CALL VZERO (LQ(L4STAX+10),NWCL)
      CALL FZOPHR
      IF (IDX(2).GT.0)             GO TO 23
      GO TO 991

C-----------------------------------------------------------
C------            IDX(2) = 0 :  up-date buffer parameters
C--                              write the buffer if complete
C-----------------------------------------------------------

C-        Note : L4CURX=L4ENDX+1 if a double-precision number spans
C-               two physical records, giving NWFREE = -1

   41 N4DONX = L4CURX - L4STAX
      NWNEW  = N4DONX - IQ(KQSP+LBPARX-2)
      N4RESX = IQ(KQSP+LBPARX-3) - NWNEW
      IQ(KQSP+LBPARX-3) = N4RESX
      IQ(KQSP+LBPARX-2) = N4DONX

   42 NWFREE = L4ENDX - L4CURX
#if defined(CERNLIB_QDEVZE)
      IF (LOGLVX.GE.3)  WRITE (IQLOG,9743) N4RESX,N4DONX,NWFREE
 9743 FORMAT (12X,'N4RESX,N4DONX,NWFREE=',3I7)
#endif
      IF (NWFREE.LT.0)             GO TO 55
      IF (N4RESX.NE.0)  THEN
          IF (NWFREE.EQ.0)         GO TO 55
          GO TO 999
        ENDIF

C--                End of LR reached

      IF (N4DONX.EQ.0)             GO TO 999
      ISTENX = IOR(ISTENX,2)
      IF (NWFREE.LT.4)             GO TO 54
      IF (JFAST.NE.0)              GO TO 54
      IF (JBIT(MSTATX,15).NE.0)    GO TO 54
      GO TO 999

C-----------------------------------------------------------
C------            IDX(2) = -1 :  flush the buffer
C-----------------------------------------------------------

   51 IF (IDX(2).LT.-1)            GO TO 71
      IF (N4RESX.NE.0)             GO TO 31
      IF (N4DONX.EQ.0)             GO TO 991
      ISTENX = IOR(ISTENX,2)
      NWFREE = MAXREX - N4DONX

C--       construct dummy LR to complete the buffer
   54 IF (NWFREE.LE.0)             GO TO 55
      IF (JFAST.EQ.0) THEN
          IF (LQ(L4STAX+6).EQ.0)  LQ(L4STAX+6)=N4DONX
        ENDIF

      L4CURX = L4STAX + N4DONX
      LQ(L4CURX)   = NWFREE - 1
      IF (NWFREE.LT.2)             GO TO 55
      LQ(L4CURX+1) = 5
      NWFILX = NWFREE - 2
      IF (NWFILX.EQ.0)             GO TO 55
#if !defined(CERNLIB_FZMEMORY)
      CALL VZERO (LQ(L4CURX+2),NWFILX)
#endif
#if defined(CERNLIB_FZMEMORY)
      IF (IFIFOX.NE.3)  CALL VZERO (LQ(L4CURX+2),NWFILX)
#endif

C--                Write the physical record

   55 MWOVSV = LQ(L4ENDX)
      CALL FZOPHR

C--                Ready buffer for next physical record

      ISTENX = 0
      IF (N4RESX.GT.0)             GO TO 57
      IF (NWFREE.LT.0)             GO TO 57
C--       no pending data for current logical record
      IQ(KQSP+LBPARX-2) = 0
      L4CURX = L4STAX
      GO TO 991

C--       current logical record continues
   57 IF (IQ(KQSP+LBPARX-6).EQ.0)  THEN
          LQ(L4STAX+4) = MAXREX
          IF (MEDIUX.LT.4)  THEN
              LQ(L4STAX+5) = IQ(KQSP+LQFX+23)
            ELSE
              LQ(L4STAX+5) = 0
            ENDIF
          LQ(L4STAX+6) = 0
          N4DONX = 8
        ELSE
          N4DONX = 0
        ENDIF
      L4CURX = L4STAX + N4DONX
      IQ(KQSP+LBPARX-2) = N4DONX
      IF (NWFREE.GE.0)             GO TO 991

C--       double-precision number spannning the physical records
      LQ(L4CURX) = MWOVSV
      N4DONX = N4DONX + 1
      L4CURX = L4CURX + 1
      IQ(KQSP+LBPARX-2) = N4DONX
      IF (N4RESX.GT.0)             GO TO 991
      GO TO 42

C-----------------------------------------------------------
C------            IDX(2) < -1 :  ENDFILE
C-----------------------------------------------------------

   71 CALL FZOPHR

  991 IDX(2) = 0

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
