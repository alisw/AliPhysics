*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:29  mclareni
* Kernlib
*
*
      SUBROUTINE GETWDF (TEXT)

C CERN PROGLIB# Z265    GETWDF          .VERSION KERNVAX  2.36  910306
C ORIG. 31/07/90, Federico Carminati
C MOD.  11/01/91, JZ, change name to GETWDF

C-    Simulate UNIX system call

      COMMON /SLATE/ NTEXT, ISTAT, ISLATE(38)
      CHARACTER    TEXT*(*)
      INTEGER      SYS$SETDDIR

      ISTAT = LIB$SYS_TRNLOG ('SYS$DISK',NDISK,TEXT,,,%VAL(0))

      IF (.NOT.ISTAT)              GO TO 99

      ISTAT = SYS$SETDDIR(%VAL(0), NDIR, TEXT(NDISK+1:))

      IF (.NOT.ISTAT)              GO TO 99

      NALL = NDISK + NDIR
      TEXT(NALL+1:) = ' '
      NTEXT = NALL
      ISTAT = 0
      RETURN

   99 TEXT = ' '
      NTEXT = 0
      RETURN
      END
