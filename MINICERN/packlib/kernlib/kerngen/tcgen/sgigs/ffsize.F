*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      SUBROUTINE FFSIZE (LUN,NWREC,JREC,ISTAT)

C-    Get the size of a disk file connected to Fortran unit LUN

      INTEGER      FSEEK, FTELL

      ISTAT = FSEEK (LUN, 0, 2)
      IF (ISTAT.NE.0)              RETURN

      NBYTES = FTELL (LUN)
      IF (NBYTES.LT.0)  THEN
          ISTAT = -NBYTES
          RETURN
        ENDIF

      NWORDS = NBYTES / 4
      JREC   = NWORDS / NWREC
      RETURN
      END
