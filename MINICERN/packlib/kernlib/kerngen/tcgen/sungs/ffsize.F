*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/02 23:16:08  thakulin
* Add support for Apogee Fortran.
*
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      SUBROUTINE FFSIZE (LUN,NWREC,JREC,ISTAT)

C-    Get the size of a disk file connected to Fortran unit LUN

#ifdef CERNLIB_QFAPOGEE

C-    This is the method we have to use for Apogee compilers

      INTEGER APOFSZ

      INQUIRE(UNIT=LUN,NAME=FNAME)
      ISTAT = APOFSZ(FNAME,NBYTES)
      IF (ISTAT.NE.0) RETURN
#else
      INTEGER      FSEEK, FTELL

      ISTAT = FSEEK (LUN, 0, 2)
      IF (ISTAT.NE.0)              RETURN

      NBYTES = FTELL (LUN)
#endif
      IF (NBYTES.LT.0)  THEN
          ISTAT = -NBYTES
          RETURN
        ENDIF

      NWORDS = NBYTES / 4
      JREC   = NWORDS / NWREC
      RETURN
      END
