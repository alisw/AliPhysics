*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:29  mclareni
* Kernlib
*
*
      INTEGER FUNCTION RENAMEF (FROM,TO)
C
C CERN PROGLIB# Z265    RENAMEF         .VERSION KERNVAX  2.36  910306
C ORIG. 06/03/91, JZ
C Fortran interface routine to rename
C
      CHARACTER*(*)  FROM, TO

      LGFR = LNBLNK (FROM)
      LGTO = LNBLNK (TO)
      ISTAT = LIB$RENAME_FILE (FROM(1:LGFR), TO(1:LGTO))

      IF (ISTAT)  ISTAT = 0
      RENAMEF = ISTAT
      RETURN
      END
