*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:29  mclareni
* Kernlib
*
*
      INTEGER FUNCTION SYSTEMF (NAME)
C
C CERN PROGLIB# Z265    SYSTEMF         .VERSION KERNVAX  2.46  940913
C ORIG. 06/03/91, JZ,  mod. march 94, JS
C Fortran interface routine to system
C
      CHARACTER    NAME*(*)
      INCLUDE      '($SSDEF)'

      LGN   = LNBLNK (NAME)
      ISTAT = LIB$SPAWN (NAME(1:LGN))
      IF(ISTAT.EQ.SS$_NORMAL) ISTAT = 0
      SYSTEMF = ISTAT
      RETURN
      END
