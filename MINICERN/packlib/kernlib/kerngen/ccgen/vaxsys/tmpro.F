*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:30  mclareni
* Kernlib
*
*
      SUBROUTINE TMPRO (TEXT)
C
C CERN PROGLIB# Z313    TMPRO           .VERSION KERNVAX  2.45  940228
C ORIG. 16/02/94, JZ
C
C     Print a prompt string to stdout
C
      CHARACTER    TEXT*(*)

      WRITE (6,9001) TEXT
 9001 FORMAT (1X,A,$)
      RETURN
      END
