*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:39  mclareni
* Kernlib
*
*
      LOGICAL FUNCTION INTRAC(DUMMY)
C
C CERN PROGLIB# Z044    INTRAC          .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
C-    The program is considered interactive if standard input
C-    is connected to the terminal

      LOGICAL      ISATTY

      INTRAC = ISATTY (5)
      RETURN
      END
