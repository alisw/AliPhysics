*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:41  mclareni
* Kernlib
*
*
      SUBROUTINE TIMED (T)
C CERN PROGLIB# Z007    TIMED           .VERSION KERNVAX  2.45  940228
C ORIG.  J. Shiers, 15/02/94
C
      INCLUDE '($JPIDEF)'

      DATA  ITOLD  / 0 /
      SAVE  ITOLD


      ISTAT = LIB$GETJPI(JPI$_CPUTIM,,,ITIME)
      T     = FLOAT(ITIME-ITOLD)/100.
      ITOLD = ITIME
      RETURN
      END
