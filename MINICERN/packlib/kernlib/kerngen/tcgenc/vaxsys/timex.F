*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:41  mclareni
* Kernlib
*
*
      SUBROUTINE TIMEX (T)
C CERN PROGLIB# Z007    TIMEX           .VERSION KERNVAX  2.45  940228
C ORIG.  J. Shiers, 15/02/94
C
      INCLUDE '($JPIDEF)'


      ISTAT = LIB$GETJPI(JPI$_CPUTIM,,,ITIME)
      T     = FLOAT(ITIME)/100.
      RETURN
      END
