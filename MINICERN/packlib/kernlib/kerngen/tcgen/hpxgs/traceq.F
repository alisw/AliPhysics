*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:57  mclareni
* Kernlib
*
*
      SUBROUTINE TRACEQ (LUNP,LEV)
C
C CERN PROGLIB# N105    TRACEQ          .VERSION KERNHPX  1.04  950928
C ORIG.  3/05/95  FR, JZ
C
*     unit number for standard output
      PARAMETER    (LUNSTP=6)

      LUN = LUNP
      IF (LUN.EQ.0)  LUN = LUNSTP
      WRITE (LUN,9001)
 9001 FORMAT (/' TRACEQ.  Trace-back printed to standard output.')

      CALL TRACEQC
      RETURN
      END
