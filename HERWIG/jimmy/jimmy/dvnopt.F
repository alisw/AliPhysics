*CMZ :          05/01/95  11.18.45  by  Jonathan Butterworth
*-- Author :
      SUBROUTINE DVNOPT
C --------------------------------------------------------------------------
C      CERNLIB INTEGRATION OPTIONS
C --------------------------------------------------------------------------

#include "jimmy.inc"


      COMMON/QUADRE/IDEG
      COMMON/PRINT/IPRINT
      COMMON/BNDLMT/FLOW,FHIGH
C
      IDEG   = 1
C
      IF (JMBUG.LT.1) THEN
        IPRINT = 0
      ELSE
        IPRINT = 1
      ENDIF
C
      FLOW = 0.0
C
      RETURN
      END
