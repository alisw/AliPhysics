*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:00  mclareni
* Kernlib
*
*
      SUBROUTINE KERNGT (LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNOS9  1.01  940718
C ORIG.  16/03/89
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF
 
      COMMON /SLATE/ MCH(40)
 
 
      CALL UCTOH1   (
     + ' KERNOS9  1.01  940929 16.30            '
     +,  MCH, 40)
      CALL GSHOWF (LUN)
      RETURN
      END
