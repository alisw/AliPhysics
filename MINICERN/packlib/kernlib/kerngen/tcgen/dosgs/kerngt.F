*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:57  mclareni
* Kernlib
*
*
      SUBROUTINE KERNGT (LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNDOS  1.00  920624
C ORIG.  16/03/89
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF

      COMMON /SLATE/ MCH(40)


      CALL UCTOH1   (
     + ' KERNDOS  1.01/5 940309 10.40           '
     +,  MCH, 40)
      CALL GSHOWF (LUN)
      RETURN
      END
