*
* $Id$
*
* $Log$
* Revision 1.1  1997/02/04 17:36:21  mclareni
* Merge Winnt and 97a versions
*
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
     + ' KERNNT  1.00/0 960510 12.23 by V.Fine '
     +,  MCH, 40)
      CALL GSHOWF (LUN)
      RETURN
      END
