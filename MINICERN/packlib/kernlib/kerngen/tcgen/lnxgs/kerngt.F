*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      SUBROUTINE KERNGT (LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNLNX  1.00  930829
C ORIG.  16/03/89
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF

      COMMON /SLATE/ MCH(40)


      CALL UCTOH1   (
     + 'KERNLNX  1.03  1995/07/18 13.54           '       
     +,  MCH, 40)
      CALL GSHOWF (LUN)
      RETURN
      END
