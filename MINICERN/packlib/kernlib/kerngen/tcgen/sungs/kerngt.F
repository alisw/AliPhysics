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
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNSUN  1.06  920511
C ORIG. 16/03/89  JZ
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF
C
      COMMON /SLATE/ MCH(40)
      CHARACTER    VIDQQ*(*)
      PARAMETER   (VIDQQ = '@(#)' //
     + 'KERNSUN  1.11  1995/09/28 19.14         >')

      CALL CLEFT (VIDQQ,1,0)
      CALL UCTOH1 (
     + 'KERNSUN  1.11  1995/09/28 19.14         '
     +,  MCH,40)
      CALL GSHOWF (LUN)
      RETURN
      END
