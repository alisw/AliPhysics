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
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNHPX  1.03  940511
C ORIG.  16/03/89
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF

      COMMON /SLATE/ MCH(40)
      CHARACTER    VIDQQ*(*)
      PARAMETER   (VIDQQ = '@(#)' //
     + ' KERNHPX  1.04  1995/09/28 19.17        >')

      CALL CLEFT (VIDQQ,1,0)
      CALL UCTOH1   (
     + 'KERNHPX  1.04  1995/09/28 19.17         '
     +,  MCH, 40)
      CALL GSHOWF (LUN)
      RETURN
      END
