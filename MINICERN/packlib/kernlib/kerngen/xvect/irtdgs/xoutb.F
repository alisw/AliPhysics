*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
      SUBROUTINE XOUTB  (LUN,BV,NB)
C
C CERN PROGLIB# Z203    XOUTB           .VERSION KERNIRT  1.06  930811
C ORIG. 15/07/78 JZ
C
      DIMENSION    BV(NB)
C
      REAL XNB
      INTEGER INB
      EQUIVALENCE (XNB,INB)
      INB = NB
      WRITE (LUN)  XNB,BV
      RETURN
      END
