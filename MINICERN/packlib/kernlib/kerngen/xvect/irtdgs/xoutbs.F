*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
      SUBROUTINE XOUTBS (LUN,AV,NA,BV,NB)
C
C CERN PROGLIB# Z203    XOUTBS          .VERSION KERNIRT  1.06  930811
C ORIG. 15/07/78 JZ
C
      DIMENSION    BV(NB), AV(NA)
C
      REAL XNB
      INTEGER INB
      EQUIVALENCE (XNB,INB)
      INB = NB
      WRITE (LUN)  XNB,AV,BV
      RETURN
      END
