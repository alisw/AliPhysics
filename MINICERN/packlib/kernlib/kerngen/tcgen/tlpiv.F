*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLPIV(A,B,IASEP,IBSEP,M,H,G)
C
C CERN PROGLIB# E230    TLPIV           .VERSION KERNFOR  2.06  740511
C ORIG. 11/05/74 WH+WM
C
C.  SUBROUTINE TLPIV (A,B,IASEP,IBSEP,M,H,G)
C.
C.        CALCULATES PIVOTING FACTORS FOR TLS
C.
C.-------------------------------------------------------------------
C
      COMMON /SLATE/ DUM(37),I,JA,JB
      DIMENSION A(*),B(*)
C
      H=0.0
      G=0.0
      JA=1
      JB=1
      DO 2 I=1,M
      G=G+A(JA)*B(JB)
      H=H+A(JA)*A(JA)
      JA=JA+IASEP
    2 JB=JB+IBSEP
      RETURN
      END
