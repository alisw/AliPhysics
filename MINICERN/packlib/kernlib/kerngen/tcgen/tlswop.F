*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:54  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TLSWOP (A,AD,N,NR)
C
C CERN PROGLIB# E230    TLSWOP          .VERSION KERNFOR  2.06  740511
C ORIG. 11/05/74 WH+WM
C
C.  SUBROUTINE TLSWOP (A,AD,N,NR)
C.
C.       INTERCHANGES THE FIRST ROW  OF MATRICES A AND AD.
C.
C.-------------------------------------------------------------------
C
      COMMON /SLATE/ DUM(37),H,I,JA
      DIMENSION      A(*), AD(*)
C
      JA = 1
      DO           10        I=1,NR
      H  = A(JA)
      A(JA)  = AD(JA)
      AD(JA) = H
   10 JA = JA + N
C
      RETURN
      END
