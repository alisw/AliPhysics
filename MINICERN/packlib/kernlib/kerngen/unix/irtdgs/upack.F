*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:07  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
@PROCESS AUTODBL(NONE)
         SUBROUTINE UPACK(A,B,NDIM)
*
* Unpack 32 bit words to 64 bit words padding rightmost significant
* bit with 0. Should work whatever A and B adresses provided.
* MR/IBM, 8-9-1993
*
         INTEGER A(2*NDIM),B(2*NDIM)
*
         IADA = LOCB(A)
         IADB = LOCB(B)
         NS   = ISHFT (LOCB(A)-LOCB(B),-2)

         DO I=1,NS
            IX2 = I
            IX1 = 2*I -1
            B(IX1)   = A(IX2)
         ENDDO

         DO I=NDIM,NS+1,-1
            IX1 = I
            IX2 = 2*I -1
            B(IX2) = A(IX1)
         ENDDO

         DO I=2,2*NDIM,2
            B(I) = 0
         ENDDO

         END
