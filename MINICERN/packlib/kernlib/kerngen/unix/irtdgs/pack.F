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
         SUBROUTINE PACK(A,B,NDIM)
         INTEGER A(2*NDIM),B(2*NDIM)
         B(1) = A(1)
         DO I=3,2*NDIM,2
            IX2 = I
            IX1 = (I+1)/2
            B(IX1) = A(IX2)
            B(IX2) = 0
         ENDDO
         END
