*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:22  mclareni
* Kernlib
*
*
@PROCESS AUTODBL(NONE)
         SUBROUTINE IE3TOS(VIEEE,VINT,NV,JCODE)
         DOUBLE PRECISION VINT(NV)
         REAL VIEEE(2*NV)
         DO I=1,NV
          VINT(I) = DBLE(VIEEE((I-1)*2+1))
         ENDDO
         END
