*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:22  mclareni
* Kernlib
*
*
@PROCESS AUTODBL(NONE)
         SUBROUTINE IE3FOS(VINT,VIEEE,NV,JCODE)
         DOUBLE PRECISION VINT(NV)
         REAL VIEEE(2*NV)
         DO I=1,NV
          I2 = 2*I - 1
          VIEEE(I2) = SNGL(VINT(I))
         ENDDO
         END
