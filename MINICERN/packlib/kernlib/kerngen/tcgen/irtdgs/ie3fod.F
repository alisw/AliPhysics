*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:57  mclareni
* Kernlib
*
*
@PROCESS AUTODBL(NONE)
         SUBROUTINE IE3FOD(VINT,VIEEE,NV,JCODE)
         INTEGER VINT(4*NV),VIEEE(4*NV)
         DO I=1,NV
          I1 = (I-1)*4 + 1
          I2 = (I-1)*4 + 2
          I3 = (I-1)*4 + 3
          VIEEE(I1) = VINT(I1)
          VIEEE(I3) = VINT(I2)
         ENDDO
         END
