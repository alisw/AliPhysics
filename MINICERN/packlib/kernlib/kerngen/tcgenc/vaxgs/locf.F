*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:39  mclareni
* Kernlib
*
*
      FUNCTION LOCF (IVAR)
C
C CERN PROGLIB# N100    LOCF            .VERSION KERNVAX  2.39  930112
C ORIG. 92/07/13 Miguel Marquina
C
      DIMENSION    IVAR(9)

      J = %LOC(IVAR)
      LOCF = ISHFT (J, -2)
      END
