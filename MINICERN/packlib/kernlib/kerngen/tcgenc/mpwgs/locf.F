*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:40  mclareni
* Kernlib
*
*
      FUNCTION LOCF (AM)
C
C CERN PROGLIB# N100    LOCF            .VERSION KERNMPW  1.00  891208
C ORIG. 13/12/83  F.Carminati
C
      DIMENSION AM(*)
 
      LOCF= %LOC (AM(1))
      LOCF= ISHFT(LOCF, -2)
      END
