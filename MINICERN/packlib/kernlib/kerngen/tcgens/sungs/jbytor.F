*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      FUNCTION JBYTOR (MZ,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    JBYTOR          .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
      JBYTOR = OR  (MZ, IBITS(IZW,IZP-1,NZB))
      RETURN
      END
