*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:21  mclareni
* Kernlib
*
*
      FUNCTION MSBIT (MZ,IZW,IZP)
C
C CERN PROGLIB# M421    MSBIT           .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C

      MSBIT = OR (AND (IZW, NOT(LSHIFT(1,IZP-1)) )
     +          , LSHIFT(AND(MZ,1),IZP-1) )
      RETURN
      END
