*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      FUNCTION JBIT (IZW,IZP)
C
C CERN PROGLIB# M421    JBIT            .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
      JBIT = AND (RSHIFT(IZW,IZP-1), 1)
      RETURN
      END
