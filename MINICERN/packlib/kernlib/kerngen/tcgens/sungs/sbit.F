*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:21  mclareni
* Kernlib
*
*
      SUBROUTINE SBIT (IT,IZW,IZP)
C
C CERN PROGLIB# M421    SBIT            .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C

      IZW = OR (AND(NOT(LSHIFT(1,IZP-1)),IZW),
     +          LSHIFT(AND(IT,1),IZP-1))
      RETURN
      END
