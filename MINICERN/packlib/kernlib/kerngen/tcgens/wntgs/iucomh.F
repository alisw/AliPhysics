*
* $Id$
*
* $Log$
* Revision 1.2  1997/09/02 14:27:05  mclareni
* WINNT correction
*
* Revision 1.1  1997/02/04 17:36:29  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1  1996/02/15 17:50:24  mclareni
* Kernlib
*
*
      FUNCTION IUCOMH (ISTR1, ISTR2, NCH)
C
C CERN PROGLIB#         IUCOMH          .VERSION KERNDOS  1.00  920624
C ORIG.  FCA+JZ, 20/02/90
C
C Compare the two Hollerith strings of NCH bytes length and return
C         -1  if  STRING1  < STRING2
C          0               =
C         +1               >
C which is determined by the first non-identical character.
C
      DIMENSION ISTR1(9), ISTR2(9)
C
      NDO = NCH
      IF (NDO.LE.0)                GO TO 30
      NWD = (NDO-1)/4 + 1

      DO 29  JW=1,NWD
      IWD1 = ISTR1(JW)
      IWD2 = ISTR2(JW)
      NL   = MIN(NDO,4)
      IF (IWD1. NE. IWD2)    THEN

        DO 28  J=1,NL
#ifndef CERNLIB_QF2C
          ICHAR1 = IAND (IWD1,255)
          ICHAR2 = IAND (IWD2,255)
#else
          ICHAR1 = AND (IWD1,255)
          ICHAR2 = AND (IWD2,255)
#endif
          IF (ICHAR1-ICHAR2)           32, 26, 31
#ifndef CERNLIB_QF2C
   26     IWD1 = ISHFT (IWD1,-8)
   28   IWD2 = ISHFT (IWD2,-8)
#else
   26   IWD1 = rshift (IWD1,8)
   28   IWD2 = rshift (IWD2,8)
#endif
      END IF
   29 NDO  = NDO - 4

   30 IUCOMH = 0
      RETURN

   31 IUCOMH = 1
      RETURN
   32 IUCOMH = -1
      RETURN
      END
