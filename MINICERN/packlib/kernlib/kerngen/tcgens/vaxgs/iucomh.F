*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:23  mclareni
* Kernlib
*
*
      FUNCTION IUCOMH (ISTR1, ISTR2, NCH)
C
C CERN PROGLIB#         IUCOMH          .VERSION KERNVAX  2.39  930112
C ORIG.  FCA+JZ, 20/02/90
C
C Compare the two Hollerith strings of NCH bytes length and return
C         -1  if  STRING1  < STRING2
C          0               =
C         +1               >
C which is determined by the first non-identical character.
C
      DIMENSION ISTR1(9), ISTR2(9)
      PARAMETER (NCHAPW=4,NBITPC=8)
C
      NDO = NCH
      IF (NDO.LE.0)                GO TO 30
      NWD = (NDO-1)/NCHAPW + 1

      DO 29  JW=1,NWD
      IWD1 = ISTR1(JW)
      IWD2 = ISTR2(JW)
      NL   = MIN(NDO,4)

      DO 28  J=1,NL
      ICHAR1 = IAND (IWD1,255)
      ICHAR2 = IAND (IWD2,255)
      IF (ICHAR1-ICHAR2)           32, 26, 31
   26 IWD1 = ISHFT (IWD1,-8)
   28 IWD2 = ISHFT (IWD2,-8)
   29 NDO  = NDO - 4

   30 IUCOMH = 0
      RETURN

   31 IUCOMH = 1
      RETURN
   32 IUCOMH = -1
      RETURN
      END
