*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      FUNCTION INCBYT (INCP,IA,JTHP,NBITS)
C
C CERN PROGLIB# M423    INCBYT          .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C
      DIMENSION INCP(9),IA(99),JTHP(9),NBITS(3)

      INC = INCP(1)
      JTH = JTHP(1)
      NZB = NBITS(1)
      NPW = NBITS(2)
      NPB = NBITS(3)

      INCBYT = 0

      JWORD= (JTH-1)/NPW + 1
      JPOS = JTH - (JWORD-1)*NPW
      NL   = (JPOS-1)*NZB
#if defined(CERNLIB_BUGLRSHFT)
      MASK = LSHIFT (ishft(-1,NZB-32),NL)
      IZW  = IA(JWORD)
      NEW  = INC + ishft (AND(MASK,IZW), -NL)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MASK = LSHIFT (lrshft(-1,32-NZB),NL)
      IZW  = IA(JWORD)
      NEW  = INC + lrshft (AND(MASK,IZW), NL)
#endif
      IF (NEW.GT.NPB)        GO TO 21
   17 IA(JWORD) = OR (AND(NOT(MASK),IZW), LSHIFT(NEW,NL))
      RETURN
C
   21 INCBYT = NEW - NPB
      IF (INCBYT.EQ.INC)     RETURN
      NEW = NPB
      GO TO 17
      END
