*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      FUNCTION LOCBYT (IT,IA,N,NW,LA,NBITS)
C
C CERN PROGLIB# M428    LOCBYT          .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C
      DIMENSION IA(9)
C
C
      IF (N.EQ.0)  GO TO 19
#if defined(CERNLIB_BUGLRSHFT)
      MASK = LSHIFT (ishft(-1,-32+NBITS),LA-1)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MASK = LSHIFT (lrshft(-1,32-NBITS),LA-1)
#endif
      ITSH = LSHIFT (IT,LA-1)

      DO 18 J=1,N,NW
      IF (AND(IA(J),MASK).EQ.ITSH)  GO TO 21
   18 CONTINUE
   19 LOCBYT= 0
      RETURN
C
   21 LOCBYT= J
      RETURN
      END
