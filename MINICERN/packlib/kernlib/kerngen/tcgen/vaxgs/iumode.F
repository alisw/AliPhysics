*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:59  mclareni
* Kernlib
*
*
      FUNCTION IUMODE (IWP)
C
C CERN PROGLIB# M506    IUMODE          .VERSION KERNVAX  2.09  840215
C ORIG. 7/12/79 JZ
C
      COMMON /CIUMOD/ JMODE, JAMBIG
      DIMENSION    IWP(9)

#include "kerngen/q_jbit.inc"

      IWORD = IWP(1)
      JAMBIG= 0
      IF (IWORD.EQ.0)        GO TO 39

      JL = JBYT (IWORD,25,8)
      IF (JMODE.EQ.0)        GO TO 31

C--                EXPECT FLOATING

      IF (JL.EQ.0)           GO TO 24
      IF (JL.NE.255)         GO TO 29
   24 J = JBYT (IWORD,8,8)
      IF (J.EQ.0)            GO TO 39
      IF (J.EQ.255)          GO TO 39
      JAMBIG = 1
   29 IUMODE = 1
      RETURN

C--                EXPECT INTEGER

   31 IF (JL.EQ.0)           GO TO 32
      IF (JL.NE.255)         GO TO 29
   32 J = JBYT (IWORD,8,8)
      IF (J.EQ.0)            GO TO 39
      IF (J.EQ.255)          GO TO 39
      JAMBIG = 1
   39 IUMODE = 0
      RETURN
      END
