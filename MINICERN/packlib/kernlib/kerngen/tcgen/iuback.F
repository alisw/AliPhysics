*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUBACK (IC,JL,JR)
C
C CERN PROGLIB# M410    IUBACK          .VERSION KERNFOR  4.08  840613
C ORIG. 27/05/67
C
      COMMON /SLATE/ ND,JEND, DUMMY(38)
      DIMENSION IC(80),NUM(11)
      DATA  NUM   /4H    ,4H0   ,4H1   ,4H2   ,4H3   ,4H4   ,4H5   ,
     +                           4H6   ,4H7   ,4H8   ,4H9    /
C
C
      ND  = 0
      JEND= JR
      IVAL= 0
      IPOW= 1
C
   12 IF (JEND.LT.JL)  GO TO 20
      K= IUCOMP (IC(JEND),NUM,11) - 1
      IF (K) 20,19,14
   14 IVAL= (K-1)*IPOW + IVAL
      IPOW= 10*IPOW
      ND= ND + 1
   19 JEND= JEND - 1
      GO TO 12
C
   20 IUBACK= IVAL
      RETURN
      END
