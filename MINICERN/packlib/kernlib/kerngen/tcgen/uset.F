*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE USET (INT,IC,JL,JR)
C
C CERN PROGLIB# M410    USET            .VERSION KERNFOR  4.08  840613
C ORIG. 15/03/68
C
      COMMON/SLATE/ND,JEND,DUMMY(38)
      DIMENSION IC(80),NUM(11)
      DATA  NUM   /4H    ,4H0   ,4H1   ,4H2   ,4H3   ,4H4   ,4H5   ,
     +                           4H6   ,4H7   ,4H8   ,4H9    /
C
      ND=0
      JEND=JR
      IVAL=INT
C
   12 IF (JEND.LT.JL) RETURN
      K=MOD (IVAL,10)
      IC(JEND)=NUM(K+2)
      ND=ND+1
      JEND=JEND-1
      IVAL=IVAL/10
      IF (IVAL.GE.1) GO TO 12
      RETURN
      END
