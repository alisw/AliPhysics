*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:13  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE PKCHAR (IGV,ICHV,NN,IPAR)
C
C CERN PROGLIB# M427    PKCHAR          .VERSION KERNFOR  4.20  881028
C ORIG. 07/07/72 JZ
C
      COMMON /SLATE/ JCH,N,NBITS,NCHAR,NZONE,IGNOR,JIN,NBDONE,NBEND
     +,              JGO,JBT,NLEFT,NRIGH,DUMMY(27)
      DIMENSION ICHV(*),IGV(*),IPAR(9)
C                .
C
#include "kerngen/wordsize.inc"
C
C
      N = NN
      IF (N.EQ.0)  RETURN
      NBITS = IPAR(1)
      NCHAR = IPAR(2)
      NZONE = IPAR(3)
      IGNOR = IPAR(4)
      NFILL = IPAR(5)
      JIN = 1
      JCH = 1
      ICHV(1)= NFILL
      IF (NZONE.NE.0)  GO TO 41
C
C----              ZONE = WORD
C
   21 JGO  = NBITPW+1 - IGNOR
      NBEND= NBITS + 1
      IF (NCHAR.EQ.0)  GO TO 22
      NBEND = MAX (NBEND,JGO-NBITS*(NCHAR-1))
C
   22 JBT = JGO
C
   24 JBT = JBT - NBITS
      CALL SBYT (IGV(JIN),ICHV(JCH),JBT,NBITS)
      IF (JIN.EQ.N)  RETURN
      JIN = JIN + 1
      IF (JBT.GE.NBEND)  GO TO 24
      JCH = JCH + 1
      ICHV(JCH)= NFILL
      GO TO 22
C
C----              ZONE NOT EQUALS A WORD
C
   41 IF (NZONE.EQ.NBITPW)  GO TO 21
      JBT   = NBITPW - IGNOR
      NBEND = NZONE+1 - NBITS
      IF (NCHAR.EQ.0)  GO TO 42
      NBEND = MIN (NBEND, IGNOR+NCHAR*NBITS)
C
   42 NBDONE= IGNOR
   43 IF (JBT.GE.0)  GO TO 44
      JBT= JBT + NBITPW
      JCH= JCH + 1
      ICHV(JCH)= NFILL
      GO TO 43
C
   44 JBT = JBT - NBITS
      IF (JBT.LT.0)  GO TO 51
      CALL SBYT (IGV(JIN),ICHV(JCH),JBT+1,NBITS)
      GO TO 58
C
C--                CHARACTER ACROSS WORD BOUNDARY
   51 NLEFT= NBITS + JBT
      NRIGH= -JBT
      JBT  = JBT + NBITPW
      IF (NLEFT.EQ.0)  GO TO 54
      CALL CBYT (IGV(JIN),NRIGH+1,ICHV(JCH),1,NLEFT)
   54 JCH= JCH + 1
      ICHV(JCH)= NFILL
      CALL SBYT (IGV(JIN),ICHV(JCH),JBT+1,NRIGH)
C
   58 IF (JIN.EQ.N)  RETURN
      JIN= JIN + 1
      NBDONE= NBDONE + NBITS
      IF (NBDONE.LT.NBEND)  GO TO 44
      JBT = JBT - (NZONE-NBDONE) - IGNOR
      GO TO 42
      END
