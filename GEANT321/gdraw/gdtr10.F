*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR10(JNODE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Plot tree structure from JNODE                           *
C.    *                                                                *
C.    *       JNODE = starting node of tree scanning (input)           *
C.    *                                                                *
C.    *    ==>Called by : GDTREE                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
C.
      J=JNODE
      IDIR=1
   10 CONTINUE
      IF (IDIR.EQ.1) THEN
         CALL GDTR11(J,IDIR,JNODE)
         GO TO 10
      ELSE IF (IDIR.EQ.2) THEN
   20    JB=IQ(JBRO1+J)
         IF (JB.EQ.0) THEN
            J=IQ(JMOT1+J)
            IF (J.EQ.0) GO TO 999
            LEVVER=LEVVER-1
            JJ=-J
            CALL GDTR12(JJ)
            IF (J.EQ.JNODE) GO TO 999
            GO TO 20
         ELSE
            J=JB
            CALL GDTR11(J,IDIR,JNODE)
            GO TO 10
         ENDIF
      ELSE
         GO TO 999
      ENDIF
  999 RETURN
      END
