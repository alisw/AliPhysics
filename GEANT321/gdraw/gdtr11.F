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
      SUBROUTINE GDTR11(J,IDIR,JNODE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Down a node on the tree                                  *
C.    *                                                                *
C.    *       J = node pointer (input/output)                          *
C.    *       IDIR = direction of search; IDIR=1 for down              *
C.    *              and IDIR=2 for right (input/output)               *
C.    *       JNODE = starting node of tree scanning (input)           *
C.    *                                                                *
C.    *    ==>Called by : GDTR10                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
C.
      JX=IQ(JXON1+J)
      IF (LEVVER.EQ.MAXV) JX=0
      IF (JX.EQ.0) THEN
         CALL GDTR12(J)
         IF (J.EQ.JNODE) THEN
            IDIR=0
         ELSE
            JB=IQ(JBRO1+J)
            IF (JB.EQ.0) THEN
               JM=IQ(JMOT1+J)
               LEVVER=LEVVER-1
               IF (JM.EQ.0) THEN
                  IDIR=0
               ELSE
                  J=JM
                  JJ=-J
                  CALL GDTR12(JJ)
                  IF (J.EQ.JNODE) THEN
                     IDIR=0
                  ELSE
                     IDIR=2
                  ENDIF
               ENDIF
            ELSE
               J=JB
               IDIR=1
            ENDIF
         ENDIF
      ELSE
         LEVVER=LEVVER+1
         MLEVV=MAX(MLEVV,LEVVER)
         J=JX
         IDIR=1
      ENDIF
      END
