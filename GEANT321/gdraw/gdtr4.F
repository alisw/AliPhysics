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
      SUBROUTINE GDTR4
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Collect multiplicity for all JVOLUM structure            *
C.    *                                                                *
C.    *    ==>Called by : GDTR0                                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
      DO 10 J=1,NUMNOD
         IQ(JSCA+J)=0
         IQ(JDUP+J)=0
   10 CONTINUE
      IQ(JDUP+1)=1
      DO 30 J=1,NUMNOD
         IF (IQ(JSCA+J).EQ.0) THEN
            CALL GDTR3(J)
            DO 20 I=J,NUMNOD
               IF (IQ(JSCA+I).EQ.-1) IQ(JSCA+I)=0
   20       CONTINUE
         ENDIF
   30 CONTINUE
      END
