*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE COMPO
C
C *** GET PARAMETERS FOR THE TARGET NUCLEUS ***
C *** NVE 11-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (07-DEC-1984)
C
#include "geant321/mxgkgh.inc"
#include "geant321/s_curpar.inc"
#include "geant321/s_prntfl.inc"
C --- GEANT COMMONS ---
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmate.inc"
      DIMENSION RNDM(1)
C
C --- CHECK FOR COMPOUND ---
      KK=ABS(Q(JMA+11))
      IF (KK .GE. 2) GO TO 10
C
C --- ELEMENTS ---
      ATNO2=A
      ZNO2=Z
      IF (NPRT(4)) WRITE(NEWBCD,1000) ATNO2,ZNO2
 1000 FORMAT(' *COMPO* INTERACTION ON NUCLEUS OF ELEMENT WITH',
     $ ' A = ',G12.5,' Z = ',G12.5)
      GO TO 9999
C
C --- COMPOUNDS ===> SELECT NUCLEUS ---
 10   CONTINUE
C
      SUM=0.0
      DO 11 I=1,KK
      AI=Q(JMIXT+I)
      ZI=Q(JMIXT+KK+I)
      WI=Q(JMIXT+2*KK+I)
      SUM=SUM+WI/AI
 11   CONTINUE
      CALL GRNDM(RNDM,1)
      TEST1=RNDM(1)*SUM
C
      TEST2=0.0
      DO 12 I=1,KK
      JCOMPO=I
      AI=Q(JMIXT+I)
      ZI=Q(JMIXT+KK+I)
      WI=Q(JMIXT+2*KK+I)
      TEST2=TEST2+WI/AI
      IF (TEST2 .GT. TEST1) GO TO 20
 12   CONTINUE
C
 20   CONTINUE
      ATNO2=AI
      ZNO2=ZI
      IF (NPRT(4)) WRITE(NEWBCD,1001) JCOMPO,ATNO2,ZNO2
 1001 FORMAT(' *COMPO* COMPONENT NO. ',I2,' HAS BEEN SELECTED WITH',
     $ ' A = ',G12.5,' Z = ',G12.5)
C
 9999 CONTINUE
      END
