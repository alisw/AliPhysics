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
*
      SUBROUTINE GDSHIF(IVOLNA,ICPOIN)
*
*     This subroutine allows to shift volumes in a more visible
*     place just for drawing, without alterating the data structure,
*     but perfectly simulating it.
*
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcgobj.inc"
*SG
#include "geant321/gcmutr.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gchiln.inc"
#include "geant321/gcspee.inc"
*SG
      COMMON /QUEST/IQUEST(100)
*SG
      SAVE IOLDVO
      DATA IOLDVO/0/
*
*    Shifting object
         IF(KSHIFT.EQ.2)THEN
           DO 333 IJ=1,100
             IF(IVOLNA.EQ.IVECVO(IJ))GOTO 554
             IF(IVECVO(IJ).EQ.0)THEN
               IVECVO(IJ)=IVOLNA
               GOTO 334
             ENDIF
 333       CONTINUE
         ENDIF
 334     CONTINUE
         IF(KSHIFT.GT.0)THEN
            IF(ISUBLI.EQ.1.AND.NIET.EQ.1)THEN
             NIET=1
             CALL CGSHIF(PORGX,PORGY,PORGZ,
     +       Q(ICPOIN))
             IF(IVOLNA.EQ.IOLDVO)GOTO 554
             GOTO 54
            ENDIF
            IF(ISUBLI.EQ.1.AND.NIET.EQ.2.AND.IVOOLD.NE.0)THEN
             NIET=1
             IF(IVOLNA.EQ.IOLDVO)THEN
              CALL CGSHIF(PORGX,PORGY,PORGZ,
     +        Q(ICPOIN))
              GOTO 554
             ENDIF
             DO 77 I=1,15
            IF((NLEVEL-I).GT.0)THEN
             IF(POX(NLEVEL-I).NE.0.OR.POY(NLEVEL-I).NE.0.OR.
     +       POZ(NLEVEL-I).NE.0)THEN
              PORGX=POX(NLEVEL-I)
              PORGY=POY(NLEVEL-I)
              PORGZ=POZ(NLEVEL-I)
              GOTO 777
             ENDIF
            ELSE
               PORGX=0.
               PORGY=0.
               PORGZ=0.
               POX(NLEVEL)=0.
               POY(NLEVEL)=0.
               POZ(NLEVEL)=0.
               GOTO 777
            ENDIF
 77          CONTINUE
 777         CONTINUE
             CALL CGSHIF(PORGX,PORGY,PORGZ,
     +       Q(ICPOIN))
             GOTO 54
            ENDIF
             IF(NIET.EQ.0.AND.IVOLNA.EQ.IOLDVO)GOTO 554
             NIET=0
 54          CONTINUE
             DO 553 KHHH=MULTRA,1,-1
                CALL UCTOH(GNVNV(KHHH),IVN,4,4)
                IF(IVN.EQ.IVOLNA)THEN
                   NIET=1
                   CALL CGSHIF(GXXXX(KHHH),GYYYY(KHHH),GZZZZ(KHHH),
     +             Q(ICPOIN))
                   IF(IVOLNA.NE.IOLDVO)THEN
                    POX(NLEVEL)=PORGX+GXXXX(KHHH)
                    POY(NLEVEL)=PORGY+GYYYY(KHHH)
                    POZ(NLEVEL)=PORGZ+GZZZZ(KHHH)
                    PORGX=PORGX+GXXXX(KHHH)
                    PORGY=PORGY+GYYYY(KHHH)
                    PORGZ=PORGZ+GZZZZ(KHHH)
                   ENDIF
                   IVOOLD=IVOLNA
                   IOLDVO=IVOLNA
                   GOTO 554
                ENDIF
  553        CONTINUE
             IOLDVO=IVOLNA
         ENDIF
  554    CONTINUE
  999    END
 
 
