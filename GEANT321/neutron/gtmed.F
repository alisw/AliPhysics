*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE GTMED(MEDGEA,MEDMOR)
#include "geant321/mmicap.inc"
#include "geant321/minput.inc"
C get MICAP material number
      DO 10 I=1,MEDIA
         IF(LD(LGE2MO+I).EQ.MEDGEA) THEN
            MEDMOR = I
            GOTO 20
         ENDIF
   10 CONTINUE
      WRITE(IOUT,'('' MICAP GTMED: GEANT Medium '',I5,
     +             '' not found ==> STOP'')') MEDGEA
      STOP
   20 RETURN
      END
