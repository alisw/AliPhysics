*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GCKMAT (ITMED,NATMED)
*
*
*    This routine check the material table coming from a Cadd system.
*    It may help to take default/savety actions.
*
*
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
*
      CHARACTER*20  NAME,NATMED
*
* *** check tracking media number
      JTM = LQ(JTMED - ITMED)
      IF (ITMED.GT.NTMED.OR.JTM.LE.0) THEN
         WRITE (CHMAIL,10200) ITMED
         CALL GMAIL (0,0)
         IMATE = 1
         CALL GSTMED (ITMED,NATMED,IMATE, 0, 0, 1., 10.,
     +                100., 0.2, 0.005, 0.001, 0, 0)
      ELSE
         CALL UHTOC (IQ(JTM+1),4,NAME,20)
         IF (NATMED.NE.NAME) THEN
            WRITE (CHMAIL,10300) ITMED
            CALL GMAIL (0,0)
         ENDIF
      ENDIF
*
10000 FORMAT(' *** GCKMAT *** Warning, material number :', I3,
     +       ' does not exist; force by default IMATE = 1' )
*
10100 FORMAT(' *** GCKMAT *** Warning, material number :', I3,
     +      ' maybe not egal to the Geant definition. KEEP the old one')
*
10200 FORMAT(' *** GCKMAT *** Warning, tracking media  :', I3,
     +    ' does not exist; It will be created with default parameters')
*
10300 FORMAT(' *** GCKMAT *** Warning, tracking media  :', I3,
     +      ' maybe not egal to the Geant definition. KEEP the old one')
*
      END
