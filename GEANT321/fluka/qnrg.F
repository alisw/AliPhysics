*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE QNRG.FOR
*COPY QNRG
*                                                                      *
*=== qnrg =============================================================*
*                                                                      *
      FUNCTION QNRG ( A1, Z1, A2, Z2 )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C--------------------------------------------------------------------
C SUBNAME = QNRG
C--------------------------------------------------------------------
#include "geant321/eva0.inc"
C--------------------------------------------------------------------
      KA1  = NINT (A1)
      KZ1  = NINT (Z1)
      KA2  = NINT (A2)
      KZ2  = NINT (Z2)
      N1= KA1 - KZ1
      N2= KA2 - KZ2
      IF (N1 .LE. 0) STOP 105
      IF (N2 .LE. 0) STOP 106
      IZZ1 = NINT ( WAPS (KA1,1) )
      IZZ2 = NINT ( WAPS (KA2,1) )
      IF ( KZ1 .LT. IZZ1 .OR. KZ1 .GT. IZZ1 + 9 ) THEN
         ENRG1 = ENRG ( A1, Z1 )
      ELSE
         IZ1   = KZ1 - IZZ1 + 2
         ENRG1 = WAPS ( KA1, IZ1 )
         IF ( ENRG1 .EQ. 0.D0 .AND. ( KA1 .NE. 12 .OR. KZ1 .NE. 6) )
     &        THEN
            ENRG1 = ENRG ( A1, Z1 )
         END IF
      END IF
      IF ( KZ2 .LT. IZZ2 .OR. KZ2 .GT. IZZ2 + 9 ) THEN
         ENRG2 = ENRG ( A2, Z2 )
      ELSE
         IZ2   = KZ2 - IZZ2 + 2
         ENRG2 = WAPS ( KA2, IZ2 )
         IF ( ENRG2 .EQ. 0.D0 .AND. ( KA2 .NE. 12 .OR. KZ2 .NE. 6) )
     &        THEN
            ENRG2 = ENRG ( A2, Z2 )
         END IF
      END IF
      QNRG = ENRG1 - ENRG2
      RETURN
      END
