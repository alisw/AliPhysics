*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE SIGFER.FOR
*COPY SIGFER
*
*=== sigfer ===========================================================*
*
      SUBROUTINE SIGFER ( KP, EKIN, POO, LFERMI )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/paprop.inc"
      LOGICAL LFERMI
*
      IF ( LFERMI ) THEN
         EKEWLL = EKIN   + VPRBIM
         EEMAX  = EKEWLL + EKFBIM + AMNUCL (ITNCMX) + AM (KP)
         PPRWLL = SQRT ( EKEWLL * ( EKEWLL + 2.D+00 * AM (KP) ) )
         IF ( PFRBIM .LT. PPRWLL ) THEN
            PPMAX  = PPRWLL + PFRBIM
            PPMIN  = PPRWLL - PFRBIM
            UMO2   = ( EEMAX + PPMAX ) * ( EEMAX - PPMAX )
            EKEMIN = 0.5D+00 * ( UMO2 - AM (KP)**2 - AMNUSQ (ITNCMX) )
     &             / AMNUCL (ITNCMX) - AM (KP)
            EKEMIN = MIN ( EKIN, EKEMIN )
            TMPEKI = 0.1666D+00 * EKIN
            EKEMIN = MAX ( EKEMIN, TMPEKI )
         ELSE
            EKEMIN = MAX ( 0.005D+00, 0.1666D+00 * EKIN )
            PPMIN  = 0.D+00
         END IF
         PPRMIN = SQRT ( EKEMIN * ( EKEMIN + 2.D+00 * AM (KP) ) )
         UMO2   = ( EEMAX + PPMIN ) * ( EEMAX - PPMIN  )
         EKEMAX =  0.5D+00 * ( UMO2 - AM (KP)**2 - AMNUSQ (ITNCMX) )
     &          / AMNUCL (ITNCMX) - AM (KP)
         PPRMAX = SQRT ( EKEMAX * ( EKEMAX + 2.D+00 * AM (KP) ) )
      ELSE
         EKEMIN = EKIN
         PPRMIN = POO
         EKEMAX = EKIN
         PPRMAX = POO
      END IF
*
  50  CONTINUE
      GO TO (  100,  200,  300,  400,  500,  600,  700,  800,  900,
     &        1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
     &        1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700,
     &        2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600,
     &        3700, 3800, 3900 ), KP
      STOP 'GEO-KP'
  100 CONTINUE
         IF ( EKEMIN .LE. 0.700D+00 ) THEN
            BETAPR = PPRMIN / ( EKEMIN + AM (KP) )
            IF ( EKEMIN .LE. 0.04D+00 ) THEN
               SIGMAN = 3.D+03 * PI / ( 1.206D+03 * EKEMIN + ( -1.86D+00
     &                + 0.09415D+03 * EKEMIN + 0.0001306D+06 * EKEMIN**2
     &                )**2 ) + 1.D+03 * PI / ( 1.206D+03 * EKEMIN
     &                + ( 0.4223D+00 + 0.13D+03 * EKEMIN )**2 )
               IF ( EKEMIN .LT. 0.02D+00 ) THEN
                  SIGMAP = 0.3333333333333333D+00 * SIGMAN
               ELSE
                  SIGMAP = 10.63D+00 / BETAPR**2 - 29.92D+00 / BETAPR
     &                   + 42.9D+00
               END IF
            ELSE
               SIGMAN = 34.10D+00 / BETAPR**2 - 82.2D+00 / BETAPR
     &                + 82.2D+00
               SIGMAP = 10.63D+00 / BETAPR**2 - 29.92D+00 / BETAPR
     &                + 42.9D+00
            END IF
         ELSE
            STOP 'Sigfer: EKE'
         END IF
      GO TO 4000
  200 CONTINUE
      GO TO 4000
  300 CONTINUE
      GO TO 4000
  400 CONTINUE
      GO TO 4000
  500 CONTINUE
      GO TO 4000
  600 CONTINUE
      GO TO 4000
  700 CONTINUE
         SIGMAN = 0.D+00
         SIGMAP = 0.D+00
      GO TO 4000
  800 CONTINUE
         IF ( EKEMIN .LE. 0.700D+00 ) THEN
            BETAPR = PPRMIN / ( EKEMIN + AM (KP) )
            IF ( EKEMIN .LE. 0.04D+00 ) THEN
               SIGMAP = 3.D+03 * PI / ( 1.206D+03 * EKEMIN + ( -1.86D+00
     &                + 0.09415D+03 * EKEMIN + 0.0001306D+06 * EKEMIN**2
     &                )**2 ) + 1.D+03 * PI / ( 1.206D+03 * EKEMIN
     &                + ( 0.4223D+00 + 0.13D+03 * EKEMIN )**2 )
               IF ( EKEMIN .LT. 0.02D+00 ) THEN
                  SIGMAN = 0.3333333333333333D+00 * SIGMAP
               ELSE
                  SIGMAN = 10.63D+00 / BETAPR**2 - 29.92D+00 / BETAPR
     &                   + 42.9D+00
               END IF
            ELSE
               SIGMAP = 34.10D+00 / BETAPR**2 - 82.2D+00 / BETAPR
     &                + 82.2D+00
               SIGMAN = 10.63D+00 / BETAPR**2 - 29.92D+00 / BETAPR
     &                + 42.9D+00
            END IF
         ELSE
            STOP 'Sigfer: EKE'
         END IF
      GO TO 4000
  900 CONTINUE
      GO TO 4000
 1000 CONTINUE
      GO TO 4000
 1100 CONTINUE
      GO TO 4000
 1200 CONTINUE
      GO TO 4000
 1300 CONTINUE
      GO TO 4000
 1400 CONTINUE
      GO TO 4000
 1500 CONTINUE
      GO TO 4000
 1600 CONTINUE
      GO TO 4000
 1700 CONTINUE
      GO TO 4000
 1800 CONTINUE
      GO TO 4000
 1900 CONTINUE
      GO TO 4000
 2000 CONTINUE
      GO TO 4000
 2100 CONTINUE
      GO TO 4000
 2200 CONTINUE
      GO TO 4000
 2300 CONTINUE
      GO TO 4000
 2400 CONTINUE
      GO TO 4000
 2500 CONTINUE
      GO TO 4000
 2600 CONTINUE
      GO TO 4000
 2700 CONTINUE
      GO TO 4000
 2800 CONTINUE
      GO TO 4000
 2900 CONTINUE
      GO TO 4000
 3000 CONTINUE
      GO TO 4000
 3100 CONTINUE
      GO TO 4000
 3200 CONTINUE
      GO TO 4000
 3300 CONTINUE
      GO TO 4000
 3400 CONTINUE
      GO TO 4000
 3500 CONTINUE
      GO TO 4000
 3600 CONTINUE
      GO TO 4000
 3700 CONTINUE
      GO TO 4000
 3800 CONTINUE
      GO TO 4000
 3900 CONTINUE
      GO TO 4000
 4000 CONTINUE
      SIGMAP = 0.1D+00 * SIGMAP
      SIGMAN = 0.1D+00 * SIGMAN
      RETURN
*=== End of subroutine sigfer =========================================*
      END
