*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE NIZLNW.FOR
*COPY NIZLNW
*                                                                      *
*=== nizlnw ===========================================================*
*                                                                      *
      SUBROUTINE NIZLNW ( IT, ZZZ, AAA, EKE, PO, SI, ZL )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C
#include "geant321/paprop.inc"
C
      PARAMETER ( AVOGMB = 1.0D+27 / AVOGAD )
C
      DIMENSION ITT (39)
      SAVE ITT
      DATA ITT/   1,   7,   0,   0,   0,   0,   0,   2,   8,   0,   0,
*              Klong pi+  pi-   K+   K- Lam  Alam Kshrt Sig- Sig+ Sig0
     &           10,   3,   4,   6,   5,   2,   8,  10,   2,   1,   2,
*               pi0   K0  AK0  pi0  res. res. res. res. Asi- Asi0 Asi+
     &            3,  10,   9,   3,   0,   0,   0,   0,   8,   8,   7,
*                X0  Ax0   X-  AX-  Om- Aom+
     &            2,   8,   2,   8,   2,   8 /
*
      SI = AZRZRZ
      ZL = AINFNT
C
C     CALCULATE THE NEW PARTICLE NUMBER IIT:   1=P,2=N,3=PI+,4=PI-,
C     5=K-,6=K+,7=P BAR,8=N BAR,9=K ZERO BAR,10=K ZERO
C
      IIT = ITT (IT)
      IF ( IIT .LE. 0 ) RETURN
*  Check the kinetic energy: no interaction below 50 MeV at present
      IF ( IIT .LE. 2 .AND. EKE .LT. 0.400D+00 .AND. AAA .GT. 1.5D+00 )
     &   THEN
         IF ( IIT .EQ. 1 ) THEN
            SI = XSEPRO ( EKE, ZZZ, AAA )
         ELSE
            SI = XSENEU ( EKE, ZZZ, AAA )
         END IF
         IF ( EKE .GT. 0.200D+00 ) THEN
            CALL NIZL ( IT, AAA, EKE, PO, SI2, ZL2 )
            WEIGH1 = ( 0.400D+00 - EKE ) / 0.200D+00
            SI = WEIGH1 * SI + ( 1.D+00 - WEIGH1 ) * SI2
         END IF
         IF ( SI .GT. 0.D+00 ) THEN
            ZL = AVOGMB * AAA / SI
         ELSE
            SI = 0.D+00
            ZL = AINFNT
         END IF
      ELSE
         CALL NIZL ( IT, AAA, EKE, PO, SI, ZL )
      END IF
      RETURN
      END
