*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE SBCOMP.FOR
*COPY SBCOMP
*
*=== sbcomp ===========================================================*
*
      SUBROUTINE SBCOMP ( SBHAL0, SBSKI0, SBCEN0, SBCEN1, SBSKI1,
     &                    SBHAL1 )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/nucgeo.inc"
*
      SAVE  R0HLA , R0HLB , R0SKA , R0SKB , R0CEA , R0CEB ,
     &      R1CEA , R1CEB , R1SKA , R1SKB , R1HLA , R1HLB ,
     &      SQRH0A, SQRH0B, SQRS0A, SQRS0B, SQRC0A, SQRC0B,
     &      SQRC1A, SQRC1B, SQRS1A, SQRS1B, SQRH1A, SQRH1B,
     &      BIMPSQ
#include "geant321/nucstf.inc"
*
      IF ( ABS (R0TRAJ) .LT. BIMPCT .OR. ABS (R1TRAJ) .LT. BIMPCT )
     &   STOP 'BIMPCT'
      BIMPSQ = BIMPCT * BIMPCT
      IF ( R0TRAJ .GE. 0.D+00 ) THEN
         SBHAL0 = 0.D+00
         SBSKI0 = 0.D+00
         SBCEN0 = 0.D+00
      ELSE
         IF ( BIMPCT .GE. RADIU1 ) THEN
            SBSKI0 = 0.D+00
            SBCEN0 = 0.D+00
            SBCEN1 = 0.D+00
            SBSKI1 = 0.D+00
            R0HLA = ABS (R0TRAJ)
            R0HLB = MAX ( BIMPCT, - R1TRAJ )
            SQRH0A = SQRT ( ( R0HLA - BIMPCT ) * ( R0HLA + BIMPCT ) )
            SQRH0B = SQRT ( ( R0HLB - BIMPCT ) * ( R0HLB + BIMPCT ) )
            SBHAL0 = RHOSKN * ( ( 1.D+00 + RADIU1 / HALODP ) * ( SQRH0A
     &             - SQRH0B ) - 0.5D+00 * ( R0HLA / HALODP * SQRH0A
     &             - R0HLB / HALODP * SQRH0B + BIMPSQ / HALODP
     &             * LOG ( ( R0HLA + SQRH0A ) / ( R0HLB + SQRH0B ) ) ) )
         ELSE IF ( R0TRAJ .GE. -RADIU0 ) THEN
            SBHAL0 = 0.D+00
            SBSKI0 = 0.D+00
            R0CEA = - R0TRAJ
            R0CEB = MAX ( BIMPCT, - R1TRAJ )
            SQRC0A = SQRT ( ( R0CEA - BIMPCT ) * ( R0CEA + BIMPCT ) )
            SQRC0B = SQRT ( ( R0CEB - BIMPCT ) * ( R0CEB + BIMPCT ) )
            SBCEN0 = RHOCEN * ( SQRC0A - SQRC0B )
         ELSE IF ( R0TRAJ .GE. -RADIU1 ) THEN
            SBHAL0 = 0.D+00
            RHELP = MAX ( BIMPCT, - R1TRAJ )
            R0SKA = ABS (R0TRAJ)
            R0SKB = MAX ( RHELP, RADIU0 )
            SQRS0A = SQRT ( ( R0SKA - BIMPCT ) * ( R0SKA + BIMPCT ) )
            SQRS0B = SQRT ( ( R0SKB - BIMPCT ) * ( R0SKB + BIMPCT ) )
            SBSKI0 = RHOCEN * ( ( 1.D+00 + RADIU0 / SKNEFF ) * ( SQRS0A
     &             - SQRS0B ) - 0.5D+00 * ( R0SKA / SKNEFF * SQRS0A
     &             - R0SKB / SKNEFF * SQRS0B + BIMPSQ / SKNEFF
     &             * LOG ( ( R0SKA + SQRS0A ) / ( R0SKB + SQRS0B ) ) ) )
            IF ( RHELP .LT. RADIU0 ) THEN
               R0CEA = RADIU0
               R0CEB = RHELP
               SQRC0A = SQRT ( ( R0CEA - BIMPCT ) * ( R0CEA + BIMPCT ) )
               SQRC0B = SQRT ( ( R0CEB - BIMPCT ) * ( R0CEB + BIMPCT ) )
               SBCEN0 = RHOCEN * ( SQRC0A - SQRC0B )
            ELSE
               SBCEN0 = 0.D+00
            END IF
         ELSE
            RHELP = MAX ( BIMPCT, - R1TRAJ )
            R0HLA = ABS (R0TRAJ)
            R0HLB = MAX ( RHELP, RADIU1 )
            SQRH0A = SQRT ( ( R0HLA - BIMPCT ) * ( R0HLA + BIMPCT ) )
            SQRH0B = SQRT ( ( R0HLB - BIMPCT ) * ( R0HLB + BIMPCT ) )
            SBHAL0 = RHOSKN * ( ( 1.D+00 + RADIU1 / HALODP ) * ( SQRH0A
     &             - SQRH0B ) - 0.5D+00 * ( R0HLA / HALODP * SQRH0A
     &             - R0HLB / HALODP * SQRH0B + BIMPSQ / HALODP
     &             * LOG ( ( R0HLA + SQRH0A ) / ( R0HLB + SQRH0B ) ) ) )
            IF ( RHELP .LT. RADIU1 ) THEN
               R0SKA = RADIU1
               R0SKB = MAX ( RHELP, RADIU0 )
               SQRS0A = SQRT ( ( R0SKA - BIMPCT ) * ( R0SKA + BIMPCT ) )
               SQRS0B = SQRT ( ( R0SKB - BIMPCT ) * ( R0SKB + BIMPCT ) )
               SBSKI0 = RHOCEN * ( ( 1.D+00 + RADIU0 / SKNEFF )
     &                * ( SQRS0A - SQRS0B ) - 0.5D+00 * ( R0SKA / SKNEFF
     &                * SQRS0A - R0SKB / SKNEFF * SQRS0B + BIMPSQ
     &                / SKNEFF * LOG ( ( R0SKA + SQRS0A ) / ( R0SKB
     &                + SQRS0B ) ) ) )
               IF ( RHELP .LT. RADIU0 ) THEN
                  R0CEA = RADIU0
                  R0CEB = MAX ( BIMPCT, - R1TRAJ )
                  SQRC0A = SQRS0B
                  SQRC0B = SQRT ( ( R0CEB - BIMPCT )*( R0CEB + BIMPCT ))
                  SBCEN0 = RHOCEN * ( SQRC0A - SQRC0B )
               ELSE
                  SBCEN0 = 0.D+00
               END IF
            ELSE
               SBSKI0 = 0.D+00
               SBCEN0 = 0.D+00
            END IF
         END IF
      END IF
      IF ( R1TRAJ .EQ. - R0TRAJ ) THEN
         R1HLA  = R0HLB
         R1HLB  = R0HLA
         SQRH1A = SQRH0B
         SQRH1B = SQRH0A
         R1SKA  = R0SKB
         R1SKB  = R0SKA
         SQRS1A = SQRS0B
         SQRS1B = SQRS0A
         R1CEA  = R0CEB
         R1CEB  = R0CEA
         SQRC1A = SQRC0B
         SQRC1B = SQRC0A
         SBCEN1 = SBCEN0
         SBSKI1 = SBSKI0
         SBHAL1 = SBHAL0
      ELSE IF ( R1TRAJ .LE. 0.D+00 ) THEN
         SBCEN1 = 0.D+00
         SBSKI1 = 0.D+00
         SBHAL1 = 0.D+00
      ELSE
         IF ( BIMPCT .GE. RADIU1 ) THEN
            SBSKI1 = 0.D+00
            SBCEN1 = 0.D+00
            R1HLB = R1TRAJ
            R1HLA = MAX ( BIMPCT, R0TRAJ )
            SQRH1A = SQRT ( ( R1HLA - BIMPCT ) * ( R1HLA + BIMPCT ) )
            SQRH1B = SQRT ( ( R1HLB - BIMPCT ) * ( R1HLB + BIMPCT ) )
            SBHAL1 = RHOSKN * ( ( 1.D+00 + RADIU1 / HALODP ) * ( SQRH1B
     &             - SQRH1A ) - 0.5D+00 * ( R1HLB / HALODP * SQRH1B
     &             - R1HLA / HALODP * SQRH1A + BIMPSQ / HALODP
     &             * LOG ( ( R1HLB + SQRH1B ) / ( R1HLA + SQRH1A ) ) ) )
         ELSE IF ( R1TRAJ .LE. RADIU0 ) THEN
            SBHAL1 = 0.D+00
            SBSKI1 = 0.D+00
            R1CEB = R1TRAJ
            R1CEA = MAX ( BIMPCT, R0TRAJ )
            SQRC1A = SQRT ( ( R1CEA - BIMPCT ) * ( R1CEA + BIMPCT ) )
            SQRC1B = SQRT ( ( R1CEB - BIMPCT ) * ( R1CEB + BIMPCT ) )
            SBCEN1 = RHOCEN * ( SQRC1B - SQRC1A )
         ELSE IF ( R1TRAJ .LE. RADIU1 ) THEN
            SBHAL1 = 0.D+00
            RHELP = MAX ( BIMPCT, R0TRAJ )
            R1SKB = R1TRAJ
            R1SKA = MAX ( RHELP, RADIU0 )
            SQRS1A = SQRT ( ( R1SKA - BIMPCT ) * ( R1SKA + BIMPCT ) )
            SQRS1B = SQRT ( ( R1SKB - BIMPCT ) * ( R1SKB + BIMPCT ) )
            SBSKI1 = RHOCEN * ( ( 1.D+00 + RADIU0 / SKNEFF ) * ( SQRS1B
     &             - SQRS1A ) - 0.5D+00 * ( R1SKB / SKNEFF * SQRS1B
     &             - R1SKA / SKNEFF * SQRS1A + BIMPSQ / SKNEFF
     &             * LOG ( ( R1SKB + SQRS1B ) / ( R1SKA + SQRS1A ) ) ) )
            IF ( RHELP .LT. RADIU0 ) THEN
               R1CEB = RADIU0
               R1CEA = RHELP
               SQRC1A = SQRT ( ( R1CEA - BIMPCT ) * ( R1CEA + BIMPCT ) )
               SQRC1B = SQRT ( ( R1CEB - BIMPCT ) * ( R1CEB + BIMPCT ) )
               SBCEN1 = RHOCEN * ( SQRC1B - SQRC1A )
            ELSE
               SBCEN1 = 0.D+00
            END IF
         ELSE
            RHELP = MAX ( BIMPCT, R0TRAJ )
            R1HLB = R1TRAJ
            R1HLA = MAX ( RHELP, RADIU1 )
            SQRH1A = SQRT ( ( R1HLA - BIMPCT ) * ( R1HLA + BIMPCT ) )
            SQRH1B = SQRT ( ( R1HLB - BIMPCT ) * ( R1HLB + BIMPCT ) )
            SBHAL1 = RHOSKN * ( ( 1.D+00 + RADIU1 / HALODP ) * ( SQRH1B
     &             - SQRH1A ) - 0.5D+00 * ( R1HLB / HALODP * SQRH1B
     &             - R1HLA / HALODP * SQRH1A + BIMPSQ / HALODP
     &             * LOG ( ( R1HLB + SQRH1B ) / ( R1HLA + SQRH1A ) ) ) )
            IF ( RHELP .LT. RADIU1 ) THEN
               R1SKB = RADIU1
               R1SKA = MAX ( RHELP, RADIU0 )
               SQRS1A = SQRT ( ( R1SKA - BIMPCT ) * ( R1SKA + BIMPCT ) )
               SQRS1B = SQRT ( ( R1SKB - BIMPCT ) * ( R1SKB + BIMPCT ) )
               SBSKI1 = RHOCEN * ( ( 1.D+00 + RADIU0 / SKNEFF )
     &                * ( SQRS1B - SQRS1A ) - 0.5D+00 * ( R1SKB / SKNEFF
     &                * SQRS1B - R1SKA / SKNEFF * SQRS1A + BIMPSQ
     &                / SKNEFF * LOG ( ( R1SKB + SQRS1B ) / ( R1SKA
     &                + SQRS1A ) ) ) )
               IF ( RHELP .LT. RADIU0 ) THEN
                  R1CEB = RADIU0
                  R1CEA = RHELP
                  SQRC1B = SQRS1A
                  SQRC1A = SQRT ( ( R1CEA - BIMPCT )*( R1CEA + BIMPCT ))
                  SBCEN1 = RHOCEN * ( SQRC1B - SQRC1A )
               ELSE
                  SBCEN1 = 0.D+00
               END IF
            ELSE
               SBCEN1 = 0.D+00
               SBSKI1 = 0.D+00
            END IF
         END IF
      END IF
      RETURN
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
      ENTRY RSCOMP ( SBHAL0, SBSKI0, SBCEN0, SBCEN1, SBSKI1, SBHAL1 )
      SBTTOT = SBHAL0 + SBSKI0 + SBCEN0 + SBCEN1 + SBSKI1 + SBHAL1
      IF ( SBUSED .GT. SBTTOT ) STOP 'RSCOMP'
      IF ( SBUSED .LT. SBHAL0 ) THEN
         SBUSEF = SBUSED
         RDLESS = R0HLA
         RDMORE = R0HLB
         SBLESS = 0.D+00
         SBMORE = SBHAL0
         DO 1000 IBIS = 1,4
            RIMPCT = 0.5D+00 * ( RDMORE + RDLESS )
            SQRIMP = SQRT ( ( RIMPCT - BIMPCT ) * ( RIMPCT + BIMPCT ) )
            SBRIMP = RHOSKN * ( ( 1.D+00 + RADIU1 / HALODP ) * ( SQRH0A
     &             - SQRIMP ) - 0.5D+00 * ( R0HLA / HALODP * SQRH0A
     &             - RIMPCT / HALODP * SQRIMP + BIMPSQ / HALODP
     &             * LOG ( ( R0HLA + SQRH0A ) / ( RIMPCT + SQRIMP ) ) ))
            IF ( SBRIMP .GE. SBUSEF ) THEN
               SBMORE = SBRIMP
               RDMORE = RIMPCT
            ELSE
               SBLESS = SBRIMP
               RDLESS = RIMPCT
            END IF
 1000    CONTINUE
         RIMPCT = RDLESS + ( SBUSEF - SBLESS ) / ( SBMORE - SBLESS )
     &          * ( RDMORE - RDLESS )
         RIMPCT = - RIMPCT
         SQRIMP = SQRT ( ( RIMPCT - BIMPCT ) * ( RIMPCT + BIMPCT ) )
      ELSE IF ( SBUSED .LT. SBHAL0 + SBSKI0 ) THEN
         SBUSEF = SBUSED - SBHAL0
         RDLESS = R0SKA
         RDMORE = R0SKB
         SBLESS = 0.D+00
         SBMORE = SBSKI0
         DO 2000 IBIS = 1,4
            RIMPCT = 0.5D+00 * ( RDMORE + RDLESS )
            SQRIMP = SQRT ( ( RIMPCT - BIMPCT ) * ( RIMPCT + BIMPCT ) )
            SBRIMP = RHOCEN * ( ( 1.D+00 + RADIU0 / SKNEFF ) * ( SQRS0A
     &             - SQRIMP ) - 0.5D+00 * ( R0SKA / SKNEFF * SQRS0A
     &             - RIMPCT / SKNEFF * SQRIMP + BIMPSQ / SKNEFF
     &             * LOG ( ( R0SKA + SQRS0A ) / ( RIMPCT + SQRIMP ) ) ))
            IF ( SBRIMP .GE. SBUSEF ) THEN
               SBMORE = SBRIMP
               RDMORE = RIMPCT
            ELSE
               SBLESS = SBRIMP
               RDLESS = RIMPCT
            END IF
 2000    CONTINUE
         RIMPCT = RDLESS + ( SBUSEF - SBLESS ) / ( SBMORE - SBLESS )
     &          * ( RDMORE - RDLESS )
         RIMPCT = - RIMPCT
         SQRIMP = SQRT ( ( RIMPCT - BIMPCT ) * ( RIMPCT + BIMPCT ) )
      ELSE IF ( SBUSED .LT. SBHAL0 + SBSKI0 + SBCEN0 ) THEN
         SBUSEF = SBUSED - SBHAL0 - SBSKI0
         SQRIMP = SQRC0A - SBUSEF / RHOCEN
         RIMPCT = SQRT ( SQRIMP**2 + BIMPSQ )
         RIMPCT = - RIMPCT
      ELSE IF ( SBUSED .LT. SBTTOT - SBSKI1 - SBHAL1 ) THEN
         SBUSEF = SBUSED - SBHAL0 - SBSKI0 - SBCEN0
         SQRIMP = SBUSEF / RHOCEN + SQRC1A
         RIMPCT = SQRT ( SQRIMP**2 + BIMPSQ )
      ELSE IF ( SBUSED .LT. SBTTOT - SBHAL1 ) THEN
         SBUSEF = SBUSED - SBHAL0 - SBSKI0 - SBCEN0 - SBCEN1
         RDLESS = R1SKA
         RDMORE = R1SKB
         SBLESS = 0.D+00
         SBMORE = SBSKI1
         DO 3000 IBIS = 1,4
            RIMPCT = 0.5D+00 * ( RDMORE + RDLESS )
            SQRIMP = SQRT ( ( RIMPCT - BIMPCT ) * ( RIMPCT + BIMPCT ) )
            SBRIMP = RHOCEN * ( ( 1.D+00 + RADIU0 / SKNEFF ) * ( SQRIMP
     &             - SQRS1A ) - 0.5D+00 * ( RIMPCT / SKNEFF * SQRIMP
     &             - R1SKA / SKNEFF * SQRS1A + BIMPSQ / SKNEFF
     &             * LOG ( ( RIMPCT + SQRIMP ) / ( R1SKA + SQRS1A ) ) ))
            IF ( SBRIMP .GE. SBUSEF ) THEN
               SBMORE = SBRIMP
               RDMORE = RIMPCT
            ELSE
               SBLESS = SBRIMP
               RDLESS = RIMPCT
            END IF
 3000    CONTINUE
         RIMPCT = RDLESS + ( SBUSEF - SBLESS ) / ( SBMORE - SBLESS )
     &          * ( RDMORE - RDLESS )
         SQRIMP = SQRT ( ( RIMPCT - BIMPCT ) * ( RIMPCT + BIMPCT ) )
      ELSE
         SBUSEF = SBUSED - SBHAL0 - SBSKI0 - SBCEN0 - SBCEN1 - SBSKI1
         RDLESS = R1HLA
         RDMORE = R1HLB
         SBLESS = 0.D+00
         SBMORE = SBHAL1
         DO 4000 IBIS = 1,4
            RIMPCT = 0.5D+00 * ( RDMORE + RDLESS )
            SQRIMP = SQRT ( ( RIMPCT - BIMPCT ) * ( RIMPCT + BIMPCT ) )
            SBRIMP = RHOSKN * ( ( 1.D+00 + RADIU1 / HALODP ) * ( SQRIMP
     &             - SQRH1A ) - 0.5D+00 * ( RIMPCT / HALODP * SQRIMP
     &             - R1HLA / HALODP * SQRH1A + BIMPSQ / HALODP
     &             * LOG ( ( RIMPCT + SQRIMP ) / ( R1HLA + SQRH1A ) ) ))
            IF ( SBRIMP .GE. SBUSEF ) THEN
               SBMORE = SBRIMP
               RDMORE = RIMPCT
            ELSE
               SBLESS = SBRIMP
               RDLESS = RIMPCT
            END IF
 4000    CONTINUE
         RIMPCT = RDLESS + ( SBUSEF - SBLESS ) / ( SBMORE - SBLESS )
     &          * ( RDMORE - RDLESS )
         SQRIMP = SQRT ( ( RIMPCT - BIMPCT ) * ( RIMPCT + BIMPCT ) )
      END IF
      IF ( RIMPCT .GT. 0.D+00 ) THEN
         XIMPCT = XBIMPC + CXIMPC * SQRIMP
         YIMPCT = YBIMPC + CYIMPC * SQRIMP
         ZIMPCT = ZBIMPC + CZIMPC * SQRIMP
      ELSE
         XIMPCT = XBIMPC - CXIMPC * SQRIMP
         YIMPCT = YBIMPC - CYIMPC * SQRIMP
         ZIMPCT = ZBIMPC - CZIMPC * SQRIMP
      END IF
      IF ( BIMPCT .GT. ANGLGB ) THEN
         TRUFAC = BIMPTR / BIMPCT
         XIMPTR = XIMPCT * TRUFAC
         YIMPTR = YIMPCT * TRUFAC
         ZIMPTR = ZIMPCT * TRUFAC
         RIMPTR = RIMPCT * TRUFAC
      ELSE
         XIMPTR = XIMPCT
         YIMPTR = YIMPCT
         ZIMPTR = ZIMPCT
         RIMPTR = RIMPCT
      END IF
      RETURN
*=== End of subroutine Sbcomp =========================================*
      END
