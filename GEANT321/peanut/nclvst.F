*
* $Id$
*
* $Log$
* Revision 1.4  1996/11/13 13:08:05  ravndal
* Paramter ZERO at least real
*
* Revision 1.3  1996/05/03 09:46:12  cernlib
* Replace 0.0 in MAX() by ZERO to get the proper type;
* add parameter statement for ZERO
*
* Revision 1.2  1996/04/26 12:20:30  ravndal
* neg. SQRT protection
*
* Revision 1.1.1.1  1995/10/24 10:22:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE NCLVST.FOR
*COPY NCLVST
*
*=== nclvst ===========================================================*
*
      SUBROUTINE NCLVST ( NA, NZ )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
       PARAMETER (ZERO=0.)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/nuclev.inc"
*
      NTAPRO = NZ
      NTANEU = NA - NZ
      NAVPRO = MOD (NTAPRO,2)
      NLSPRO = 2 - NAVPRO
      NAVPRO = NTAPRO / 2 + NAVPRO
      NAVNEU = MOD (NTANEU,2)
      NLSNEU = 2 - NAVNEU
      NAVNEU = NTANEU / 2 + NAVNEU
      PROFAC = 2.D+00
*
      CUMRAD (0,1) = 0.D+00
      CUMRAD (0,2) = 0.D+00
      DO 2000 INC = 1,2
         DDNAV = NAVNUC(INC)
         DO 1000 JNC = 1, NAVNUC (INC)
            DDJNC = JNC
            PCEJNC = PFRCEN (INC) * ( DDJNC / DDNAV
     &             )**0.3333333333333333D+00
            EKCJNC = SQRT ( PCEJNC**2 + AMNUSQ (INC) ) - AMNUCL (INC)
            EKFJNC = EKFCEN (INC) - EKCJNC
            PFRJNC = SQRT ( 
     &      MAX( EKFJNC * ( EKFJNC + 2.D+00 * AMNUCL (INC) ), ZERO) )
            RHOJNC = RHOCEN * ( PFRJNC / PFRCEN (INC) )**3
            RADJNC = FRADNC (RHOJNC)
            IF ( JNC .LT. NAVNUC (INC) ) THEN
               CUMRAD (JNC,INC) = CUMRAD (JNC-1,INC) + PROFAC
     &                          / RADJNC**3
            ELSE
               CUMRAD (JNC,INC) = CUMRAD (JNC-1,INC) + NLSNUC (INC)
     &                          / RADJNC**3
            END IF
            IF ( RADJNC .LE. RADIU0 ) THEN
               NCONUC (INC) = JNC
            ELSE IF ( RADJNC .LE. RADIU1 ) THEN
               NSKNUC (INC) = JNC
            ELSE
               NHANUC (INC) = JNC
            END IF
 1000    CONTINUE
         DO 1500 K = 1,8
            IF ( NTANUC (INC) .LE. MAGNUM (K) ) THEN
               MAGNUC (INC) = K
               GO TO 2000
            END IF
 1500    CONTINUE
 2000 CONTINUE
      RETURN
*=== End of subroutine nclvst =========================================*
      END
