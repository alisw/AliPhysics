*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE NCLVIN.FOR
*COPY NCLVIN
*
*=== nclvin ===========================================================*
*
      SUBROUTINE NCLVIN
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/nuclev.inc"
#include "geant321/nucgeo.inc"
*
      INUCLV = INUCLV + 1
      IPRNUC (1) = 0
      IPRNUC (2) = 0
      JPRNUC (1) = 0
      JPRNUC (2) = 0
      NPRNUC = 0
      RUSPRO = RADTOT
      RUSNEU = RADTOT
      NUSPRO = 0
      NUSNEU = 0
      NMXPRO = NAVPRO
      NMXNEU = NAVNEU
      DO 2000 I = 1,2
         IF ( NLSNUC (I) .LE. 1 ) JUSNUC (NAVNUC(I),I) = -1
         JMXNUC (I) = NTANUC (I)
         DEFRMI (I) = 0.6666666666666666D+00 * EKFCEN (I) / NTANUC (I)
         DEFMAG (I) = 0.D+00
         DO 1000 K = 1, MAGNUC (I)
            NUNFLL = MAGNUM (K) - NTANUC (I)
            IF ( NUNFLL .GT. 0 ) THEN
               MGSNUC (K,I) = -2 + MOD ( NUNFLL, 1 )
            ELSE
               MGSNUC (K,I) = 0
               IF ( NUNFLL .EQ. 0 ) DEFMAG (I) = MAX ( ZERZER,
     &              SHENUC (MAGNUM(K)+1,I) - SHENUC (MAGNUM(K),I)
     &            + PAENUC (NTANUC(I),I) )
            END IF
 1000    CONTINUE
 2000 CONTINUE
      RETURN
*=== End of subroutine nclvin =========================================*
      END
