*
* $Id$
*
* $Log$
* Revision 1.2  1996/09/30 14:25:47  ravndal
* Windows NT related modifications
*
* Revision 1.1.1.1  1995/10/24 10:22:02  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE PFNCLV.FOR
*COPY PFNCLV
*
*=== pfnclv ===========================================================*
*
      FUNCTION PFNCLV ( INC, LINNEW )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/eva0.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/nuclev.inc"
*
      LOGICAL LINNEW, LOCSEL
      DIMENSION RADSCB (2), DEFPAI (2), DEFSHE (2), NMSHLL (2)
      REAL RNDM(1)
      SAVE RADSCB, DEFPAI, DEFSHE, NMSHLL, NLEVEL, LOCSEL
      DATA LOCSEL / .FALSE. /
*
*======================================================================*
*======================================================================*
*
      I = INC
      IF ( ABS (RIMPCT) .GT. RUSNUC (I) ) THEN
         PFNCLV = - AINFNT
         RETURN
      END IF
      IF ( LINNEW ) THEN
         EKFCHC = EKFIMP
      ELSE
         EKFCHC = EKFIM2
      END IF
      EKCENT = EKFCEN (I) - EKFCHC
      EKCENT = EKFCEN (I) - EKFCHC
      EKCENT = MAX ( EKCENT, ZERZER )
      JMIN = INT ( NAVNUC (I) * ( SQRT ( EKCENT * ( EKCENT + 2.D+00
     &     * AMNUCL (I) ) ) / PFRCEN (I) )**3 ) + 1
      JMAX = NAVNUC (I)
      PROBTT = JMAX - JMIN + 1
      CALL GRNDM(RNDM,1)
      RNDJTA = RNDM (1) * PROBTT
      JMAX   = MIN ( INT (RNDJTA) + JMIN, JMAX )
      CUMMAX = CUMRAD (JMAX,I)
      IF ( LINNEW ) THEN
         NPRNUC = 1
         IF ( JUSNUC (JMAX,I) .EQ. INUCLV ) THEN
            P_FNCLV = - AINFNT
            PFNCLV = P_FNCLV
            RETURN
         ELSE IF ( JUSNUC (JMAX,I) .EQ. -INUCLV ) THEN
            CALL GRNDM(RNDM,1)
            IF ( RNDM (1) .LT. 0.5D+00 ) THEN
               P_FNCLV = - AINFNT
               PFNCLV = P_FNCLV
               RETURN
            END IF
            DEFPAI (I) = 0.D+00
         ELSE
            IF ( JUSNUC (JMAX,I) .NE. -1 ) THEN
               DEFPAI (I) = PAENUC ( NTANUC (I), I )
            ELSE
               DEFPAI (I) = 0.D+00
            END IF
         END IF
         PFCHCK = PFRIMP
         PCJMAX = PFRCEN (I) * RMASS (JMAX) / RMASS (NAVNUC(I))
         EKCJMX = SQRT ( PCJMAX**2 + AMNUSQ (I) ) - AMNUCL (I)
         EKFERM = EKFCHC + EKCJMX - EKFCEN (I)
         P_FNCLV = SQRT ( EKFERM * ( EKFERM + 2.D+00 * AMNUCL (I) ) )
         IO = 2 - I/2
         DEFPAI (IO) = 0.D+00
         DEFSHE (IO) = DEFMAG (IO)
      ELSE
         NPRNUC = NPRNUC + 1
         IF ( JUSNUC (JMAX,I) .EQ. INUCLV ) THEN
            P_FNCLV = - AINFNT
            PFNCLV = P_FNCLV
            RETURN
         ELSE IF ( I .EQ. IPRNUC (1) .AND. JMAX .EQ. JPRNUC (1) ) THEN
            IF ( JUSNUC (JMAX,I) .EQ. -INUCLV .OR. JUSNUC (JMAX,I) .EQ.
     &         -1 ) THEN
               P_FNCLV = - AINFNT
               PFNCLV = P_FNCLV
               RETURN
            END IF
         ELSE IF ( JUSNUC (JMAX,I) .EQ. -INUCLV ) THEN
            CALL GRNDM(RNDM,1)
            IF ( RNDM (1) .LT. 0.5D+00 ) THEN
               P_FNCLV = - AINFNT
               PFNCLV = P_FNCLV
               RETURN
            END IF
         ELSE
            IF ( JUSNUC (JMAX,I) .NE. -1 ) DEFPAI (I) = DEFPAI (I)
     &                                     + PAENUC ( NTANUC (I), I )
         END IF
         PFCHCK = PFRIM2
         PCJMAX = PFRCEN (I) * RMASS (JMAX) / RMASS (NAVNUC(I))
         EKCJMX = SQRT ( PCJMAX**2 + AMNUSQ (I) ) - AMNUCL (I)
         EKFER2 = EKFCHC + EKCJMX - EKFCEN (I)
         P_FNCLV = SQRT ( EKFER2 * ( EKFER2 + 2.D+00 * AMNUCL (I) ) )
      END IF
      PFNCLV = MIN ( P_FNCLV, PFCHCK )
      IF ( JMAX .EQ. NAVNUC (I) ) THEN
         RADSCB (NPRNUC) = NLSNUC (I) / ( CUMMAX - CUMRAD (JMAX-1,I) )
         NLEVEL = NLSNUC (I)
         JNUCLN = 2 * ( JMAX - 1 ) + NLEVEL
      ELSE
         RADSCB (NPRNUC) = 2.D+00 / ( CUMMAX - CUMRAD (JMAX-1,I) )
         NLEVEL = 2
         JNUCLN = 2 * JMAX
      END IF
      IPRNUC (NPRNUC) = I
      JPRNUC (NPRNUC) = JMAX
      DO 3000 MG = MAGNUC (I), 2
         IF ( MAGNUM (MG-1) .LT. JNUCLN ) GO TO 4000
 3000 CONTINUE
      MG = 1
 4000 CONTINUE
      NMSHLL (NPRNUC) = MG
      IF ( MGSNUC (MG,I) .EQ. 0 .AND. NTANUC (I) .NE. MAGNUM (MG) ) THEN
         DEFSHE (I) = SHENUC ( MAGNUM (MG) + 1, I )
     &              - SHENUC ( MAGNUM (MG), I ) + PAENUC (MAGNUM(MG),I)
     &              + DEFMAG (I)
      ELSE
         DEFSHE (I) = DEFMAG (I)
      END IF
      IF ( NUSCIN .EQ. 0 ) THEN
         DEFNUC (1) = MAX ( DEFSHE (1), ZERZER )
         DEFNUC (2) = MAX ( DEFSHE (2), ZERZER )
      ELSE
         DEFNUC (1) = MAX ( DEFPAI (1) + DEFRMI (1), DEFSHE (1), ZERZER)
         DEFNUC (2) = MAX ( DEFPAI (2) + DEFRMI (2), DEFSHE (2), ZERZER)
      END IF
      RETURN
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
      ENTRY NCLVFX
      NCLVFX = PIPIPI
      DO 5000 N=1,NPRNUC
         I = IPRNUC (N)
         J = JPRNUC (N)
         IF ( MOD (MGSNUC(NMSHLL(N),I),2) .EQ. 0 )
     &      MGSNUC (NMSHLL(N),I) = MGSNUC (NMSHLL(N),I) + 1
         IF ( JUSNUC (J,I) .EQ. -INUCLV .OR. JUSNUC (J,I) .EQ. -1 ) THEN
            JUSNUC (J,I) = INUCLV
            IF ( J .GE. JMXNUC (I) ) THEN
               JMXNUC (I) = JMXNUC (I) - NLEVEL
               RUSNUC (I) = RADSCB (N)**0.3333333333333333D+00
            END IF
         ELSE
            JUSNUC (J,I) = -INUCLV
         END IF
         NUSNUC (I) = NUSNUC (I) + 1
 5000 CONTINUE
      RETURN
*=== End of function pfnclv ===========================================*
      END
