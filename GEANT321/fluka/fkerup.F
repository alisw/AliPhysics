*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:06  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*=== erup =============================================================*
*                                                                      *
      SUBROUTINE FKERUP (JFISS)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created  on   15 may 1990     by     Alfredo & Paola Sala        *
*                                              INFN - Milan            *
*     Last change  on   10-apr-93   by     Alfredo Ferrari, INFN-Milan *
*                                                                      *
*     Derived from the ERUP routine of EVAP-V, HERMES, KFA-Julich      *
*                                                                      *
*----------------------------------------------------------------------*
*
C*****MODIFIED TO OBTAIN APR,ZPR AFTER CAS + EVAP (8-68,T.W.A.)
#include "geant321/eva1.inc"
#include "geant321/forcn.inc"
#include "geant321/inpflg.inc"
#include "geant321/hetc5.inc"
#include "geant321/hetc7.inc"
#include "geant321/hettp.inc"
#include "geant321/higfis.inc"
*     COMMON / AZ /    LOWAZ
      LOGICAL LOPPAR
      DIMENSION  FPART (6)
C     -------------------------------------- CHECK PARAMETER
      FISINH=.FALSE.
*  +-------------------------------------------------------------------*
*  |  Check the excitation energy
      IF ( EX .LE. ZERZER ) THEN
*  |  No excitation energy:
         IF ( JFISS .LE. 0 ) THEN
            DO 201 I=1,6
               NPART(I)=0
  201       CONTINUE
            HEVSUM = ZERZER
         END IF
*        UU = ZERZER
         UU = EX
         RETURN
*  |
*  +-------------------------------------------------------------------*
*  |  Positive excitation energy:
      ELSE
*  |  Try evaporation
         M2 = NINT (APR)
         M3 = NINT (ZPR)
 8801    CONTINUE
         CALL FKDRES (M2,M3,EX,UU,EREC,LOPPAR,JFISS)
         FPARTT = ZERZER
*  |  +----------------------------------------------------------------*
*  |  |  No previous evaporation for this event
         IF ( JFISS .LE. 0 ) THEN
            DO 801 I=1,6
               FPART(I) = NPART(I)
               FPARTT   = FPARTT + FPART (I)
  801       CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Other evaporation trials already performed for this event
         ELSE
            DO 802 I=1,6
               FPART(I) = NPART(I)-NPARTF(I,JFISS-1)
               FPARTT   = FPARTT + FPART (I)
  802       CONTINUE
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |  No particle evaporated and pairing corrections accounted for
         IF ( FPARTT + FKEY .LT. ANGLGB ) THEN
            IF ( .NOT. LOPPAR ) GO TO 8802
            FKEY = ONEONE
            GO TO 8801
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 8802    CONTINUE
         FKEY = ZERZER
         ZPR  = ZPR - FPART(2) - FPART(3) - TWOTWO * ( FPART(5)
     &        + FPART(6)) - FPART(4)
         APR  = APR - FPART(1) - FPART(2) - TWOTWO * FPART(3)
     &        - THRTHR * ( FPART(4) + FPART(5) ) - FOUFOU * FPART(6)
         IF (IANG .GT. 0) THEN
         ELSE
            DO 440 K=1,6
               NP = NPART(K)
               IF ( JFISS .GT. 0 ) THEN
                  NP0 = NPARTF(K,JFISS-1) + 1
               ELSE
                  NP0 = 1
               END IF
               DO 410 J=NP0,NP
                  CALL RACO(COSEVP(1,J,K),COSEVP(2,J,K),COSEVP(3,J,K))
  410          CONTINUE
  440       CONTINUE
         END IF
      END IF
*  |
*  +-------------------------------------------------------------------*
      RETURN
*=== End of subroutine Erup ===========================================*
      END
