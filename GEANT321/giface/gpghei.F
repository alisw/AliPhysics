*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:14  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.37  by  S.Giani
*-- Author :
      SUBROUTINE GPGHEI
C
C *** COMPUTE DISTANCE TO NEXT HADRONIC INTERACTION POINT ***
C *** THIS ROUTINE IS AN INTERFACE TO GHEISHA8 ***
C *** NVE 06-APR-1988 CERN GENEVA ***
C
C CALLED BY : GUPHAD (USER ROUTINE)
C ORIGIN : F.CARMINATI
C
#include "geant321/gcflag.inc"
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcjloc.inc"
C --- GHEISHA COMMONS ---
#include "geant321/s_prntfl.inc"
C
C --- INITIALIZE RELEVANT GHEISHA VARIABLES AT FIRST PASS ---
      IF (IFINIT(4) .EQ. 0) CALL GHEINI
C
      IF (Z .LT. 1.0) GO TO 1000
C
      KK=ABS(Q(JMA+11))
      IF (KK .GT. 1) GO TO 10
C
      SIG=GHESIG(VECT(7),GEKIN,A,A,Z,1.0,1,DENS,0.0,IPART)
      GO TO 20
C
 10   CONTINUE
      QCOR=0.0
      IF(JTM.GT.0) THEN
         LNVE=LQ(JTM)
         IF (LNVE .GT. 0) QCOR=Q(LNVE+26)
      ENDIF
      SIG=GHESIG(VECT(7),GEKIN,A,Q(JMIXT+1),Q(JMIXT+KK+1),
     $           Q(JMIXT+2*KK+1),KK,DENS,QCOR,IPART)
C
 20   CONTINUE
      IF (SIG .LE. 0.0) GO TO 1000
      SHADR=ZINTHA/SIG
      IF (NPRT(9)) PRINT 2000,KK,SIG,SHADR
 2000 FORMAT(' *GPGHEI* KK,SIG,SHADR = ',I3,1X,2(G12.5,1X))
      GO TO 9999
C
C --- ENSURE NO INTERACTION IN CURRENT MEDIUM ---
C
 1000 CONTINUE
      SHADR=BIG
      IF (NPRT(9)) PRINT 2001,KK,SIG,SHADR
 2001 FORMAT(' *GPGHEI* === NO INTERACTION IN CURRENT MEDIUM ==='/
     $ ' KK,SIG,SHADR = ',I3,1X,2(G12.5,1X))
C
 9999 CONTINUE
      END
