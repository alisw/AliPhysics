*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:25  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GFSHDC(IELEM,Z)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Fetch Shell Decay Constants                                   *
C.    *                                                                *
C.    *    ==>CALLED BY : GPHXSI                                       *
C.    *       AUTHOR    : J. Chwastowski                               *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
      REAL ONEEV
      PARAMETER (ONEEV = 1.E-9)
      DIMENSION PRB(4),PRBR(92),ER(92),PRBNR(92),ENR(92)
      DIMENSION ESHL(24),NRAD(24),NONRAD(24)
      PARAMETER (NFNBIN = 17, NSHELL = 4)
C
C Push NZ JPFN banks which will contain constants for each Z
C
      JPHXS = LQ(JPHOT-1)
C
C Get Z, the shell potentials and the decay modes
C
      DO 10 I = 1,24
         ESHL(I) = 0.0
         NRAD(I) = 0
         NONRAD(I) = 0
   10 CONTINUE
      DO 20 I = 1,4
         PRB(I) = 0.0
   20 CONTINUE
      DO 30 I = 1,92
         PRBNR(I) = 0.0
         PRBNR(I) = 0.0
         ENR(I) = 0.0
         ER(I) = 0.0
   30 CONTINUE
      CALL GFSHLS(Z,ESHL,NSHLL)
      CALL GFRDT(Z,ESHL,NSHELL,NWR,NRAD,PRBR,ER)
      CALL GFNRDT(Z,ESHL,NSHELL,NWNR,NONRAD,PRBNR,ENR)
C Calculate how many words are needed for the final state bank JPHFN
      NWORD = 0
      DO 40 J = 1,NSHELL
         IF(NRAD(J).GT.0) NWORD = NWORD+2*NRAD(J)+1
         IF(NONRAD(J).GT.0) NWORD = NWORD+2*NONRAD(J)+1
   40 CONTINUE
      NBOOK = NWORD+NFNBIN
      JPHFN = LQ(JPHXS-IELEM)
C Push bank to store final state parameters
      CALL MZPUSH(IXCONS,JPHFN,0,NBOOK,'R')
      NUSED = 5*Q(JPHFN+1)+1
      JPHFN = JPHFN+NUSED
      Q(JPHFN+1) = NSHELL
C Get probability of the shell radiative decay
      CALL GFSDPR(Z,NSHELL,PRB)
C
C Copy potentials and radiative decay probabilities
C
      DO 50 J = 1,NSHELL
         IF(ESHL(J).GT.0.0) THEN
            Q(JPHFN+1+J) = ESHL(J)*ONEEV
            Q(JPHFN+1+J+NSHELL) = PRB(J)
         ELSE
C if the shell potential is zero set it to -1
            Q(JPHFN+1+J) = -1.
            Q(JPHFN+1+J+NSHELL) = -1.
         ENDIF
   50 CONTINUE
C
C Now configurations of the final state
C
      K = 18
      KR = 1
      KNR = 1
      IF(NWORD.GT.0) THEN
         Q(JPHFN+10) = 18+NUSED
         DO 100 J = 1,NSHELL
            IF(ESHL(J).GT.0.0) THEN
               IF(NRAD(J).GT.0) THEN
                  IF(J.GT.1) Q(JPHFN+9+J) = K+NUSED
                  Q(JPHFN+K) = NRAD(J)
                  K = K+1
                  KER = KR+NRAD(J)-1
                  DO 60 L = KR,KER
                     Q(JPHFN+K) = PRBR(L)
                     K = K+1
   60             CONTINUE
                  DO 70 L = KR,KER
                     Q(JPHFN+K) = ER(L)
                     K = K+1
   70             CONTINUE
                  KR = KR+NRAD(J)
               ENDIF
               IF(NONRAD(J).GT.0) THEN
                  Q(JPHFN+13+J) = K+NUSED
                  Q(JPHFN+K) = NONRAD(J)
                  K = K+1
                  KNER = KNR+NONRAD(J)-1
                  DO 80 L = KNR,KNER
                     Q(JPHFN+K) = PRBNR(L)
                     K = K+1
   80             CONTINUE
                  DO 90 L = KNR,KNER
                     Q(JPHFN+K) = ENR(L)
                     K = K+1
   90             CONTINUE
                  KNR = KNR+NONRAD(J)
               ENDIF
            ENDIF
  100    CONTINUE
      ELSE
C You should never land here unless Z < 6
         IF(Z.GT.5.) THEN
C               CALL MZDROP(IXCONS,JPHFN,'L')
            WRITE(CHMAIL,'(A25,I3)') ' GFSHDC. JPHFN Z > 5. Z = ',Z
            CALL GMAIL(0,0)
         ENDIF
      ENDIF
      END
