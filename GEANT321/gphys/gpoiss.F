*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:32  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GPOISS(AMVEC,NPVEC,LEN)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Generates  a vector NPVEC of LEN random numbers               *
C.    *  POISSON distribued with mean values AMVEC                     *
C.    *                                                                *
C.    *  If the mean value A greater than PLIM, N is calculated        *
C.    *  according to the Gaussian approximation of the Poisson        *
C.    *  distribution.                                                 *
C.    *                                                                *
C.    *     ==> Called by : GLANDZ,GMCOUL                              *
C.    *                                                                *
C.    *  Author      : L.Urban                                         *
C.    *  Date        : 28.04.1988           Last update :  1.02.1990   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      REAL AMVEC(*),RNDM(2), N
      INTEGER NPVEC(*)
      PARAMETER (PLIM=16.,HMXINT=2E+9)
*
      DO 30 I=1,LEN
*     Protection against negative mean values
         N=0.
         IF(AMVEC(I).GT.0.) THEN
            IF(AMVEC(I).LE.PLIM) THEN
               CALL GRNDM(RNDM,1)
               R=RNDM(1)
               P=EXP(-AMVEC(I))
               S=P
               IF(R.LE.S) GOTO 20
   10          N=N+1.
               P=P*AMVEC(I)/N
               S=S+P
               IF(S.LT.R.AND.P.GT.1.E-30) GOTO 10
            ELSE
               CALL GRNDM(RNDM,2)
               RR=SQRT(-2.*LOG(RNDM(1)))
               PHI=TWOPI*RNDM(2)
               X=RR*COS(PHI)
               N=MIN(MAX(AMVEC(I)+X*SQRT(AMVEC(I)),0.),HMXINT)
            ENDIF
         ENDIF
*
   20    NPVEC(I) = N
   30 CONTINUE
*
      END
