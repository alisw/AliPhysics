*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:38  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_ASHO)
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
C
      SUBROUTINE GOSCMG
 
#include "geant321/gcmate.inc"
#include "geant321/gcasho.inc"
C-----------------------------------------------------------------------
C      This program calculates the parameters of the oscillators.
C      The input data are in the GCASHO. The lowest poten-
C      tials are convoluted to escape unreasonably high accuracy.
C-----------------------------------------------------------------------
C      Update history: 10-8-93
C
C-----------------------------------------------------------------------
      REAL KSI0
      PARAMETER (ENERR=0.02,TWOELM=2*511)
C-----------------------------------------------------------------------
C      Convolution of the first shells
C-----------------------------------------------------------------------
      C1 = 153.6*STEP*DENS/(BE2*AMED)
      C2 = 2.*LOG(PLIN)-BE2
      KSI0 = C1*ZMED
      X = 2.*LOG(E0MED/PLASM) + 1.
      CALL GSTCOR(X,E0MED,PLOG,DENS,D)
      A0 = LOG(TWOELM*TRNSMA/E0MED**2)+C2-D
      EMEAN = KSI0*A0
      DO 10 I = 1,NSMED
         IF (ESMED(I).GT.EMEAN*ENERR/ALFA) GOTO 20
   10 CONTINUE
   20 J = MAX(I - 1,1)
      S1 = 0.
      Z1 = 0.
      DO 30 I = 1,J
         S1 = S1 + LOG(ESMED(I))*ZSMED(I)
         Z1 = Z1 + ZSMED(I)
   30 CONTINUE
      S1  = EXP(S1/Z1)
      NOSC = NSMED - J + 1
      DO 40 I = 2,NOSC
         X = ESMED(I + J - 1)/S1
         IOSC(I) = X + 0.5
         EOSC(I) = S1*IOSC(I)*ALFA
         ZOSC(I) = ZSMED(I + J - 1)
   40 CONTINUE
      EOSC(1) = S1*ALFA
      IOSC(1) = 1
      ZOSC(1) = Z1
C-----------------------------------------------------------------------
C      Calculation of the oscillator parameters
C-----------------------------------------------------------------------
      C3 = C2 + 0.423
      DO 50 I = 1,NOSC
         X = C1*ZOSC(I)
         BOSC(I) = X/EOSC(I)
         X = 2.*LOG(EOSC(I)/PLASM)+1.
         CALL GSTCOR(X,EOSC(I),PLOG,DENS,D)
         AOSC(I) = LOG(TWOELM/EOSC(I)) +C3 -D
   50 CONTINUE
      END
#endif
