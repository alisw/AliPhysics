*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:09  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE WBSJA(Z,A,NMAX,ND,CB)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE CBSJA(Z,A,NMAX,ND,CB)
#endif
#include "gen/imp64.inc"
#include "gen/defc64.inc"
     +    Z,ZA,CB,I,CBA,RR,F,C,S,R,GCONJG,ZZ
      REAL R0,D,T,Q,U,V,TC(11)
      CHARACTER*80 ERRTXT
      CHARACTER NAME*(*)
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CBSJA/WBSJA')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CBSJA')
#endif
      LOGICAL LER
      DIMENSION CB(0:*),CBA(0:100),RR(0:100)

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif

#if defined(CERNLIB_DOUBLE)
      PARAMETER (I = (0D0,1D0))
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (I = (0E0,1E0))
#endif
      PARAMETER (Z0 = 0, Z1 = 1, HF = Z1/2, Z10 = 10)

      DATA TC / 5.7941 E-5,-1.76148E-3, 2.08645E-2,-1.29013E-1,
     1          8.5777 E-1, 1.0125 E+0, 7.75   E-1, 2.3026 E+0,
     2          1.3863 E+0, 7.3576 E-1, 1.3591 E+0/

#if defined(CERNLIB_DOUBLE)
      GCONJG(ZZ)=DCONJG(ZZ)
#endif
#if !defined(CERNLIB_DOUBLE)
      GCONJG(ZZ)=CONJG(ZZ)
#endif
#include "c344cod1.inc"
#if defined(CERNLIB_DOUBLE)
      F=EXP(YA+I*(A*ATAN2(YA,X)-X))*(HF*AZ)**A/DGAMMA(1+A)
#endif
#if !defined(CERNLIB_DOUBLE)
      F=EXP(YA+I*(A*ATAN2(YA,X)-X))*(HF*AZ)**A/GAMMA(1+A)
#endif
#include "c344cod2.inc"
