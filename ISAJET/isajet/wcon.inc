#ifndef CERNLIB_ISAJET_WCON_INC
#define CERNLIB_ISAJET_WCON_INC
*
*
* wcon.inc
*
      COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),
     1MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),RBRWW(12,4,2),EZ,
     2AQDP(12,4),BQDP(12,4),EZDP,WFUDGE
      SAVE /WCON/
#if defined(CERNLIB_DOUBLE)
#include "isajet/wcon2.inc"
* Ignoring t=pass
#endif
#if defined(CERNLIB_SINGLE)
#include "isajet/wcon1.inc"
* Ignoring t=pass
#endif
      INTEGER   MATCH
      REAL      SIN2W,WMASS,WGAM,AQ,BQ,COUT,WCBR,CUTOFF,CUTPOW,TBRWW,
     +          RBRWW,EZ,WFUDGE
      COMMON/WCON2/CUMWBR(25,3)
      REAL CUMWBR

#endif
