#include "isajet/pilot.h"
      SUBROUTINE MGINIT
C
C          Initialize common blocks for MadGraph code in ISAJET
C          Note the QCD coupling constant is g=1.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
C
#include "isajet/itapes.inc"
#include "isajet/sstype.inc"
#include "isajet/mgcoms.inc"
C
      INTEGER I
      REAL AMGMW
      REAL*8 SW2
C
C          Fermion masses and widths
C
      FMASS(1) = AMGMW(IDE,1)
      FMASS(2) = 0D0
      FMASS(3) = AMGMW(IDUP,1)
      FMASS(4) = AMGMW(IDDN,1)
      FMASS(5) = AMGMW(IDMU,1)
      FMASS(6) = 0D0
      FMASS(7) = AMGMW(IDCH,1)
      FMASS(8) = AMGMW(IDST,1)
      FMASS(9) = AMGMW(IDTAU,1)
      FMASS(10)= 0D0
      FMASS(11)= AMGMW(IDTP,1)
      FMASS(12)= AMGMW(IDBT,1)
      DO 100 I=1,12
        FWIDTH(I)=0D0
100   CONTINUE
C
C          Boson masses and widths
C
      AMASS=0D0
      AWIDTH=0D0
      WMASS=AMGMW(IDW,1)
      WWIDTH=AMGMW(IDW,2)
      ZMASS=AMGMW(IDZ,1)
      ZWIDTH=AMGMW(IDZ,2)
      HMASS=AMGMW(IDH,1)
      HWIDTH=AMGMW(IDH,2)
      SW2=AMGMW(1,3)
C
C          Calls to Helas routines to set couplings
C
      CALL COUP1X(SW2,GW,GWWA,GWWZ)
      CALL COUP2X(SW2,GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1)
      CALL COUP3X(SW2,ZMASS,HMASS,GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH)
      DO 110 I=1,12
         CALL COUP4X(SW2,ZMASS,FMASS(I),GCHF(1,I))
110   CONTINUE
C
C          QCD couplings
C
      G = 1D0
      GG(1)=-G
      GG(2)=-G
      RETURN
      END
