#ifndef CERNLIB_ISAJET_MGCOMS_INC
#define CERNLIB_ISAJET_MGCOMS_INC
*
*
* mgcoms.inc
*
C=====     Begin common blocks used by MadGraph
      REAL*8            GW, GWWA, GWWZ
      COMMON /COUP1/    GW, GWWA, GWWZ
      SAVE /COUP1/
      REAL*8            GAL(2),GAU(2),GAD(2),GWF(2)
      COMMON /COUP2A/   GAL,   GAU,   GAD,   GWF
      SAVE /COUP2A/
      REAL*8            GZN(2),GZL(2),GZU(2),GZD(2),G1(2)
      COMMON /COUP2B/   GZN,   GZL,   GZU,   GZD,   G1
      SAVE /COUP2B/
      REAL*8            GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMMON /COUP3/    GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      SAVE /COUP3/
      COMPLEX*16        GCHF(2,12)
      COMMON /COUP4/    GCHF
      SAVE /COUP4/
      REAL*8            WMASS,WWIDTH,ZMASS,ZWIDTH
      COMMON /VMASS1/   WMASS,WWIDTH,ZMASS,ZWIDTH
      SAVE /VMASS1/
      REAL*8            AMASS,AWIDTH,HMASS,HWIDTH
      COMMON /VMASS2/   AMASS,AWIDTH,HMASS,HWIDTH
      SAVE /VMASS2/
      REAL*8            FMASS(12), FWIDTH(12)
      COMMON /FERMIONS/ FMASS,     FWIDTH
      SAVE /FERMIONS/
      REAL*8            GG(2), G
      COMMON /COUPQCD/  GG,    G
      SAVE /COUPQCD/
C=====     End common blocks used by MadGraph

#endif
