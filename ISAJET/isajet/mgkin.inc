#ifndef CERNLIB_ISAJET_MGKIN_INC
#define CERNLIB_ISAJET_MGKIN_INC
*
*
* mgkin.inc
*
C          Double precision PJETS; MXJETS defined in /JETLIM/
C          Format matches MadGraph
      COMMON/MGKIN/PJETS8(0:3,MXLIM+2),AMJET8(MXLIM+2)
      REAL*8 PJETS8,AMJET8
      SAVE /MGKIN/

#endif
