#ifndef CERNLIB_ISAJET_MGSIGS_INC
#define CERNLIB_ISAJET_MGSIGS_INC
*
*
* mgsigs.inc
*
C
C          Running totals for MadGraph cross sections
C          WTTOT8/NWTTOT  = total cross section
C          WTSUM8/NWT8    = channel cross section
C          IFUNC8, IDENT8 = MadGraph function code channel flavors
C
      INTEGER MXSIG8
      PARAMETER (MXSIG8=1000)
      COMMON /MGSIGS/WTTOT8,WTSUM8(MXSIG8),WTMAX8(MXSIG8),NSIG8,
     $NWTTOT,NWT8(MXSIG8),IFUNC8(MXSIG8),IDENT8(MXLIM+2,MXSIG8),
     $ISORT8(MXSIG8)
      REAL*8 WTTOT8,WTSUM8,WTMAX8
      INTEGER NSIG8,NWTTOT,NWT8,IFUNC8,IDENT8,ISORT8
      SAVE /MGSIGS/

#endif
