#include "isajet/pilot.h"
        REAL FUNCTION SSGX11(ET)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> ziss + tp + tb
C          Baer's XT11
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL MG,MT,MZ,MST1,MST2,ET
        DOUBLE PRECISION DET,DMG,DMT,DMZ,DMT1,DMT2,TOP,BOT,DXT11
        DOUBLE PRECISION XT,MUT,MUZ,XMIN,XMAX,EMIN,EMAX,SSDLAM
        DOUBLE PRECISION PI
        DATA PI/3.14159265D0/
        MG=TMP(1)
        MT=TMP(2)
        MZ=TMP(3)
        MST1=TMP(4)
        MST2=TMP(5)
        DET=ET
        DMG=TMP(1)
        DMT=TMP(2)
        DMZ=TMP(3)
        DMT1=TMP(4)
        DMT2=TMP(5)
        XT=2*ET/MG
        MUT=(MT/MG)**2
        MUZ=(MZ/MG)**2
        XMIN=((2.D0-XT)*(1.D0+2*MUT-MUZ-XT)-DSQRT(DMAX1(0.D0,
     $   (XT**2-4*MUT)*SSDLAM((1.D0+MUT-XT),MUT,MUZ))))
     $   /2.D0/(1.D0-XT+MUT)
        XMAX=((2.D0-XT)*(1.D0+2*MUT-MUZ-XT)+DSQRT(DMAX1(0.D0,
     $   (XT**2-4*MUT)*SSDLAM((1.D0+MUT-XT),MUT,MUZ))))
     $   /2.D0/(1.D0-XT+MUT)
        EMIN=XMIN*MG/2.D0
        EMAX=XMAX*MG/2.D0
        TOP=DMG**2-2*DMG*EMAX+DMT**2-DMT2**2
        BOT=DMG**2-2*DMG*EMIN+DMT**2-DMT2**2
        DXT11=-PI**2*DET*DLOG(TOP/BOT)/2.D0/(DMG**2-2*DMG*DET+DMT**2
     $         -DMT1**2)
        SSGX11=DXT11
        RETURN
        END
