#include "isajet/pilot.h"
        REAL FUNCTION SSGX8(ET)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> ziss + tp + tb
C          Baer's XT8 - Y- eq. a.6.h of prd45,142 (1992)
C          Modified for t_1 and t_2 eigenstates
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstmp.inc"
        REAL ET
        DOUBLE PRECISION DET,DMG,DMT,DMZ,DMT1,DMT2,TOP,BOT,DXT8
        DOUBLE PRECISION XT,MUT,MUZ,XMIN,XMAX,EMIN,EMAX,SSDLAM
        DOUBLE PRECISION PI
        DATA PI/3.14159265D0/
        DET=ET
        DMG=TMP(1)
        DMT=TMP(2)
        DMZ=TMP(3)
        DMT1=TMP(4)
        DMT2=TMP(5)
        XT=2*DET/DMG
        MUT=(DMT/DMG)**2
        MUZ=(DMZ/DMG)**2
        XMIN=((2.D0-XT)*(1.D0+2*MUT-MUZ-XT)-DSQRT(DMAX1(0.D0,
     $   (XT**2-4*MUT)*SSDLAM((1.D0+MUT-XT),MUT,MUZ))))
     $   /2.D0/(1.D0-XT+MUT)
        XMAX=((2.D0-XT)*(1.D0+2*MUT-MUZ-XT)+DSQRT(DMAX1(0.D0,
     $   (XT**2-4*MUT)*SSDLAM((1.D0+MUT-XT),MUT,MUZ))))
     $   /2.D0/(1.D0-XT+MUT)
        EMIN=XMIN*DMG/2.D0
        EMAX=XMAX*DMG/2.D0
        TOP=DMG**2-2*DMG*EMAX+DMT**2-DMT2**2
        BOT=DMG**2-2*DMG*EMIN+DMT**2-DMT2**2
        DXT8=((EMAX-EMIN)*(DMG**2-2*DMG*DET+DMT**2)+
     $        (DMZ**2*DMG**2-DMT2**2*DMG**2+2*DMT2**2*DMG*DET+DMT**4
     $         -DMT2**2*DMT**2)/2.D0/DMG*DLOG(TOP/BOT))*PI**2/2.D0/
     $         (DMG**2+DMT**2-2*DMG*DET-DMT1**2)
        SSGX8=DXT8
        RETURN
        END
