*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE NUPREL.FOR
*COPY NUPREL
*
*=== nuprel ===========================================================*
*
      SUBROUTINE NUPREL(KPROJ,EKIN,PLAB,CTX,CTY,CTZ)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C
C     Last change on  16  september  1991  by  Alfredo Ferrari
C
C  HJM  24/10/88
C
C                 NUCLEON-PROTON ELASTIC SCATTERING EVENTS
C                 PARAMETRIZATION FOR COS(THETA) SAMPLING
C                                 FROM HETKFA2 / CLOTH VERSION
C
C ( Now also antiproton, antineutron - proton scattering )
C-------------------------------------------------------------------
C
#include "geant321/paprop.inc"
#include "geant321/finuc.inc"
C
      AMP=AM(KPROJ)
      AMP2=AMP**2
      AMT=AM(1)
      AMT2=AMT**2
      ELAB=EKIN + AMP
      S=AMT2 + AMP2 + 2.D0*AMT*ELAB
      WS=SQRT(S)
      BGAM=PLAB/WS
      GAM=(ELAB + AMT)/WS
      ECMP=(S + AMP2 - AMT2)*0.5D0/WS
      PCM=SQRT((ECMP + AMP)*(ECMP - AMP))
      ECMT=WS - ECMP
*  +-------------------------------------------------------------------*
*  |  Horrible patch for ap, an
      IF ( KPROJ .EQ. 2 ) THEN
         KSAMC = 8
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE IF ( KPROJ .EQ. 9 ) THEN
         KSAMC = 1
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         KSAMC = KPROJ
      END IF
*  |
*  +-------------------------------------------------------------------*
C
      CALL SAMCST(KSAMC,EKIN,CST)
C
      PCPL=PCM*CST
      EPF=GAM*ECMP + BGAM*PCPL
      PPLF=BGAM*ECMP + GAM*PCPL
      ETF=GAM*ECMT - BGAM*PCPL
      PTLF=BGAM*ECMT - GAM*PCPL
C
      PPF=SQRT((EPF+AMP)*(EPF-AMP))
      CTPR=PPLF/PPF
      STPR=SQRT((1.D0-CTPR)*(1.D0+CTPR))
      PTF=SQRT((ETF+AMT)*(ETF-AMT))
      CTTA=PTLF/PTF
      STTA=SQRT((1.D0-CTTA)*(1.D0+CTTA))
C
      CALL SFECFE(SFE,CFE)
      CALL TTRANS(CTX,CTY,CTZ,CTPR,STPR,SFE,CFE,CXR(1),CYR(1),CZR(1))
      SFE=-SFE
      CFE=-CFE
      CALL TTRANS(CTX,CTY,CTZ,CTTA,STTA,SFE,CFE,CXR(2),CYR(2),CZR(2))
C
      NP=2
      KPART(1)=KPROJ
      KPART(2)=1
      TKI(1)=EPF - AMP
      PLR(1)=PPF
      TKI(2)=ETF - AMT
      PLR(2)=PTF
      TV=0.D0
C
      RETURN
      END
