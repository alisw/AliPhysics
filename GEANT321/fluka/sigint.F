*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_OLDNAME)
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*=== sigint ===========================================================*
*
      SUBROUTINE SIGINT (IRE ,PLAB,N,IE ,AMT ,AMN,ECM ,SI ,ITAR)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/reac.inc"
#include "geant321/redver.inc"
*
**** *,IEII,IKII,NURE
      COMMON / FKABLT / AM   (110), GA   (110), TAU  (110), ICH   (110),
     &                  IBAR (110), K1   (110), K2   (110)
*
      IE=IEFUN(PLAB,IRE)
      IF (IE.LE.IEII(IRE)) IE=IE+1
      AMT=AM(ITAR)
      AMN=AM(N)
      AMN2=AMN*AMN
      AMT2=AMT*AMT
      ECM=SQRT(AMN2+AMT2+2.D0*AMT*SQRT(AMN2+PLAB**2))
C*** INTERPOLATION PREPARATION
      ECMO=UMO(IE)
      ECM1=UMO(IE-1)
      DECM=ECMO-ECM1
      DEC=ECMO-ECM
      IIKI=IKII(IRE)+1
      EKLIM=-THRESH(IIKI)
      WOK=SIIN(IE)
      WDK=WOK-SIIN(IE-1)
      IF (ECM.GT.ECMO) WDK=0.D0
C*** INTERPOLATION IN CHANNEL WEIGHTS
      IELIM=IEFUN(EKLIM,IRE)
      DELIM=UMO(IELIM)+EKLIM
     *+1.D-16
      DETE=(ECM-(ECMO-EKLIM)*.5D0)*2.D0
      IF (DELIM*DELIM-DETE*DETE) 112,112,113
  113 DECC=DELIM
      GO TO 114
  112 DECC=DECM
  114 CONTINUE
      WKK=WOK-WDK*DEC/(DECC+1.D-9)
      IF (WKK.LT.0.D0) WKK=0.D0
      SI=WKK+1.D-12
      IF (-EKLIM.GT.ECM) SI=1.D-14
      RETURN
      END
#endif
