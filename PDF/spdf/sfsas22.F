C-----------------------------------------------------------------------
#include "pdf/pilot.h"
      SUBROUTINE SFSAS22(DX,DQ2,DP2,IP2,
     +                     DUPV,DDNV,DSEA,DSEAD,DSTR,DCHM,DBOT,DTOP,DGL)
C
C   ********************************************************************
C   *                                                                  *
C   *        Interface to SASset of structure functions                *
C   *                                                                  *
C   *        Author:    H. Plothow-Besch (CERN-PPE)                    *
C   *                                                                  *
C   ********************************************************************
C
C :::::::::::: Structure functions from the SAS group version 2
C :::::::::::: Lambda = 0.200 GeV, Q**2 = 0.36 GeV**2 (MS_bar)
C
#include "pdf/expdp.inc"
     +          DX,DQ2,DP2,
     +          DUPV,DDNV,DSEA,DSEAD,DSTR,DCHM,DBOT,DTOP,DGL
      DIMENSION XPDFGM(-6:6)
      REAL X, Q, Q2, P2, F2GAM, XPDFGM
      PARAMETER (ISET=2)
#include "pdf/w50511.inc"
C
      X  = DX
      Q  = SQRT(DQ2)
      Q2 = DQ2
      P2 = DP2
      IP = IP2
      IF(IP.LT.0 .OR. IP.GT.7) IP=0
C
C     generate the individual structure fcn calls
C
      CALL SASGAM2(ISET,X,Q2,P2,IP,F2GAM,XPDFGM)
      UPV = XPDFGM(2)
      DUPV = UPV
      DNV = XPDFGM(1)
      DDNV = DNV
      SEAU = XPDFGM(-2)
      DSEA = SEAU
      SEAD = XPDFGM(-1)
      DSEAD = SEAD
      STR = XPDFGM(-3)
      DSTR = STR
      CHM = XPDFGM(-4)
      DCHM = CHM
      BOT = XPDFGM(-5)
      DBOT = BOT
      TOP = 0.
C      IF (DSCAL.GT.TMAS) TOP = XPDFGM(6)
      DTOP = TOP
      GL = XPDFGM(0)
      DGL = GL
C
      RETURN
      END
