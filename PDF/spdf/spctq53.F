#include "pdf/pilot.h"
C-----------------------------------------------------------------------
      SUBROUTINE SPCTQ53(DX,DSCAL,
     +                     DUPV,DDNV,DSEA,DSEAD,DSTR,DCHM,DBOT,DTOP,DGL)
C
C   ********************************************************************
C   *                                                                  *
C   *        Interface to CTEQset of structure functions               *
C   *                                                                  *
C   *        Author:    H. Plothow-Besch (CERN-PPE)                    *
C   *                                                                  *
C   ********************************************************************
C
C :::::::::::: Structure functions from CETQ5L Parametrization 
C :::::::::::: Lambda = 0.192 GeV (LO) (CTEQ5L)
C
#include "pdf/impdp.inc"
C
      PARAMETER (ISET=3)
#include "pdf/w50511.inc"
C
      X = DX
      Q = DSCAL
C
C     generate the individual structure fcn calls
C
      UPV = Ctq5Lp(1,X,Q)
      DNV = Ctq5Lp(2,X,Q)
      SEAU = Ctq5Lp(-1,X,Q)
      DSEA = SEAU * X
      DUPV = (UPV-SEAU) * X
      SEAD = Ctq5Lp(-2,X,Q)
      DSEAD = SEAD * X
      DDNV = (DNV-SEAD) * X
      STR = Ctq5Lp(-3,X,Q)
      DSTR = STR * X
      CHM = Ctq5Lp(-4,X,Q)
      DCHM = CHM * X
      BOT = Ctq5Lp(-5,X,Q)
      DBOT = BOT * X
      TOP = 0.
      IF (DSCAL.GT.TMAS) TOP = Ctq5Lp(-6,X,QL)
      DTOP = TOP * X
      GL = Ctq5Lp(0,X,Q)
      DGL = GL * X
C
      RETURN
      END
