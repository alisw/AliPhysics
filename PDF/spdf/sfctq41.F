C-----------------------------------------------------------------------
#include "pdf/pilot.h"
      SUBROUTINE SFCTQ41(DX,DSCAL,
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
C :::::::::::: Structure functions from the CTEQ collaboration version 4
C :::::::::::: Lambda = 0.181 GeV (LO)
C
#include "pdf/impdp.inc"
      PARAMETER (ISET=3)
#include "pdf/w50511.inc"
C
      X = DX
      Q = DSCAL
C
C     generate the individual structure fcn calls
C
      UPV = Ctq4Fn(ISET, 1,X,Q)
      DNV = Ctq4Fn(ISET, 2,X,Q)
      SEAU = Ctq4Fn(ISET,-1,X,Q)
      DSEA = SEAU*X
      DUPV = (UPV-SEAU)*X
      SEAD = Ctq4Fn(ISET,-2,X,Q)
      DSEAD = SEAD*X
      DDNV = (DNV-SEAD)*X
      STR = Ctq4Fn(ISET,-3,X,Q)
      DSTR = STR*X
      CHM = Ctq4Fn(ISET,-4,X,Q)
      DCHM = CHM*X
      BOT = Ctq4Fn(ISET,-5,X,Q)
      DBOT = BOT*X
      TOP = 0.
C      IF (DSCAL.GT.TMAS) TOP = Ctq4Fn(ISET,-6,X,Q)
      DTOP = TOP*X
      GL = Ctq4Fn(ISET, 0,X,Q)
      DGL = GL*X
C
      RETURN
      END
