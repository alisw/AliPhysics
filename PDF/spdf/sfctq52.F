C-----------------------------------------------------------------------
#include "pdf/pilot.h"
      SUBROUTINE SFCTQ52(DX,DSCAL,
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
C :::::::::::: Structure functions from the CTEQ collaboration version 5
C :::::::::::: Lambda = 0.326 GeV (DIS)
C
#include "pdf/impdp.inc"
      PARAMETER (ISET=2)
#include "pdf/w50511.inc"
C
      X = DX
      Q = DSCAL
C
C     generate the individual structure fcn calls
C
      UPV = Ctq5Pdf(ISET, 1,X,Q)
      DNV = Ctq5Pdf(ISET, 2,X,Q)
      SEAU = Ctq5Pdf(ISET,-1,X,Q)
      DSEA = SEAU*X
      DUPV = (UPV-SEAU)*X
      SEAD = Ctq5Pdf(ISET,-2,X,Q)
      DSEAD = SEAD*X
      DDNV = (DNV-SEAD)*X
      STR = Ctq5Pdf(ISET,-3,X,Q)
      DSTR = STR*X
      CHM = Ctq5Pdf(ISET,-4,X,Q)
      DCHM = CHM*X
      BOT = Ctq5Pdf(ISET,-5,X,Q)
      DBOT = BOT*X
      TOP = 0.
C      IF (DSCAL.GT.TMAS) TOP = Ctq5Pdf(ISET,-6,X,Q)
      DTOP = TOP*X
      GL = Ctq5Pdf(ISET, 0,X,Q)
      DGL = GL*X
C
      RETURN
      END
