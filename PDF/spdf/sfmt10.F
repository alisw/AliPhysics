*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:30:29  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:45  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
      SUBROUTINE SFMT10(DX,DSCAL,
     +                     DUPV,DDNV,DSEA,DSEAD,DSTR,DCHM,DBOT,DTOP,DGL)
C
C   ********************************************************************
C   *                                                                  *
C   *        Interface to M-T set of structure functions               *
C   *                                                                  *
C   *        Author:    H. Plothow-Besch (CERN-PPE)                    *
C   *                                                                  *
C   ********************************************************************
C
C :::::::::::: Structure functions from Morfin and Wu-Ki Tung
C :::::::::::: Lambda = 0.237 GeV
C
#include "pdf/expdp.inc"
     +       DX,DSCAL,DUPV,DDNV,DSEA,DSEAD,DSTR,DCHM,DBOT,DTOP,DGL
      REAL   X, Q
      PARAMETER (ISET=10)
#include "pdf/w50511.inc"
C
      X = DX
      Q = DSCAL
C
C     generate the individual structure fcn calls
C
      UPV = PDZXMT(ISET, 1,X,Q,IFL)
      DUPV = UPV
      DNV = PDZXMT(ISET, 2,X,Q,IFL)
      DDNV = DNV
      SEAU = PDZXMT(ISET,-1,X,Q,IFL)
      DSEA = SEAU
      SEAD = PDZXMT(ISET,-2,X,Q,IFL)
      DSEAD = SEAD
      STR = PDZXMT(ISET,-3,X,Q,IFL)
      DSTR = STR
      CHM = PDZXMT(ISET,-4,X,Q,IFL)
      DCHM = CHM
      BOT = PDZXMT(ISET,-5,X,Q,IFL)
      DBOT = BOT
      TOP = 0.
      IF (DSCAL.GT.TMAS) TOP = PDZXMT(ISET,-6,X,Q,IFL)
      DTOP = TOP
      GL = PDZXMT(ISET, 0,X,Q,IFL)
      DGL = GL
C
      RETURN
      END
