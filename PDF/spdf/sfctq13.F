*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:30  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:32  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C-----------------------------------------------------------------------
      SUBROUTINE SFCTQ13(DX,DSCAL,
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
C :::::::::::: Structure functions from Morfin and Wu-Ki Tung
C :::::::::::: Lambda = 0.231 GeV (MSb - singular gluon) (CTEQ1MS)
C
#include "pdf/expdp.inc"
     +       DX,DSCAL,DUPV,DDNV,DSEA,DSEAD,DSTR,DCHM,DBOT,DTOP,DGL
      REAL   X, Q
      PARAMETER (ISET=2)
#include "pdf/w50511.inc"
C
      X = DX
      Q = DSCAL
C
C     generate the individual structure fcn calls
C
      UPV = Ctq1Pd(ISET, 1,X,Q,IFL)
      DUPV = UPV
      DNV = Ctq1Pd(ISET, 2,X,Q,IFL)
      DDNV = DNV
      SEAU = Ctq1Pd(ISET,-1,X,Q,IFL)
      DSEA = SEAU
      SEAD = Ctq1Pd(ISET,-2,X,Q,IFL)
      DSEAD = SEAD
      STR = Ctq1Pd(ISET,-3,X,Q,IFL)
      DSTR = STR
      CHM = Ctq1Pd(ISET,-4,X,Q,IFL)
      DCHM = CHM
      BOT = Ctq1Pd(ISET,-5,X,Q,IFL)
      DBOT = BOT
      TOP = 0.
      IF (DSCAL.GT.TMAS) TOP = Ctq1Pd(ISET,-6,X,Q,IFL)
      DTOP = TOP
      GL = Ctq1Pd(ISET, 0,X,Q,IFL)
      DGL = GL
C
      RETURN
      END
