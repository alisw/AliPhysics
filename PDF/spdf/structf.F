*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:31:55  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:06  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STRUCTF(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,TOP,GL)
C
C  *********************************************************************
C  *                                                                   *
C  *             Interface to main steering routine                    *
C  *             for all sets of structure functions                   *
C  *                                                                   *
C  *                                                                   *
C  *   Input:    X     = x value of parton                             *
C  *             SCALE = QCD scale in GeV                              *
C  *                                                                   *
C  *   Output:   UPV   = up valence quark                              *
C  *             DNV   = down valence quark                            *
C  *             SEA   = sea (up and down)                             *
C  *             STR   = strange quark                                 *
C  *             CHM   = charm quark                                   *
C  *             BOT   = bottom quark                                  *
C  *             TOP   = top quark                                     *
C  *             GL    = gluon                                         *
C  *                                                                   *
C  *                                                                   *
C  *   The variable  MODE  should be,                                  *
C  *   the variables NFL,LO,TMAS,QCDL4,QCDL5,XMIN,XMAX,Q2MIN,Q2MAX     *
C  *   could be provided by the user via a call to the                 *
C  *   subroutine PDFSET at the initialization phase, where            *
C  *                                                                   *
C  *             MODE   = number of selected structure function        *
C  *                      (Default: MODE = 45)                         *
C  *             NFL    = desired number of flavours for alpha(s)      *
C  *                      (Default: NFL = 5)                           *
C  *             LO     = order of alpha(s) calculation                *
C  *                      (Default: LO = 2)                            *
C  *             TMAS   = top quark mass in GeV/c**2  (optional)       *
C  *                      (Default: TMAS = 100.0D0)                    *
C  *             QCDL4  = QCD scale in GeV for four flavours           *
C  *             QCDL5  = QCD scale in GeV for five flavours           *
C  *                      corresponding to QCDL4                       *
C  *             XMIN   = minimal allowed x value                      *
C  *             XMAX   = maximal allowed x value                      *
C  *             Q2MIN  = minimal allowed Q**2 value                   *
C  *             Q2MAX  = maximal allowed Q**2 value                   *
C  *                                                                   *
C  *   for each set of structure fuction.                              *
C  *                                                                   *
C  *   The internal COMMON blocks                                      *
C  *                                                                   *
C  *   COMMON/W50511/ MODE,NFL,LO,NDUMMY,TMAS                          *
C  *   COMMON/W50512/ QCDL4,QCDL5                                      *
C  *   COMMON/W50513/ XMIN,XMAX,Q2MIN,Q2MAX                            *
C  *                                                                   *
C  *   are filled then by the subroutine PDFSET.                       *
C  *                                                                   *
C  *                                                                   *
C  *     Note: STRUCTF returns X * parton distribution function !      *
C  *                                                                   *
C  *                                                                   *
C  *     Author:   H. Plothow-Besch                                    *
C  *               CERN-PPE, CH - 1211 Geneva 23, Switzerland          *
C  *                                                                   *
C  *     Please return any problems, questions, suggestions            *
C  *     to the author                                                 *
C  *                                                                   *
C  *********************************************************************
C
C
#include "pdf/impdp.inc"
C
        CALL STRUCTM(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
        SEA = (USEA + DSEA)/2.
C
      RETURN
      END
