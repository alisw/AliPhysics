*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCLOS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Close last view bank opened and reset IVIEW to 0,        *
C.    *       meaning screen; a release of unused space in             *
C.    *       view banks is also done by negative pushing              *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
C             If current view is screen don't do anything
C
      IF (IDVIEW.EQ.0) GO TO 999
C
C             Close view number IDVIEW
C
      JV=LQ(JDRAW-IDVIEW)
C
C             Store parameters for that view, i.e. logical (I..)
C             and phisical (MAX..) length of view banks
C
      Q(JV+1)=IGU
C---  Q(JV+2)=MAXGU
      Q(JV+4)=IGS
C---  Q(JV+5)=MAXGS
      Q(JV+7)=ITU
C---  Q(JV+8)=MAXTU
      Q(JV+10)=ITS
C---  Q(JV+11)=MAXTS
C
      Q(JV+22)=ICUT
C
C             Release unused space
C
      ND1=MAXGU-IGU
      ND2=ND1
      ND3=MAXTU-ITU
      ND4=MAXGS-IGS
      ND5=ND4
      ND6=MAXTS-ITS
      JV = LQ(JV-1)
      CALL MZPUSH(IXCONS,JV,0,-ND1,'I')
      JV = LQ(JDRAW-IDVIEW)
      JV = LQ(JV-2)
      CALL MZPUSH(IXCONS,JV,0,-ND2,'I')
      JV = LQ(JDRAW-IDVIEW)
      JV = LQ(JV-3)
      CALL MZPUSH(IXCONS,JV,0,-ND3,'I')
      JV = LQ(JDRAW-IDVIEW)
      JV = LQ(JV-4)
      CALL MZPUSH(IXCONS,JV,0,-ND4,'I')
      JV = LQ(JDRAW-IDVIEW)
      JV = LQ(JV-5)
      CALL MZPUSH(IXCONS,JV,0,-ND5,'I')
      JV = LQ(JDRAW-IDVIEW)
      JV = LQ(JV-6)
      CALL MZPUSH(IXCONS,JV,0,-ND6,'I')
C
C             Reset view number to 0 (screen)
C
      IDVIEW=0
C
  999 RETURN
      END
