*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:36  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE GSCANU
*
************************************************************************
*             SCAN utility called at the end of each tracking step     *
*             INWVOL is different from 0 when the track has reached    *
*                    a volume boundary                                 *
*                                                                      *
*          ==>Called by : <USER>,GUSTEP                                *
*             Author    R.Brun  *********                              *
************************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcscal.inc"
#include "geant321/gcmate.inc"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcscan.inc"
      SAVE ISLOLD,IDOLD,IDNUM
*
*     -----------------------------------------------------------------
*
      IF(SLENG.LE.0.)THEN
         SX0=0.
         SABS=0.
         ITETA=ITRA
         ISLOLD=0
         IDOLD=0
      ENDIF
      SX0=SX0+STEP/RADL
      SABS=SABS+STEP/ABSL
      IF(INWVOL.NE.2)THEN
*
*            Compute new ISL,IDET
         ISL=0
         INEW=0
         DO 20 IL=NLEVEL,1,-1
            DO 10 I=1,NSLIST
               IF(NAMES(IL).EQ.ISLIST(I))THEN
                  ISL=I
                  IDNUM=NUMBER(IL)
                  IF(IL.EQ.NLEVEL.OR.IDNUM.NE.IDOLD)INEW=1
                  GO TO 30
               ENDIF
  10        CONTINUE
  20     CONTINUE
  30     CONTINUE
         IF(INWVOL.EQ.3.OR.ISL.NE.ISLOLD.OR.INEW.NE.0)THEN
            CALL GSSCAN(ISLOLD,IDOLD)
            IDOLD=IDNUM
            ISLOLD=ISL
         ENDIF
      ENDIF
*
      CALL GDEBUG
*
      END
