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
      SUBROUTINE GSCANK
*
************************************************************************
*             Generates Kinematics for the SCAN utility                *
*                                                                      *
*          ==>Called by : <USER>,GUKINE                                *
*             Author    R.Brun  *********                              *
************************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcscal.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcflag.inc"
#include "geant321/gconsp.inc"
*
      DIMENSION VERTEX(3),PLAB(3),BUF(1)
      DATA BUF/0./
*
*     -----------------------------------------------------------------
*
      IF(LSCAN.EQ.0)CALL GSCANI
*
      VERTEX(1)=VSCAN(1)
      VERTEX(2)=VSCAN(2)
      VERTEX(3)=VSCAN(3)
      IKSCAN = 48
      DPHI   = (PHIMAX-PHIMIN)/NPHI
      DTETA  = (TETMAX-TETMIN)/NTETA
      CALL GSVERT(VERTEX,0,0,UBUF,0,NVERT)
*
      DO 20  IPHI=IPHI1,IPHIL
         IF(LQ(LSCAN-IPHI).NE.0)GO TO 20
         PHI=DEGRAD*(PHIMIN+(IPHI-0.5)*DPHI)
         DO 10 ITETA=1,NTETA
            TETA=TETMIN+(ITETA-0.5)*DTETA
            IF(MODTET.EQ.1)THEN
               THETA=2.*ATAN(EXP(-TETA))
            ELSE IF(MODTET.EQ.2) THEN
               THETA=DEGRAD*TETA
            ELSE IF(MODTET.EQ.3) THEN
               THETA=ACOS(TETA)
            ENDIF
            PLAB(1) = SIN(THETA)*COS(PHI)
            PLAB(2) = SIN(THETA)*SIN(PHI)
            PLAB(3) = COS(THETA)
            CALL GSKINE(PLAB,IKSCAN,NVERT,BUF,0,NT)
   10    CONTINUE
         GO TO 30
   20 CONTINUE
*
*              Initialize LSTEMP
*
   30 CONTINUE
      ISCUR=NTETA+1
      IF(LSTEMP.EQ.0)THEN
         CALL MZBOOK(IXCONS,LSTEMP,LSTEMP,1,'LPHI',0,0,30000,1,0)
      ENDIF
*
*              Kinematic debug (controled by ISWIT(1))
*
      IF(IDEBUG.EQ.1.AND.ISWIT(1).EQ.1) THEN
         CALL GPRINT('VERT',0)
         CALL GPRINT('KINE',0)
      ENDIF
*
      END
