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
      SUBROUTINE GSCANI
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
*
*     -----------------------------------------------------------------
*
      IF(NPHI.LE.0)NPHI=1
      IF(IPHI1.LE.0)IPHI1=1
      IF(IPHI1.GT.NPHI)IPHI1=NPHI
      IF(IPHIL.LT.IPHI1)IPHIL=IPHI1
      IF(IPHIL.GT.NPHI)IPHIL=NPHI
      TETMIN = MAX(TETMID(MODTET), TETMIN)
      TETMAX = MIN(TETMAD(MODTET), TETMAX)
      IF(TETMIN.GT.TETMAX)THEN
         TETMIN = TETMID(MODTET)
         TETMAX = TETMAD(MODTET)
      ENDIF
*             Create top LSCAN bank
*
      IF(LSCAN.NE.0)THEN
         CALL MZDROP(IXCONS,LSCAN,' ')
      ENDIF
      CALL MZFORM('SCAN','10I 10F -H'      ,IOSCAN)
      NL=NPHI
      CALL MZBOOK(IXCONS,LSCAN,LSCAN,1,'SCAN',NL,NL,52,IOSCAN,0)
      IQ(LSCAN+1)=NPHI
      IQ(LSCAN+2)=NTETA
      IQ(LSCAN+3)=MODTET
      IQ(LSCAN+4)=NSLIST
      IQ(LSCAN+5)=NSLMAX
      Q (LSCAN+11)=PHIMIN
      Q (LSCAN+12)=PHIMAX
      Q (LSCAN+13)=TETMIN
      Q (LSCAN+14)=TETMAX
      Q (LSCAN+15)=VSCAN(1)
      Q (LSCAN+16)=VSCAN(2)
      Q (LSCAN+17)=VSCAN(3)
      Q (LSCAN+18)=FACTX0
      Q (LSCAN+19)=FACTL
      Q (LSCAN+20)=FACTR
      DO 10 I=1,NSLIST
         IQ(LSCAN+20+I)=ISLIST(I)
   10 CONTINUE
      END
