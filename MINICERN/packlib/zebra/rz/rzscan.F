*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZSCAN(CHPATH,UROUT)
*
************************************************************************
*
*       Scan all directories from CHPATH and call user routine UROUT
*       for directory CHPATH and for every subdirectory.
* Input:
*   CHPATH  Character variable specifying the directory pathname.
*           '  ' means  the CWD.
*   UROUT   EXTERNAL user routine to be called for each directory
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 17.07.91
*  Last mod: 22.10.92 JDS - return immediately if cannot cd to CHPATH
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
      CHARACTER *(*) CHPATH
      EXTERNAL UROUT
      DIMENSION ISD(NLPATM),NSD(NLPATM),IHDIR(4)
*
*-----------------------------------------------------------------------
*
      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
*
*            General case
*
      IF(LCDIR.EQ.0)GO TO 99
      CALL RZCDIR(CHWOLD,'R')
      CALL RZCDIR(CHPATH,' ')
      IF(IQUEST(1).NE.0) GOTO 99
      CALL RZPAFF(CHPAT,NLPAT,CHL)
      NLPAT0=NLPAT
      ITIME=0
*
*
*            Set CWD to the current level
*
  10  CONTINUE
      IF(ITIME.NE.0)THEN
         CALL RZPAFF(CHPAT,NLPAT,CHL)
      IF(IQUEST(1).NE.0)THEN
         NLPAT=NLPAT-1
         GO TO 20
      ENDIF
         CALL RZCDIR(CHL,' ')
      ENDIF
      IF(IQUEST(1).NE.0)THEN
         NLPAT=NLPAT-1
         GO TO 20
      ENDIF
      ISD(NLPAT)=0
      NSD(NLPAT)=IQ(KQSP+LCDIR+KNSD)
*
*            List current directory
*
      CALL UROUT(CHL)
*
*            Process possible down directories
*
  20  ISD(NLPAT)=ISD(NLPAT)+1
      IF(ISD(NLPAT).LE.NSD(NLPAT))THEN
         NLPAT=NLPAT+1
         LS=IQ(KQSP+LCDIR+KLS)
         IH=LS+7*(ISD(NLPAT-1)-1)
         CALL ZITOH(IQ(KQSP+LCDIR+IH),IHDIR,4)
         CALL UHTOC(IHDIR,4,CHPAT(NLPAT),16)
         ITIME=ITIME+1
         GO TO 10
      ELSE
         NLPAT=NLPAT-1
         IF(NLPAT.GE.NLPAT0)THEN
            LUP=LQ(KQSP+LCDIR+1)
            CALL MZDROP(JQPDVS,LCDIR,' ')
            LCDIR=LUP
            GO TO 20
         ENDIF
      ENDIF
*
*            Reset CWD
*
  90  CALL RZCDIR(CHWOLD,' ')
*
  99  RETURN
      END
