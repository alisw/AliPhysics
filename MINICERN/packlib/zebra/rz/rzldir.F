*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZLDIR(CHPATH,CHOPT)
*
************************************************************************
*
*       List directory CHPATH
* Input:
*   CHPATH  Character variable specifying the directory pathname.
*           '  ' means list information for the CWD.
*   CHOPT   Character variable specifying the options:
*           default
*           'A' List also Keys entered with option 'A'
*           'T' List also subdirectory tree
*           'X' eXtended listing - show also OFFSET,REC1,REC2
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 06.04.86
*  Last mod: 18.11.92 JDS - add option X
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
      CHARACTER *(*) CHPATH,CHOPT
      DIMENSION    IOPTV(3)
      EQUIVALENCE (IOPTA,IOPTV(1)), (IOPTT,IOPTV(2)), (IOPTX,IOPTV(3))
      DIMENSION ISD(NLPATM),NSD(NLPATM),IHDIR(4)
*
*-----------------------------------------------------------------------
*
      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
*
*            List TOP directory only
*
      IF(CHPATH.EQ.'//')THEN
         WRITE(IQPRNT,7000)
         LRZ=LQRS
   2     IF(LRZ.NE.0)THEN
            LUNT   = IQ(KQSP+LRZ-5)
            IF(LUNT.EQ.0)GO TO 99
            CALL RZDATE(IQ(KQSP+LRZ+KDATEC),IDATEC,ITIMEC,1)
            CALL RZDATE(IQ(KQSP+LRZ+KDATEM),IDATEM,ITIMEM,1)
            CALL ZITOH(IQ(KQSP+LRZ+1),IHDIR,4)
            WRITE(IQPRNT,8000)IHDIR,LUNT,IDATEC,ITIMEC,IDATEM,ITIMEM
            LRZ=LQ(KQSP+LRZ)
            GO TO 2
         ENDIF
         GO TO 99
      ENDIF
*
*            General case
*
      IF(LCDIR.EQ.0)GO TO 99
      CALL UOPTC (CHOPT,'ATX',IOPTV)
      CALL RZCDIR(CHWOLD,'R')
      CALL RZCDIR(CHPATH,' ')
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
      CALL RZLDI1(IOPTA,IOPTX)
      IF(IOPTT.EQ.0)GO TO 90
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
 7000 FORMAT(///,' ************** List of TOP directories',/)
 8000 FORMAT(' //',4A4,'   LUN=',I6,
     +' Created ',I6,'/',I4,'  Modified ',I6,'/',I4)
  99  RETURN
      END
