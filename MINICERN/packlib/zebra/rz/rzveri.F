*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:19  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZVERI(CHPATH,CHOPT)
*
************************************************************************
*
*           Routine to build bit map of records used in an RZ file
* Input:
*   CHOPT   Character variable specifying the selected options.
*
*           'B' - rebuild bit map in memory
*           'C' - compare bit map in memory against file
*           'O' - check for overwriting on a word by word basis
*                 This requires a suitably dimensioned array
*                 in sequence RZBMAP. KDMAX = NRECS * LRECL / 32
*           'P' - print directories and objects pointing to overwritten
*                 records. Implies O.
*
* Called by RZFILE, RZCLOS
* Based on RZVERI program of Rene Brun
*
*  Author  : J.Shiers CN/AS/DL
*  Written : 23.03.92
*  Last mod: 05.11.92 - IQUEST(2) = number of records in use
*                       but marked as free
*
************************************************************************
*
      CHARACTER*(*) CHOPT
      CHARACTER*10  CHOPTT
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
#include "zebra/rzbmap.inc"
#include "zebra/rzover.inc"
      CHARACTER *(*) CHPATH
      DIMENSION ISD(NLPATM),NSD(NLPATM),IHDIR(4)
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"

      IQUEST(1)= 0
      IRET     = 0
      JRET     = 0
      NBAD     = 0
      NPASS    = 0
      LOGLV    = JBYT(IQ(KQSP+LTOP),15,3)-3
      NCHO     = LENOCC(CHOPT)
      IOPTB    = INDEX(CHOPT(1:NCHO),'B')
      IOPTC    = INDEX(CHOPT(1:NCHO),'C')
      IOPTO    = INDEX(CHOPT(1:NCHO),'O')
      IOPTP    = INDEX(CHOPT(1:NCHO),'P')
      IF(IOPTP.NE.0) IOPTO = 1
*
*     Option B: clear existing bitmap
*
      IF(IOPTB.NE.0) THEN
         LB=IQ(KQSP+LTOP+KLB)
         CALL VZERO(IQ(KQSP+LTOP+LB+3),IQ(KQSP+LTOP+LB))
      ENDIF
 
      NCHO   = 0
      CHOPTT = ' '
 
      IF(IOPTB.NE.0) THEN
         NCHO = NCHO + 1
         CHOPTT(NCHO:NCHO) = 'B'
      ENDIF
 
      IF(IOPTC.NE.0) THEN
         NCHO = NCHO + 1
         CHOPTT(NCHO:NCHO) = 'C'
      ENDIF
 
      IF(IOPTO.NE.0) THEN
         NCHO = NCHO + 1
         CHOPTT(NCHO:NCHO) = 'O'
      ENDIF
*
      IF(IOPTO.NE.0) THEN
         IWORD = 0
         DO 10 I=1,32
            CALL SBIT1(IWORD,I)
   10    CONTINUE
         DO 20 I=1,KDMAX
            IDATA(I) = IWORD
   20    CONTINUE
      ENDIF
      IF(LQRS.EQ.0)GOTO 70
*
*            General case
*
      IF(LCDIR.EQ.0)GOTO 70
   30 CONTINUE
      NPASS = NPASS + 1
      CALL RZCDIR(CHWOLD,'R')
      CALL RZCDIR(CHPATH,' ')
      CALL RZPAFF(CHPAT,NLPAT,CHL)
      NLPAT0=NLPAT
      ITIME=0
*
*
*            Set CWD to the current level
*
   40 CONTINUE
      IF(ITIME.NE.0)THEN
         CALL RZPAFF(CHPAT,NLPAT,CHL)
         IF(IQUEST(1).NE.0)THEN
            IF(LOGLV.GE.1) THEN
               LCHL = LENOCC(CHL)
               WRITE(IQPRNT,*) 'RZVERI. error setting directory to ',
     +         CHL(1:LCHL)
            ENDIF
            IRET = IRET + IQUEST(1)
            NLPAT=NLPAT-1
            GOTO 50
         ENDIF
         CALL RZCDIR(CHL,' ')
      ENDIF
      IF(IQUEST(1).NE.0)THEN
         IF(LOGLV.GE.1) THEN
            LCHL = LENOCC(CHL)
            WRITE(IQPRNT,*) 'RZVERI. error setting directory to ',
     +         CHL(1:LCHL)
         ENDIF
         IRET = IRET + IQUEST(1)
         NLPAT=NLPAT-1
         GOTO 50
      ENDIF
      ISD(NLPAT)=0
      NSD(NLPAT)=IQ(KQSP+LCDIR+KNSD)
*
*            Check current directory
*
      CALL RZVER1(CHL,CHOPTT,ISTAT)
      IRET = IRET + ISTAT
      JRET = JRET + IQUEST(2)
*
*            Process possible down directories
*
   50 ISD(NLPAT)=ISD(NLPAT)+1
      IF(ISD(NLPAT).LE.NSD(NLPAT))THEN
         NLPAT=NLPAT+1
         LS=IQ(KQSP+LCDIR+KLS)
         IH=LS+7*(ISD(NLPAT-1)-1)
         CALL ZITOH(IQ(KQSP+LCDIR+IH),IHDIR,4)
         CALL UHTOC(IHDIR,4,CHPAT(NLPAT),16)
         ITIME=ITIME+1
         GOTO 40
      ELSE
         NLPAT=NLPAT-1
         IF(NLPAT.GE.NLPAT0)THEN
            LUP=LQ(KQSP+LCDIR+1)
            CALL MZDROP(JQPDVS,LCDIR,' ')
            LCDIR=LUP
            GOTO 50
         ENDIF
      ENDIF
 
   60 CONTINUE
*
*     Print directories and objects using overwritten records
*
      IF(IOPTP.NE.0.AND.NPASS.EQ.1.AND.IRET.NE.0) THEN
         CHOPTT = 'P'
         WRITE(IQPRNT,*)
     +      'RZVERI. List of suspect directories/objects'
         GOTO 30
      ENDIF
*
*     Mark top directory as modified
*
      IF(IOPTB.NE.0) CALL SBIT1(IQ(KQSP+LTOP),2)
*
*     Reset CWD
*
      CALL RZCDIR(CHWOLD,' ')
*
   70 IQUEST(1) = IRET
      IQUEST(2) = JRET
      END
