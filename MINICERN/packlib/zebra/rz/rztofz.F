*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:14  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZTOFZ(LUNFZ,CHOPT)
*
************************************************************************
*
*        Copy the CWD tree to a sequential FZ file
*        The FZ file must have been declared with FZOPEN
* Input:
*   LUNFZ   Logical unit number of the FZ sequential access file
*   CHOPT   default save only the highest cycle to LUNFZ
*           'C' save all cycles
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 14.05.86
*  Last mod: 26.06.92 JDS - protect against RZPAFF problems
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
      CHARACTER*(*) CHOPT
      DIMENSION ISD(NLPATM),NSD(NLPATM),IHDIR(4)
*
*-----------------------------------------------------------------------
*
      IQUEST(1)=0
      IQ1=0
      IF(LQRS.EQ.0)GO TO 99
*
      CALL UOPTC(CHOPT,'C',IOPTC)
      NLPAT0=NLPAT
      DO 5 I=1,NLPAT0
         CHPAT(I)=CHCDIR(I)
   5  CONTINUE
      ITIME=0
      CALL RZCDIR(CHWOLD,'R')
*
*        Garbage collection in user short range divisions
*        in primary store
*
      CALL MZGARB(21,0)
*
*            Write general header
*
      IHDIR(1)=12345
      IHDIR(2)=NLPAT0
      CALL FZOUT(LUNFZ,JQPDVS,0,1,'Z',1,2,IHDIR)
      IF(IQUEST(1).NE.0)THEN
         IQ1=IQUEST(1)
         GO TO 90
      ENDIF
*
*            Set CWD to the current level
*
  10  CONTINUE
      IF(ITIME.NE.0)THEN
         CALL RZPAFF(CHPAT,NLPAT,CHL)
      IF(IQUEST(1).NE.0)THEN
         IQ1=IQUEST(1)
         NLPAT=NLPAT-1
         GO TO 20
      ENDIF
         CALL RZCDIR(CHL,' ')
      ENDIF
      ISD(NLPAT)=0
      NSD(NLPAT)=IQ(KQSP+LCDIR+KNSD)
*
*            Write current directory
*
      CALL RZTOF1(LUNFZ,IOPTC)
      IF(IQUEST(1).NE.0)THEN
         IQ1=IQUEST(1)
         NLPAT=NLPAT-1
         GO TO 20
      ENDIF
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
*             Write final trailer
*
      NLPAT=NLPAT0
      CALL FZOUT(LUNFZ,JQPDVS,0,1,'Z',1,1,99)
      IF(IQUEST(1).NE.0)THEN
         IQ1=IQUEST(1)
         GO TO 90
      ENDIF
      LCORD=LQ(KQSP+LTOP-4)
      IF(LCORD.NE.0)THEN
         CALL MZDROP(JQPDVS,LCORD,'L')
         LCORD=0
      ENDIF
*
*            Reset CWD
*
  90  CONTINUE
      CALL RZCDIR(CHWOLD,' ')
      IF(IQ1.NE.0.AND.IQUEST(1).EQ.0)IQUEST(1)=1
*
  99  RETURN
      END
