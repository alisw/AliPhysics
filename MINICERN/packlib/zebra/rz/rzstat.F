*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZSTAT(CHPATH,NLEVEL,CHOPT)
*
************************************************************************
*
*        To print statistics for tree CHPATH
* Input:
*   CHPATH  The pathname of the directory  about which information has to
*           be provided.
*   NLEVEL  Number of levels  below CHPATH about which  space information
*           has to be accumulated.
*   CHOPT   Character variable specifying the options desired
*           default
*                 Print the statistics
*            'Q'  Return the statistics in  the user communication vector
*                 IQUEST
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 16.05.86
*  Last mod: 24.06.92 JDS - protect against RZPAFF problems
*          : 02.05.95 JMM - set iquest(13-14)
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
      CHARACTER*(*) CHPATH,CHOPT
      DIMENSION ISD(NLPATM),NSD(NLPATM),IHDIR(4)
      DIMENSION NRUSED(NLPATM),NWUSED(NLPATM)
      DIMENSION    IOPTV(2)
      EQUIVALENCE (IOPTQ,IOPTV(1)), (IOPTX,IOPTV(2))
*
*-----------------------------------------------------------------------
*
*
      IQUEST(1)=0
      IRUSED    = 0
      IWUSED    = 0
      CALL UOPTC(CHOPT,'Q',IOPTV)
      IF(LQRS.EQ.0)GO TO 99
*
      ITIME=0
      NLEV=NLEVEL
      IF(NLEV.LE.0)NLEV=99
      CALL RZCDIR(CHWOLD,'R')
      IF(IOPTQ.EQ.0)WRITE(IQPRNT,1000)
 1000 FORMAT(//'      NREC    NWORDS    QUOTA(%)  FILE(%)   DIR. NAME')
*
*            Set CWD to the current level
*
   10 CONTINUE
      IF(ITIME.NE.0)THEN
         CALL RZPAFF(CHPAT,NLPAT,CHL)
         IF(IQUEST(1).NE.0)THEN
            NLPAT=NLPAT-1
            GO TO 20
         ENDIF
         CALL RZCDIR(CHL,' ')
      ELSE
         CALL RZCDIR(CHPATH,' ')
         IF(IQUEST(1).NE.0)GO TO 99
         NLPAT0=NLPAT
         CALL VZERO(NRUSED,NLPATM)
         CALL VZERO(NWUSED,NLPATM)
      ENDIF
      IF(IQUEST(1).NE.0)THEN
         NLPAT=NLPAT-1
         GO TO 20
      ENDIF
      ISD(NLPAT)=0
      NSD(NLPAT)=IQ(KQSP+LCDIR+KNSD)
*
*            Process possible down directories
*
   20 ISD(NLPAT)=ISD(NLPAT)+1
      IF(ISD(NLPAT).LE.NSD(NLPAT))THEN
         NLPAT=NLPAT+1
         LS=IQ(KQSP+LCDIR+KLS)
         IH=LS+7*(ISD(NLPAT-1)-1)
         CALL ZITOH(IQ(KQSP+LCDIR+IH),IHDIR,4)
         CALL UHTOC(IHDIR,4,CHPAT(NLPAT),16)
         ITIME=ITIME+1
         GO TO 10
      ELSE
*
*            Write information on current directory
*
         CALL RZPAFF(CHPAT,NLPAT,CHL)
         DO 30 I=100,1,-1
            IF(CHL(I:I).NE.' ')GO TO 40
   30    CONTINUE
   40    CONTINUE
         DO 50 K=NLPAT0,NLPAT
            NRUSED(K)=NRUSED(K)+IQ(KQSP+LCDIR+KRUSED)
            NWUSED(K)=NWUSED(K)+IQ(KQSP+LCDIR+KWUSED)
     +                         +1000000*IQ(KQSP+LCDIR+KMEGA)
   50    CONTINUE
         PERD=100.*FLOAT(NRUSED(NLPAT))/FLOAT(IQ(KQSP+LCDIR+KQUOTA))
         PERF=100.*FLOAT(NRUSED(NLPAT))/FLOAT(IQ(KQSP+LTOP +KQUOTA))
         IF(NLPAT.LE.(NLPAT0+NLEV))THEN
            IF(IOPTQ.EQ.0) WRITE(IQPRNT,2000)NRUSED(NLPAT),
     +                     NWUSED(NLPAT),PERD,PERF,CHL(1:I)
            IRUSED=NRUSED(NLPAT)
            IWUSED=NWUSED(NLPAT)
         ENDIF
 2000    FORMAT(I9,I11,2X,F8.2,1X,F8.2,3X,A)
         NRUSED(NLPAT)=0
         NWUSED(NLPAT)=0
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
   90 CALL RZCDIR(CHWOLD,' ')
*
   99 CONTINUE
      IQUEST(11)=IRUSED
      IQUEST(12)=IWUSED
*--------- Next two lines added    JMM  2 May 1995 ---------------------------
      iquest(13)=iq(kqsp+lcdir+kquota)
      iquest(14)=iq(kqsp+ltop +kquota)
*-----------------------------------------------------------------------------
      RETURN
      END
