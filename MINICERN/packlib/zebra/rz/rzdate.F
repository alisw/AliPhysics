*
* $Id$
*
* $Log$
* Revision 1.3  1997/08/22 08:40:50  jamie
* fix rzdate for 2000-2016. Will FAIL in 2017
*
* Revision 1.2  1996/04/24 17:26:44  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZDATE(IWORD,IDATE,ITIME,ICASE)
*
***********************************************************************
*
*       If ICASE=1 returns IDATE and ITIME as stored in IWORD
*       If ICASE=2 pack current time and date into IWORD
*       If ICASE=3 same as ICASE=2 but take IDATE,ITIME from
*                  input parameters instead of calling DATIME
*
* Called by RZCOP1,RZLDIR,RZLLOK,RZLOCK,RZMAKE,RZMDIR,RZOUT,
*           RZSAVE,RZVOUT
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 01.12.86
*            22.08.97 JDS - fix for 2000-2016. RZDATE will FAIL in 2017
*
*
***********************************************************************
#include "zebra/q_jbyt.inc"
*
      IF(ICASE.EQ.1)THEN
         ICONT = JBYT(IWORD,9,24)
         IMINUT= MOD(ICONT,60)
         IM1   = ICONT-IMINUT
         IHOUR = MOD(IM1/60,24)
         ITIME = 100*IHOUR+IMINUT
         IM2   = IM1-60*IHOUR
         IDAY  = MOD(IM2/1440,31)
         IF(IDAY.EQ.0)IDAY=31
         IM3   = IM2-1440*IDAY
         IMONTH= MOD(IM3/44640,12)
         IF(IMONTH.EQ.0)IMONTH=12
         IYEAR = (IM3-44640*IMONTH)/535680
*
*        IYEAR = 14 is the year 2000
*
         IF(IYEAR.GE.14) THEN
            IDATE = 10000*(IYEAR-14)+100*IMONTH+IDAY
         ELSE
            IDATE = 860000+10000*IYEAR+100*IMONTH+IDAY
         ENDIF
*
      ELSE
         IF(ICASE.NE.3)CALL DATIME(IDATE,ITIME)
*
*     Handle dates > 2000 (IDATE=00MMDD)
*     Note: /SLATE/ is ignored. 000101 is assumed to be 1st January 2000
*
         IF(IDATE.GE.860000) THEN
            IDAT2 = IDATE - 860000
         ELSE
            IDAT2 = IDATE + 140000
         ENDIF

         IYEAR = IDAT2/10000
         IMONTH= (IDAT2-10000*IYEAR)/100
         IDAY  = MOD(IDAT2,100)
         IHOUR = ITIME/100
         IMINUT= MOD(ITIME,100)
         ICONT2= IDAY+31*(IMONTH+12*IYEAR)
         ICONT = IMINUT+60*(IHOUR+24*ICONT2)
         CALL SBYT(ICONT,IWORD,9,24)
      ENDIF
*
      END
