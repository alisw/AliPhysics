*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:08  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_MONITOR)
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GEAMON(ICASE,CHINFO)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine for GEANT monitoring                             *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUOUT                                *
C.    *       Author    R.Brun, F.Carminati ********                   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
C.
C.    ------------------------------------------------------------------
C.
      COMMON/CWK/IWK
      CHARACTER*5  ENTRY, EXIT, BATINT
      LOGICAL FIRST, INTRAC
      CHARACTER*(*) CHINFO
      SAVE IWK1,IM1,TIM1,FIRST,ENTRY,EXIT,BATINT
      DATA ENTRY  /'GINIT'/
      DATA EXIT   /'GLAST'/
      DATA BATINT /'BATCH'/
      DATA FIRST  /.TRUE./
*
*________________________________________________________
*
      IF(ICASE.LE.1)THEN
         IF(FIRST) THEN
            CALL TIMEL(TIM1)
            IWK1=IWK
            CALL DATIME(IDAT1,ITIM1)
            IH1=ITIM1/100
            IM1=ITIM1-100*IH1+IH1*60
            IF(INTRAC()) BATINT='INTER'
            FIRST=.FALSE.
         ENDIF
      ENDIF
      IF(ICASE.EQ.0) THEN
         ENTRY='GXINT'
         EXIT ='GXINT'
      ELSE
         IF(ICASE.EQ.1) THEN
            WRITE(CHMAIL,10000)IWK1,CHINFO,ENTRY,BATINT
         ELSEIF(ICASE.EQ.2) THEN
            IF(ENTRY.EQ.'GXINT'.AND.CHINFO.EQ.'GLAST') THEN
               EXIT='GLAST'
               GOTO 999
            ENDIF
            CALL TIMEL(TIM2)
            CPTIME=TIM1-TIM2
            CALL DATIME(IDAT2,ITIM2)
            IH2=ITIM2/100
            IM2=ITIM2-100*IH2+IH2*60
            IMD=IM2-IM1
            IF(IMD.LT.0) IMD=IMD+24*60
            IRTIME=MIN(9999,MAX(IMD,1))
            CPTIME=MIN(99999.,CPTIME)
            IF(IRTIME.GT.999)THEN
               WRITE(CHMAIL,10101)IWK1,IRTIME,CPTIME,EXIT
            ELSE
               WRITE(CHMAIL,10100)IWK1,IRTIME,CPTIME,EXIT
            ENDIF
         ENDIF
*
         CALL UMLOG('GEANTMON',CHMAIL(1:LNBLNK(CHMAIL)))
*
      ENDIF
*
10000 FORMAT('LOG321 WTYP=',I6,1X,A,2(1X,A5))
10100 FORMAT('LOGOUT WTYP=',I6,' RT=',I3,' min CP=',F9.3,' sec ',A5)
10101 FORMAT('LOGOUT WTYP=',I6,' RT', I4,' min CP=',F9.3,' sec ',A5)
*
  999 END
#endif
