*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZWRT(V,N,IOPTB,IEVENT)
*
************************************************************************
*
*      To write array V of length N into the CWD
*      RZ internal routine called by RZVOUT,RZOHN,RZOBKN
*
* Called by RZOBKN,RZOHN,RZVOUT
*
*  Author  : R.Brun, B.Holl
*  Written : 26.04.86
*  Last mod: 14.04.94 JDS. Correct(?) handling of end of records
*            29.07.94 JDS. In case of termination on EOR, last
*                          fragment of buffer was converted/copied twice
*
************************************************************************
#if !defined(CERNLIB_FQXISN)
#include "zebra/zunit.inc"
#endif
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzcout.inc"
#include "zebra/rzk.inc"
#include "zebra/mzioc.inc"
      DIMENSION V(N)
*
*-----------------------------------------------------------------------
*
*
*     Event continued?
*
      IF(IEVENT.EQ.1.AND.IP1.GT.LREC) THEN
         IF(IRLOUT.EQ.IR1)THEN
            IR=IR2
         ELSE
            IR=IRLOUT+1
         ENDIF
         IP1    = 1
         IROUT  = IR
         IRLOUT = IR
      ENDIF
 
      NREC = N + IP1 - 1
      NLEFT=LREC-IP1+1
      IF(N.LE.NLEFT)THEN
         NP1=N
      ELSE
         NP1=NLEFT
      ENDIF
#if !defined(CERNLIB_FQXISN)
      IF(IMODEX.GT.0.AND.IOPTB.NE.1)THEN
         NWFOTT = N
         NWFODN = 0
         NWFOAV = NP1
         CALL FZOCV(V,IQ(KQSP+LROUT+IP1))
         IF(NWFOAV.GT.0.OR.IFOCON(1).LT.0)GO TO 95
         IQUEST(1)=0
      ELSE
#endif
         CALL UCOPY(V,IQ(KQSP+LROUT+IP1),NP1)
#if !defined(CERNLIB_FQXISN)
      ENDIF
#endif
      IF(IP1.EQ.1)THEN
         IROUT =IRLOUT
      ENDIF
      IP1=IP1+NP1
      IF(IP1.GT.LREC)THEN
         CALL RZIODO(LUN,LREC,IRLOUT,IQ(KQSP+LROUT+1),2)
         IF(IQUEST(1).NE.0)GO TO 99
         IF(LRIN.NE.0)THEN
            IF(IQ(KQSP+LTOP+KIRIN).EQ.IRLOUT)THEN
               IQ(KQSP+LTOP+KIRIN)=0
            ENDIF
         ENDIF
         IF(IRLOUT.EQ.IR1)THEN
            IR=IR2
         ELSE
            IR=IRLOUT+1
         ENDIF
* JDS 14/04/94
         IF(IR.EQ.0) GOTO 99
   10    IP1=1
         NW=N-NP1
         IF(NW.GT.LREC)NW=LREC
         IF(NW.LT.LREC.AND.IOPTRR.NE.0)THEN
            CALL RZIODO(LUN,LREC,IR,IQ(KQSP+LROUT+1),1)
            IF(IQUEST(1).NE.0)GO TO 99
         ENDIF
#if !defined(CERNLIB_FQXISN)
         IF(IMODEX.GT.0.AND.IOPTB.NE.1)THEN
 
            IF(NWFOAV.LT.0)THEN
*
*      Case when output buffer overflowed (double precision conversion)
*
               NWFOAV=NW-1
               IQ(KQSP+LROUT+1)=IQ(KQSP+LROUT+LREC+1)
               CALL FZOCV(V,IQ(KQSP+LROUT+2))
            ELSE
*
*      Normal case
*
               NWFOAV=NW
               CALL FZOCV(V,IQ(KQSP+LROUT+1))
            ENDIF
 
            IF(NWFOAV.GT.0.OR.IFOCON(1).LT.0)GO TO 95
            IQUEST(1)=0
            IF(NW.EQ.LREC) THEN
               CALL RZIODO(LUN,NW,IR,IQ(KQSP+LROUT+1),2)
               IF(IQUEST(1).NE.0)GO TO 99
            ENDIF
         ELSE
#endif
            IF(NW.EQ.LREC) THEN
               CALL RZIODO(LUN,NW,IR,V(NP1+1),2)
               IF(IQUEST(1).NE.0)GO TO 99
            ELSE
               CALL UCOPY(V(NP1+1),IQ(KQSP+LROUT+IP1),NW)
            ENDIF
#if !defined(CERNLIB_FQXISN)
         ENDIF
#endif
         IRLOUT=IR
         IROUT =IR
         IP1=IP1+NW
         NP1=NP1+NW
*
* JDS 29/07/94
*
*        IF(NP1.GE.N) GOTO 99
*
* JDS 31/10/94
*
*     All data written - fill output buffer if on a record boundary
*
         IF(NP1.GE.N) THEN
            IF(MOD(NREC,LREC).EQ.0)THEN
#if !defined(CERNLIB_FQXISN)
               IF(IMODEX.GT.0.AND.IOPTB.NE.1)THEN
*                 NWFODN = 0
*                 IF(NWFOAV.LT.0)THEN
*                    NWFOAV=NW-1-NP1
*                    IQ(KQSP+LROUT+IP1-NW)=IQ(KQSP+LROUT+LREC+1)
*                    CALL FZOCV(V,IQ(KQSP+LROUT+IP1-NW+1))
*                 ELSE
*                    NWFOAV=NW-NP1
*                    CALL FZOCV(V,IQ(KQSP+LROUT+IP1-NW))
*                 ENDIF
*                 IF(NWFOAV.GT.0.OR.IFOCON(1).LT.0)GO TO 95
*                 IQUEST(1)=0
               ELSE
#endif
*
*     Fill output buffer for the case when vector ends on
*     a record boundary (already filled if data conversion
*     performed).
*
                  CALL UCOPY(V(NP1-NW+1),IQ(KQSP+LROUT+IP1-NW),NW)
#if !defined(CERNLIB_FQXISN)
               ENDIF
#endif
            ENDIF
            GO TO 99
         ENDIF
         IR =IR+1
         GO TO 10
      ENDIF
#if !defined(CERNLIB_FQXISN)
      GO TO 99
*
*           Conversion problem
*
   95 IQUEST(1) =5
      IQUEST(11)=NWFOTT
      IQUEST(12)=NWFORE
      IQUEST(13)=NWFOAV
      IQUEST(14)=NWFODN
      IF(JBYT(IQ(KQSP+LTOP),15,3)-3.GE.-2) WRITE(IQLOG,1000)
 1000 FORMAT(' RZOUT/RZVOUT. Error during conversion into',
     +       ' exchange format detected by RZWRT')
*1000 FORMAT(' RZWRT. Error during conversion into exchange format')
#endif
*
   99 RETURN
      END
