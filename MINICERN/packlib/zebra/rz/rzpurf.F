*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZPURF(N,KPURG)
*
************************************************************************
*
*         Update list of purged records
* Input:
*    N      Number of consecutive records purged
*    KPURG  First of the N records purged
*
* Called by RZPURG
*
*  Author  : R.Brun DD/US/PD
*  Written : 06.04.86
*  Last mod: 22.08.90
*
************************************************************************
*
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
*
*-----------------------------------------------------------------------
*
      IF(LQ(KQSP+LTOP-5).EQ.0)THEN
         CALL MZBOOK(JQPDVS,LPURG,LTOP,-5,'RZPU',0,0,21,2,0)
         IQ(KQSP+LPURG-5)=LUN
      ENDIF
*
      IF(KPURG.EQ.0)GO TO 99
      IPURG=KPURG
      IPURL=IPURG+N-1
      DO 2 I=IPURG,IPURL
         IF(I.EQ.IQ(KQSP+LCDIR+KRLOUT))THEN
            IQ(KQSP+LCDIR+KRLOUT)=0
            IQ(KQSP+LCDIR+KIP1)=1
         ENDIF
   2  CONTINUE
      NPURG=IQ(KQSP+LPURG+1)
      IF(NPURG.EQ.0)THEN
         IQ(KQSP+LPURG+1)=1
         IQ(KQSP+LPURG+2)=IPURG
         IQ(KQSP+LPURG+3)=IPURL
         GO TO 99
      ENDIF
*
      DO 5 I=1,NPURG
         IRF=IQ(KQSP+LPURG+2*I)
         IRL=IQ(KQSP+LPURG+2*I+1)
         IF(IPURG.GE.IRF.AND.IPURG.LE.IRL)IPURG=IRL+1
         IF(IPURL.GE.IRF.AND.IPURL.LE.IRL)IPURL=IRF-1
         IF(IPURG.GT.IPURL)GO TO 99
   5  CONTINUE
*
      IPU=2
  10  IF(IPURL.LT.IQ(KQSP+LPURG+IPU))THEN
         IF(IPURL.EQ.IQ(KQSP+LPURG+IPU)-1)THEN
            IQ(KQSP+LPURG+IPU)=IPURG
         ELSE
            IF(IPU.GT.2.AND.IPURG.LE.IQ(KQSP+LPURG+IPU-1))GO TO 99
            NDATA=IQ(KQSP+LPURG-1)
            IF(NDATA.LT.2*NPURG+3)THEN
               CALL MZPUSH(JQPDVS,LPURG,0,10,'I')
            ENDIF
            NLEFT=2*NPURG-IPU+2
            IF(NLEFT.GT.0)THEN
               CALL UCOPY2(IQ(KQSP+LPURG+IPU),
     +                     IQ(KQSP+LPURG+IPU+2),NLEFT)
            ENDIF
            NPURG=NPURG+1
            IQ(KQSP+LPURG+1)=NPURG
            IQ(KQSP+LPURG+IPU)=IPURG
            IQ(KQSP+LPURG+IPU+1)=IPURL
         ENDIF
         GO TO 99
      ENDIF
*
      IF(IPURG.EQ.IQ(KQSP+LPURG+IPU+1)+1)THEN
         IQ(KQSP+LPURG+IPU+1)=IPURL
         IF(IPU+2.LT.2*NPURG)THEN
            IF(IQ(KQSP+LPURG+IPU+1).EQ.IQ(KQSP+LPURG+IPU+2))THEN
               IQ(KQSP+LPURG+IPU+1)=IQ(KQSP+LPURG+IPU+3)
               NLEFT=2*NPURG-IPU-2
               IF(NLEFT.GT.0)THEN
                  CALL UCOPY2(IQ(KQSP+LPURG+IPU+4),
     +                        IQ(KQSP+LPURG+IPU+2),NLEFT)
               ENDIF
               NPURG=NPURG-1
               IQ(KQSP+LPURG+1)=NPURG
            ENDIF
         ENDIF
         GO TO 99
      ENDIF
*
      IPU=IPU+2
      IF(IPU.LE.2*NPURG)GO TO 10
*
      NDATA=IQ(KQSP+LPURG-1)
      IF(NDATA.LT.2*NPURG+3)THEN
         CALL MZPUSH(JQPDVS,LPURG,0,2,'I')
      ENDIF
      IQ(KQSP+LPURG+IPU)=IPURG
      IQ(KQSP+LPURG+IPU+1)=IPURL
      NPURG=NPURG+1
      IQ(KQSP+LPURG+1)=NPURG
*
  99  RETURN
      END
