*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZSAVE
*
************************************************************************
*
*        Write all directories which have been modified in memory
*        Write current output buffer
*        Update list of used/unused records in top-directory
*
* Called by <USER>,RZCDIR,RZCOPY,RZEND,RZFILE,RZMAKE
*
*  Author  : R.Brun DD/US/PD
*  Written : 02.04.86
*  Last mod: 04.10.90
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbit.inc"
      IF(LQRS.EQ.0)GO TO 99
      IF(LTOP.EQ.0)GO TO 99
*
*           Mark used records in BITMAP
*
      IF(JBIT(IQ(KQSP+LTOP),2).NE.0)THEN
         IF(ISAVE.NE.2)THEN
            IDTIME=0
            CALL RZDATE(IDTIME,IDATE,ITIME,2)
            IQ(KQSP+LTOP+KDATEM)=IDTIME
         ENDIF
         LUNC= IQ(KQSP+LTOP-5)
         LB  = IQ(KQSP+LTOP+KLB)
         LREK= IQ(KQSP+LTOP+LB+1)
         LUS = LQ(KQSP+LTOP-3)
         IF(LUS.NE.0)THEN
            NUSED=IQ(KQSP+LUS+1)
            IF(NUSED.GT.0)THEN
               DO 40 I=1,NUSED
                  IR1=IQ(KQSP+LUS+2*(I-1)+2)
                  IRL=IQ(KQSP+LUS+2*(I-1)+3)
                  DO 30 J=IR1,IRL
                     IWORD = (J-1)/32 + 1
                     IBIT  = J-32*(IWORD-1)
                     CALL SBIT1(IQ(KQSP+LTOP+LB+2+IWORD),IBIT)
  30              CONTINUE
  40           CONTINUE
               IQ(KQSP+LUS+1)=0
            ENDIF
         ENDIF
*
*           Mark purged records in BITMAP
*
         LPU = LQ(KQSP+LTOP-5)
         IF(LPU.NE.0)THEN
            NPURG=IQ(KQSP+LPU+1)
            IF(NPURG.GT.0)THEN
               DO 60 I=1,NPURG
                  IR1=IQ(KQSP+LPU+2*(I-1)+2)
                  IRL=IQ(KQSP+LPU+2*(I-1)+3)
                  DO 50 J=IR1,IRL
                     IWORD = (J-1)/32 + 1
                     IBIT  = J-32*(IWORD-1)
                     CALL SBIT0(IQ(KQSP+LTOP+LB+2+IWORD),IBIT)
  50              CONTINUE
  60           CONTINUE
               IQ(KQSP+LPU+1)=0
            ENDIF
         ENDIF
*
*               Write current buffer
*
         LROUT=LQ(KQSP+LTOP-6)
         IF(LROUT.NE.0)THEN
            IROUT=IQ(KQSP+LTOP+KIROUT)
            IF(IROUT.NE.0)THEN
               CALL RZIODO(LUNC,LREK,IROUT,IQ(KQSP+LROUT+1),2)
               IF(IQUEST(1).NE.0)GO TO 99
            ENDIF
         ENDIF
*
*               Write TOP directory
*
         LDS =IQ(KQSP+LTOP+KLD)
         NRD =IQ(KQSP+LTOP+LDS)
         IF(ISAVE.NE.2)THEN
            IF(LTOP.EQ.LCDIR)IQ(KQSP+LTOP+KDATEM)=IDTIME
         ENDIF
         CALL SBIT0(IQ(KQSP+LTOP),2)
         DO 70 J=NRD,1,-1
            IREC=IQ(KQSP+LTOP+LDS+J)
            L=(J-1)*LREK+1
            CALL RZIODO(LUNC,LREK,IREC,IQ(KQSP+LTOP+L),2)
            IF(IQUEST(1).NE.0)THEN
               CALL SBIT1(IQ(KQSP+LTOP),2)
               GO TO 99
            ENDIF
  70     CONTINUE
*
*               Write current directory if modified
*
         IF(LCDIR.EQ.0.OR.LTOP.EQ.LCDIR)GO TO 99
         IF(JBIT(IQ(KQSP+LCDIR),2).NE.0)THEN
            LDS =IQ(KQSP+LCDIR+KLD)
            NRD =IQ(KQSP+LCDIR+LDS)
            IF(ISAVE.NE.2)THEN
               IQ(KQSP+LCDIR+KDATEM)=IDTIME
            ENDIF
            CALL SBIT0(IQ(KQSP+LCDIR),2)
            DO 80 J=NRD,1,-1
               IREC=IQ(KQSP+LCDIR+LDS+J)
               L=(J-1)*LREK+1
               CALL RZIODO(LUNC,LREK,IREC,IQ(KQSP+LCDIR+L),2)
               IF(IQUEST(1).NE.0)THEN
                  CALL SBIT1(IQ(KQSP+LCDIR),2)
                  GO TO 99
               ENDIF
  80        CONTINUE
         ENDIF
      ENDIF
*
  99  RETURN
      END
