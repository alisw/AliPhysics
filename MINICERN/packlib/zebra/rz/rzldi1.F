*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:58  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZLDI1(IOPTA,IOPTX)
*
************************************************************************
*
*       List current directory
* Input:
*   IOPTA   =1 to list All objects
*   IOPTX   =1 for eXtended listing (default 80 columns)
*
* Called by <RZLDIR>
*
*  Author  : R.Brun DD/US/PD
*  Written : 19.01.89
*  Last mod: 02.12.92 JDS. Increase IOFF2 to 50
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*          : 27.09.95 JDS. Increase width in format 10200
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzdir.inc"
#include "zebra/rzch.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*160 CHLD
      CHARACTER*1 CK
      DIMENSION    IHTAG(40),IHDIR(4)
#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

*
*-----------------------------------------------------------------------
*
 
      IOFF1 = 6
      IF(IOPTX.EQ.0) THEN
         IOFF2 = 50
         IOFF3 = 70
         IOFF4 = 79
      ELSE
         IOFF2 = 57
         IOFF3 = 120
         IOFF4 = 132
      ENDIF
*
*            General case
*
*
      LS=IQ(KQSP+LCDIR+KLS)
      IF(IQUEST(1).EQ.0)THEN
         CALL RZDATE(IQ(KQSP+LCDIR+KDATEC),IDATEC,ITIMEC,1)
         CALL RZDATE(IQ(KQSP+LCDIR+KDATEM),IDATEM,ITIMEM,1)
         DO 10 I=LEN(CHL),1,-1
            IF(CHL(I:I).NE.' ')GO TO 20
   10    CONTINUE
   20    WRITE(IQPRNT,10000)CHL(1:I),IDATEC,ITIMEC,IDATEM,ITIMEM
*
         NSDIR=IQ(KQSP+LCDIR+KNSD)
         IF(NSDIR.GT.0)THEN
            WRITE(IQPRNT,10100)
            DO 30 I=1,NSDIR
               IH=LS+7*(I-1)
               IF (KVSCYC.EQ.0) THEN
                  IRECS = JBYT(IQ(KQSP+LCDIR+IH+5), 1,18)
               ELSE
                  IRECS = IQ(KQSP+LCDIR+IH+5)
               ENDIF
               CALL RZDATE(IQ(KQSP+LCDIR+IH+6),IDATE,ITIME,1)
               CALL ZITOH(IQ(KQSP+LCDIR+IH),IHDIR,4)
               WRITE(IQPRNT,10200)IHDIR,IDATE,ITIME,IRECS
   30       CONTINUE
         ENDIF
*
         XNUS = 0.
         NKEYS= IQ(KQSP+LCDIR+KNKEYS)
         NWKEY= IQ(KQSP+LCDIR+KNWKEY)
         LK   = IQ(KQSP+LCDIR+KLK)
         NWK1 = NWKEY
         IF(NWKEY.GT.20)THEN
            WRITE(IQPRNT,11400) NWKEY
            NWKEY = 20
         ENDIF
         NWK2 = 2*NWKEY
         KTAGS= KKDES+(NWK1-1)/10+1
         CALL ZITOH(IQ(KQSP+LCDIR+KTAGS),IHTAG,NWK2)
         WRITE(IQPRNT,10300)
         CHLD=' '
         IF( NWKEY .LE. 10 ) THEN
            WRITE(CHLD,10400)(IHTAG(I),I=1,NWK2)
         ELSE
            WRITE(CHLD,10400)(IHTAG(I),I=1,NWKEY)
            WRITE(IQPRNT,'(A)')CHLD(1:106)
            CHLD = ' '
            WRITE(CHLD,10400)(IHTAG(I),I=NWKEY+1,NWK2)
         ENDIF
         IC1 = LENOCC(CHLD) + 2
         IF(IC1.GT.IOFF2) THEN
            WRITE(IQPRNT,'(A)') CHLD(1:IC1)
            CHLD = ' '
         ENDIF
         IC1 = IOFF2
         IF(IOPTX.EQ.0) THEN
            WRITE(CHLD(IC1:IOFF4),10500)
         ELSE
            WRITE(CHLD(IC1:IOFF4),10600)
         ENDIF
         WRITE(IQPRNT,'(A)') CHLD(1:IOFF4)
         CHLD=' '
         IF(NKEYS.GT.0)THEN
            DO 60 I=1,NKEYS
               LKC=LK+(NWK1+1)*(I-1)
               LCYC=IQ(KQSP+LCDIR+LKC)
   40          IF (KVSCYC.EQ.0) THEN
                  LCOLD = JBYT(IQ(KQSP+LCDIR+LCYC+KPPCYC), 1,16)
                  IR1   = JBYT(IQ(KQSP+LCDIR+LCYC+KFRCYC),17,16)
                  IR2   = JBYT(IQ(KQSP+LCDIR+LCYC+KSRCYC),17,16)
                  IP1   = JBYT(IQ(KQSP+LCDIR+LCYC+KORCYC), 1,16)
                  NW    = JBYT(IQ(KQSP+LCDIR+LCYC+KNWCYC), 1,20)
               ELSE
                  LCOLD = IQ(KQSP+LCDIR+LCYC+KPPCYC)
                  IR1   = IQ(KQSP+LCDIR+LCYC+KFRCYC)
                  IR2   = IQ(KQSP+LCDIR+LCYC+KSRCYC)
                  IP1   = JBYT(IQ(KQSP+LCDIR+LCYC+KORCYC), 1,20)
                  NW    = IQ(KQSP+LCDIR+LCYC+KNWCYC)
               ENDIF
               CALL RZDATE( IQ(KQSP+LCDIR+LCYC+KFLCYC),IDATE,ITIME,1)
               IBIT4 = JBIT(IQ(KQSP+LCDIR+LCYC+KFLCYC),4)
               IBIT5 = JBIT(IQ(KQSP+LCDIR+LCYC+KFLCYC),5)
               IC    = JBYT(IQ(KQSP+LCDIR+LCYC+KCNCYC),21,12)
               IF(IBIT5.EQ.0)THEN
                  CK=' '
               ELSE
                  CK='*'
               ENDIF
               XNUS  = XNUS+NW
               NLEFT=LREC-IP1+1
               IF(NW.LE.NLEFT)THEN
                  NR=0
               ELSE
                  NR=(NW-NLEFT-1)/LREC+1
               ENDIF
               CHLD=' '
               IC1=IOFF1
               IFMOLD=0
               DO 50 J=1,NWKEY
                  IF(IC1.GT.IOFF3) THEN
                     WRITE(IQPRNT,'(A)')CHLD(1:IC1)
                     CHLD=' '
                     IC1=IOFF1
                     IFMOLD=0
                  ENDIF
                  IKDES=(J-1)/10
                  IKBIT1=3*J-30*IKDES-2
                  IFORM=JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3)
                  IF(IFORM.LT.3)THEN
                     KEY(J)=IQ(KQSP+LCDIR+LKC+J)
                     IF(IFMOLD.NE.0)IC1=IC1+10
                     IF(IFMOLD.EQ.4)IC1=IC1+6
                     IF(IABS(KEY(J)).LT.100000)THEN
                        WRITE(CHLD(IC1:),'(I6)')KEY(J)
                     ELSE
                        WRITE(CHLD(IC1-4:),'(I10)')KEY(J)
                     ENDIF
                  ELSE
                     CALL ZITOH(IQ(KQSP+LCDIR+LKC+J),KEY(J),1)
                     IF(IFORM.EQ.3)THEN
                        IF(IFMOLD.NE.0)IC1=IC1+10
                        IF(IFMOLD.EQ.4)IC1=IC1+6
                        WRITE(CHLD(IC1:),10700)KEY(J)
                     ELSE
                        IF(IFMOLD.EQ.4)IC1=IC1+4
                        IF(IFMOLD.NE.4.AND.IFMOLD.NE.0)IC1=IC1+10
                        WRITE(CHLD(IC1:),'(A4)')KEY(J)
                     ENDIF
                  ENDIF
                  IFMOLD=IFORM
   50          CONTINUE
               IC1=LENOCC(CHLD) + 2
               IF(IC1.GT.IOFF2) THEN
                  WRITE(IQPRNT,'(A)') CHLD(1:IC1)
                  CHLD = ' '
               ENDIF
               IC1 = IOFF2
               IF(IOPTX.EQ.0) THEN
                  WRITE(CHLD(IC1:),10800)IC,CK,IDATE,ITIME,NW
               ELSE
                  IF(IR2.EQ.0)THEN
                     WRITE(CHLD(IC1:),10900)IC,CK,IDATE,ITIME,NW,IP1,
     +               IR1
                  ENDIF
                  IF(NR.EQ.1)THEN
                     WRITE(CHLD(IC1:),11000)IC,CK,IDATE,ITIME,NW,IP1,
     +               IR1, IR2
                  ENDIF
                  IF(NR.GT.1)THEN
                     IRL=IR2+NR-1
                     WRITE(CHLD(IC1:),11100)IC,CK,IDATE,ITIME,NW,IP1,
     +               IR1, IR2,IRL
                  ENDIF
               ENDIF
               IF(IOPTA.NE.0.OR.IBIT4.EQ.0)
     +            WRITE(IQPRNT,'(A)')CHLD(1:IOFF4)
               IF(LCOLD.NE.0.AND.LCOLD.NE.LCYC)THEN
                  LCYC=LCOLD
                  GO TO 40
               ENDIF
   60       CONTINUE
         ENDIF
*
*            Print statistics
*
   70    WRITE(IQPRNT,11200)IQ(KQSP+LCDIR+KRUSED),IQ(KQSP+LCDIR+KMEGA),
     +   IQ(KQSP+LCDIR+KWUSED)
         XTOT =IQ(KQSP+LCDIR+KRUSED)*LREC
         IF(XTOT.GT.0.)THEN
            PERB =100.*XNUS/XTOT
         ELSE
            PERB=0.
         ENDIF
         PERD =100.*(FLOAT(IQ(KQSP+LCDIR+KRUSED)))/
     +         FLOAT(IQ(KQSP+LCDIR+KQUOTA))
         PERF =100.*(FLOAT(IQ(KQSP+LCDIR+KRUSED)))/
     +         FLOAT(IQ(KQSP+LTOP +KQUOTA))
         WRITE(IQPRNT,11300)PERD,PERF,PERB
      ENDIF
*
10000 FORMAT(///,' ************** Directory ===> ',A,' <===',//,
     +17X,' Created ',I6,'/',I4,'  Modified ',I6,'/',I4,/)
10100 FORMAT(/,' ===> List of subdirectories ')
10200 FORMAT(1X,4A4,' Created ',I6,'/',I4,' at record ',I10)
10300 FORMAT(/,' ===> List of objects ')
10400 FORMAT(5X,10(A4,A4,2X))
10500 FORMAT('   CYCLE    DATE/TIME   NDATA')
10600 FORMAT('   CYCLE    DATE/TIME   NDATA  OFFSET    REC1    REC2')
10700 FORMAT(2X,A4)
10800 FORMAT(4X,I4,A,2X,I6,'/',I4.4,1X,I6)
10900 FORMAT(4X,I4,A,2X,I6,'/',I4.4,1X,I6,2X,I6,2X,I6)
11000 FORMAT(4X,I4,A,2X,I6,'/',I4.4,1X,I6,2X,I6,2X,I6,2X,I6)
11100 FORMAT(4X,I4,A,2X,I6,'/',I4.4,1X,I6,2X,I6,2X,I6,2X,I6,' ==>',I6)
11200 FORMAT(/,'  Number of records =',I5,'  Number of megawords =',
     +       I3,' +',I6,' words')
11300 FORMAT('  Per cent of directory quota used =',F8.3,/
     +       '  Per cent of file used            =',F8.3,/
     +       '  Blocking factor                  =',F8.3)
11400 FORMAT(' RZLDI1. Cannot list more than 20 keys',/
     +       ' Number of keys found:',I10)
   80 RETURN
      END
