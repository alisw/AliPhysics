*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:16  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZVER1(CHL,CHOPT,IRET)
************************************************************************
*
*     Slave routine to RZVERI
*
*
*  Author  : J.Shiers CN/AS/DL
*  Written : 23.03.92
*  Last mod: 18.04.94 - set ISTAT2
*          : 02.02.95 - cater for new RZ file format
*
************************************************************************
*
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzbmap.inc"
#include "zebra/rzover.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*(*) CHL,CHOPT
      CHARACTER*255 DIRNAM
      DIMENSION KEY(KNMAX)
*...............................................................

#include "zebra/q_jbit.inc"
#include "zebra/q_jbyt.inc"

      NCHL   = LENOCC(CHL)
      LOGLV  = JBYT(IQ(KQSP+LTOP),15,3)-3
      IF(LOGLV.GE.1) WRITE(IQPRNT,*)
     +   'Processing directory: ',CHL(1:NCHL)
      DIRNAM = CHL
      IRET   = 0
      JRET   = 0
      NKEYS  = IQ(KQSP+LCDIR+KNKEYS)
      NWKEY  = IQ(KQSP+LCDIR+KNWKEY)
      LB     = IQ(KQSP+LTOP+KLB)
      LK     = IQ(KQSP+LCDIR+KLK)
      LDS    = IQ(KQSP+LCDIR+KLD)
      LREC   = IZRECL
      NCHO   = LENOCC(CHOPT)
 
      IOPTB  = INDEX(CHOPT(1:NCHO),'B')
      IOPTC  = INDEX(CHOPT(1:NCHO),'C')
      IOPTO  = INDEX(CHOPT(1:NCHO),'O')
      IOPTP  = INDEX(CHOPT(1:NCHO),'P')
*
*     Process all records of this directory
*
      NRD      = IQ(KQSP+LCDIR+LDS)
      DO 10 I  = 1,NRD
 
         ISTAT = 0
         IREC  = IQ(KQSP+LCDIR+LDS+I)
         IWORD = (IREC-1)/32 + 1
         IBIT  = IREC-32*(IWORD-1)
*
*     Print directory name and record number if rec. no in list
*     of bad records
*
          IF(IOPTP.NE.0.AND.NBAD.GT.0) THEN
             IF(IUFIND(IREC,IBAD,1,NBAD).LE.NBAD) WRITE(IQPRNT,*)
     +          'Directory: ',CHL(1:LENOCC(CHL)),' uses record ',IREC
          ENDIF
*
*     Is this record marked as free?
*
         IF(IOPTC.NE.0.AND.
     +      JBIT(IQ(KQSP+LTOP+LB+2+IWORD),IBIT).EQ.0) THEN
            WRITE(IQPRNT,*) 'RZVER1. warning - record ',IREC,
     +         ' is in use but is marked as free in bit map'
            JRET = JRET + 1
         ENDIF
*
*     Set bit to mark record as used
*
         IF(IOPTB.NE.0) CALL SBIT1(IQ(KQSP+LTOP+LB+2+IWORD),IBIT)
         IW1=(IREC-1)*LREC+1
         IW2=IW1+LREC-1
*
*     Check for overwriting at the directory level
*
         IF(IOPTO.NE.0) CALL RZVER2(IW1,IW2,ISTAT)
         IRET=IRET+ISTAT
   10 CONTINUE
      IF(IRET.NE.0)THEN
         DIRNAM=CHL
         NCHL=LENOCC(CHL)
         PRINT 10000, DIRNAM(1:NCHL),(IQ(KQSP+LCDIR+LDS+I),I=1,NRD)
*
*     Store record numbers for second pass
*
         IF(NBAD+NRD.LE.MAXBAD) THEN
            DO 20 I=1,NRD
               IBAD(NBAD+I) = IQ(KQSP+LCDIR+LDS+I)
   20       CONTINUE
            NBAD            = NBAD + NRD
         ENDIF
 
      ENDIF
*
*     Check records used for objects in this directory
*
      IF(NKEYS.GT.0)THEN
         DO 60 I=1,NKEYS
            ISTAT = 0
            LKC   = LK+(NWKEY+1)*(I-1)
            LCYC  = IQ(KQSP+LCDIR+LKC)
 
   30       CONTINUE
 
            IF(KVSCYC.EQ.0) THEN
               LCOLD = JBYT(IQ(KQSP+LCDIR+LCYC), 1,16)
               IR1   = JBYT(IQ(KQSP+LCDIR+LCYC+2),17,16)
               IR2   = JBYT(IQ(KQSP+LCDIR+LCYC ),17,16)
               IP1   = JBYT(IQ(KQSP+LCDIR+LCYC+2), 1,16)
               NW    = JBYT(IQ(KQSP+LCDIR+LCYC+3), 1,20)
            ELSE
               LCOLD = IQ(KQSP+LCDIR+LCYC)
               IR1   = IQ(KQSP+LCDIR+LCYC+2)
               IR2   = IQ(KQSP+LCDIR+LCYC+5)
               IP1   = JBYT(IQ(KQSP+LCDIR+LCYC+3),1,20)
               NW    = IQ(KQSP+LCDIR+LCYC+4)
            ENDIF
 
            NLEFT=LREC-IP1+1
            IW1=(IR1-1)*LREC+IP1
            IW2=(IR1-1)*LREC+IP1+MIN(NLEFT,NW)-1
*
*     Check for overwriting at the object level
*
            IF(IOPTO.NE.0) THEN
               CALL RZVER2(IW1,IW2,ISTAT)
               IF(ISTAT.NE.0.AND.NBAD.LT.MAXBAD) THEN
                  IBAD(NBAD+1) = IR1
                  NBAD         = NBAD + 1
                ENDIF
            ENDIF
 
            IWORD = (IR1-1)/32 + 1
            IBIT  = IR1-32*(IWORD-1)
*
*     Is this record marked as free?
*
            IF(IOPTC.NE.0.AND.JBIT(IQ(KQSP+LTOP+LB+2+IWORD),IBIT)
     +      .EQ.0) THEN
               WRITE(IQPRNT,*) 'RZVER1. warning - record ',IR1,' is in '
     +         //'use but is marked as free in bit map'
               JRET = JRET + 1
            ENDIF
*
*     Set bit to mark record as used
*
            IF(IOPTB.NE.0) CALL SBIT1(IQ(KQSP+LTOP+LB+2+IWORD),IBIT)
*
*     Is this record in the list of overwritten records?
*
            IF(IOPTP.NE.0.AND.NBAD.GT.0) THEN
               IF(IUFIND(IR1,IBAD,1,NBAD).LE.NBAD) THEN
                  NCHL = LENOCC(CHL)
                  WRITE(IQPRNT,*) 'Directory: ',CHL(1:NCHL),
     +            ' has objects in record ',IR1
                  CALL RZPRNK(CHL(1:NCHL),I,LCYC,' ')
               ENDIF
            ENDIF
*
*     Any more records for this object?
*
            IF(NW.GT.NLEFT)THEN
               NR=(NW-NLEFT-1)/LREC
               IRN=IR2+NR
               DO 40 J=IR2,IRN
                  IF(IOPTP.NE.0.AND.NBAD.GT.0) THEN
                     IF(IUFIND(J,IBAD,1,NBAD).LE.NBAD) THEN
                        NCHL = LENOCC(CHL)
                        WRITE(IQPRNT,*) 'Directory: ',
     +                    CHL(1:NCHL),' has objects in record ',IR1
                          CALL RZPRNK(CHL(1:NCHL),I,LCYC,' ')
                     ENDIF
                  ENDIF
 
                  IWORD = (J-1)/32 + 1
                  IBIT  = J-32*(IWORD-1)
*
*     Is this record marked as free?
*
                  IF(IOPTC.NE.0.AND.JBIT(IQ(KQSP+LTOP+LB+2+IWORD),
     +            IBIT).EQ.0) THEN
                     WRITE(IQPRNT,*) 'RZVER1. warning - record ',J,
     +               ' is in use but is marked as free in bit map'
                     JRET = JRET + 1
                  ENDIF
*
*     Set bit to mark record as used
*
                  IF(IOPTB.NE.0) CALL SBIT1(IQ(KQSP+LTOP+LB+2+IWORD),
     +            IBIT)
   40          CONTINUE
               IW1=(IR2-1)*LREC+1
               IW2=IW1+NW-NLEFT-1
               ISTAT2=0
               IF(IOPTO.NE.0) CALL RZVER2(IW1,IW2,ISTAT2)
               ISTAT=ISTAT+ISTAT2
            ENDIF
            IF(ISTAT2.NE.0.AND.NW.GT.NLEFT)THEN
               DO 50 J=IR2,IRN
                  IF(NBAD.LT.MAXBAD) THEN
                     IBAD(NBAD+1) = J
                     NBAD         = NBAD + 1
                  ENDIF
   50          CONTINUE
               ICYC = JBYT(IQ(KQSP+LCDIR+LCYC+3),21,12)
               DIRNAM=CHL
               NCHL=LENOCC(CHL)
               IRET=IRET+1
*
*     Get and print key of corrupted object
*
               PRINT 10100, DIRNAM(1:NCHL),IR1,((NW-1)/LREC)+1
               CALL RZPRNK(DIRNAM(1:NCHL),I,LCYC,' ')
            ENDIF
            IF(LCOLD.NE.0.AND.LCOLD.NE.LCYC)THEN
               LCYC=LCOLD
               GO TO 30
            ENDIF
   60    CONTINUE
      ENDIF
 
      IQUEST(2) = JRET
*
10000 FORMAT(' **** WARNING: Directory ',A,' possibly overwritten ****',
     +       /,' records numbers: ',/10(1X,I6))
10100 FORMAT(' **** WARNING: Object in directory ',A,' corrupted ****',
     +       /,' start record: ',I6,' number of records: ',I6)
*
      END
