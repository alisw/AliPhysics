*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GPSETS(IUSET,IUDET)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints set and detector parameters                       *
C.    *                                                                *
C.    *        IUSET  user set identifier                              *
C.    *               If * prints all detectors of all sets            *
C.    *        IUDET  user detector identifier                         *
C.    *               If * prints all detectors of set IUSET           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT   ,<GXINT> GINC4              *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 IUSET,IUDET
C.
C.    ------------------------------------------------------------------
C.
      IF(IUSET(1:1).EQ.'*')THEN
         WRITE(CHMAIL,1000)
         CALL GMAIL(0,0)
      ENDIF
      IF(JSET.LE.0)GO TO 999
      NSET=IQ(JSET-1)
      NS1=1
      NS2=NSET
      IF(IUSET(1:1).NE.'*')THEN
         CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
         IF(ISET.LE.0)GO TO 999
         NS1=ISET
         NS2=ISET
      ENDIF
C
C             loop on all selected sets
C
      DO 130 ISET=NS1,NS2
      JS=LQ(JSET-ISET)
      IF(JS.LE.0)GO TO 130
      NDET=IQ(JS-1)
      ND1=1
      ND2=NDET
      IF(IUDET(1:1).NE.'*')THEN
         CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
         IF(IDET.EQ.0)GO TO 130
         ND1=IDET
         ND2=IDET
      ENDIF
C
C             loop on selected detectors for this set
C
      DO 120 IDET=ND1,ND2
      JD=LQ(JS-IDET)
      IF(JD.LE.0)GO TO 120
      NV=IQ(JD+2)
      NWHI=IQ(JD+7)
      NWDI=IQ(JD+8)
C
      WRITE(CHMAIL,2000)IQ(JSET+ISET),IQ(JS+IDET),NWHI,NWDI
      CALL GMAIL(0,0)
C
      IF(NV.NE.0)THEN
         DO 10 I=1,NV
            WRITE(CHMAIL,3000)IQ(JD+2*I+9),IQ(JD+2*I+10)
            CALL GMAIL(0,0)
  10     CONTINUE
      ENDIF
C
      JDH=LQ(JD-1)
      IF(JDH.GT.0)THEN
         NH=IQ(JD+4)
         IF(NH.GT.0)THEN
            DO 30 I=1,NH
               WRITE(CHMAIL,4000)IQ(JDH+4*I-3),IQ(JDH+4*I-2),
     +                           Q(JDH+4*I-1),Q(JDH+4*I)
               CALL GMAIL(0,0)
  30        CONTINUE
         ENDIF
      ENDIF
C
      JDD=LQ(JD-2)
      IF(JDD.GT.0)THEN
         ND=IQ(JD+6)
         IF(ND.GT.0)THEN
            DO 50 I=1,ND
               WRITE(CHMAIL,5000)IQ(JDD+2*I-1),IQ(JDD+2*I)
               CALL GMAIL(0,0)
  50        CONTINUE
         ENDIF
      ENDIF
C
      JDU=LQ(JD-3)
      IF(JDU.GT.0)THEN
         NU=IQ(JDU-1)
         IF(NU.GT.0)THEN
            WRITE(CHMAIL,6000)
            CALL GMAIL(0,0)
            DO 60 I=1,NU,10
               L2=I+9
               IF(L2.GT.NU)L2=NU
               WRITE(CHMAIL,7000)(Q(JDU+L),L=I,L2)
               CALL GMAIL(0,0)
  60        CONTINUE
         ENDIF
      ENDIF
C
 120  CONTINUE
 130  CONTINUE
C
 1000 FORMAT('0',51('='),3X,'SETS AND DETECTORS',3X,50('='))
 2000 FORMAT('   SET ',A4,'  DETECTOR ',A4,'  NWHI=',I6,'  NWDI=',I6)
 3000 FORMAT(10X,'VOLUME ',A4,'  NBITSD=',I4)
 4000 FORMAT(10X,'HIT ELEMENT =',A4,'  NBITSH=',I4,
     +'  ORIG =',E12.4,'  FACT =',E12.4)
 5000 FORMAT(10X,'DIGIT ELEMENT =',A4,'  NBITSD=',I4)
 6000 FORMAT(10X,'USER PARAMETERS')
 7000 FORMAT(5X,10E12.4)
 999  RETURN
      END
