*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR1(JM,J1,J2,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Scan one level of JVOLUM structure                       *
C.    *                                                                *
C.    *       JM = mother node (input)                                 *
C.    *       J1 = starting node (input)                               *
C.    *       J2 = ending node (output)                                *
C.    *       IER  = error flag to detect nodes overflow (output)      *
C.    *                                                                *
C.    *    ==>Called by : GDTR2                                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 NAME, ISON
      DIMENSION N(50),IH(50)
      SAVE NUMCHK,MAXCHK
      DATA NUMCHK/4/,MAXCHK/50/
C.
C.    ------------------------------------------------------------------
C.
      IER=0
C
      IF (NUMCHK.GT.MAXCHK) GO TO 140
C
      CALL UHTOC (IQ(JNAM+JM), 4, NAME, 4)
      CALL GDNSON (NAME, NXONS, MUL)
      IF (NXONS.EQ.0) THEN
         IQ(JXON+JM)=0
         J2=J1
      ELSE
         J2=J1+1
         IQ(JXON+JM)=J2
         CALL UHTOC (IQ(JNAM+JM), 4, NAME, 4)
         CALL GDSON (1, NAME, ISON)
         CALL UCTOH (ISON, LXON, 4, 4)
         LOXON=LXON
         IQ(JNAM+J2)=LXON
         IQ(JMOT+J2)=JM
         IF (NXONS.GT.1) THEN
            IF (J2+NXONS.GE.MAXNOD) GO TO 130
C
            JSAV=J2
C
            DO 10 I=2,NXONS
               CALL UHTOC (IQ(JNAM+JM), 4, NAME, 4)
               CALL GDSON (I, NAME, ISON)
               CALL UCTOH (ISON, LXON, 4, 4)
               IF (LXON.NE.LOXON) THEN
                  LOXON=LXON
                  J=J2+1
                  IQ(JBRO+J2)=J
                  IQ(JNAM+J)=LXON
                  IQ(JMOT+J)=JM
                  J2=J
               ELSE
                  IQ(JPSM+J2)=IQ(JPSM+J2)+1
               ENDIF
   10       CONTINUE
C
            DO 20 K=JSAV,J2
               CALL UHTOC (IQ(JNAM+K), 4, NAME, 4)
               CALL GFATT (NAME, 'SEEN', KVAL)
               IF (KVAL.EQ.-3) GO TO 30
   20       CONTINUE
C
            GO TO 120
C
   30       I=JSAV
C
   40       CONTINUE
C
            DO 50  K=1,NUMCHK
               N(K)=IQ(JNAM+I+K-1)
               IH(K)=N(K)/65536
#if defined(CERNLIB_VAX)||defined(CERNLIB_MSDOS)||defined(CERNLIB_WINNT)
               IH(K)=N(K)-IH(K)*65536
#endif
               IF (IH(1).NE.IH(K)) GO TO 90
   50       CONTINUE
            DO 70  K=1,NUMCHK
               DO 60  KK=K+2,NUMCHK
                  IF (IQ(JNAM+K).EQ.IQ(JNAM+KK)) GO TO 90
   60          CONTINUE
   70       CONTINUE
            DO 80  K=I+2,J2
               IQ(JNAM+K-1)=IQ(JNAM+K)
   80       CONTINUE
            J2=J2-1
            IQ(JPSM+I+1)=IQ(JPSM+I+1)+1
            GO TO 100
   90       CONTINUE
            I=I+1
  100       CONTINUE
C
            IF (J2.GT.I+NUMCHK-2) GO TO 40
C
            DO 110 K=JSAV+1,J2-1
            IF (IQ(JPSM+K).GT.1) CALL UCTOH('    ',IQ(JNAM+K),4,4)
  110       CONTINUE
C
  120       CONTINUE
C
         ENDIF
         IQ(JBRO+J2)=0
      ENDIF
      IQ(JSCA+JM)=1
      GO TO 999
  130 WRITE (CHMAIL,1000) MAXNOD
      CALL GMAIL(0,0)
      IER=1
      GO TO 999
  140 WRITE (CHMAIL,1100)
      CALL GMAIL(0,0)
      IER=1
 1000 FORMAT (' GDTR1 : MAXNOD = ',I5,' - TREE OVERFLOW')
 1100 FORMAT (' GDTR1 : NUMCHK > MAXCHK - TREE NOT EXECUTED')
  999 RETURN
      END
