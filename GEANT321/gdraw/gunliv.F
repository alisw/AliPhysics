*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GUNLIV(IVO,NR,IVOMOT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Unlink a volume from its mother                          *
C.    *                                                                *
C.    *       IVO = volume pointer (input)                             *
C.    *       NR = volume number (input)                               *
C.    *       IVOMOT = mother volume pointer (input)                   *
C.    *                                                                *
C.    *    ==>Called by : GEDITV                                       *
C.    *       Author    P.Zanarini  *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
C.
      JVOMOT=LQ(JVOLUM-IVOMOT)
      NIN=Q(JVOMOT+3)
      IF (NIN.LT.0) THEN
         Q(JVOMOT+3)=0
         CALL MZDROP(IXCONS,LQ(JVOMOT-1),' ')
         JV = LQ(JVOLUM-IVOMOT)
         CALL MZPUSH(IXCONS,JV,-1,0,'I')
      ELSE IF (NIN.EQ.1) THEN
         Q(JVOMOT+3)=0
         CALL MZDROP(IXCONS,LQ(JVOMOT-1),' ')
         JV = LQ(JVOLUM-IVOMOT)
         CALL MZPUSH(IXCONS,JV,-1,0,'I')
         JVOMOT=LQ(JVOLUM-IVOMOT)
         NL=IQ(JVOMOT-2)
         IF (NL.EQ.1) THEN
            CALL MZDROP(IXCONS,LQ(JVOMOT-1),' ')
            JV = LQ(JVOLUM-IVOMOT)
            CALL MZPUSH(IXCONS,JV,-1,0,'I')
         ENDIF
      ELSE IF (NIN.GT.1) THEN
         NL=IQ(JVOMOT-2)
         DO 10 IN=1,NL
            JIN=LQ(JVOMOT-IN)
            IF (Q(JIN+2).EQ.IVO.AND.Q(JIN+3).EQ.NR) GO TO 20
   10    CONTINUE
         GO TO 999
   20    Q(JVOMOT+3)=NIN-1
         DO 30 I=IN,NIN
            CALL DZSWAP(IXCONS,LQ(JVOMOT-I),LQ(JVOMOT-I-1),' ')
   30    CONTINUE
         CALL MZDROP(IXCONS,LQ(JVOMOT-NIN-1),' ')
         JV = LQ(JVOLUM-IVOMOT)
         CALL MZPUSH(IXCONS,JV,-1,0,'I')
         JVOMOT=LQ(JVOLUM-IVOMOT)
         JIN=LQ(JVOMOT-NIN)
         NIN=NIN-1
         DO 40 I=1,NIN
            Q(JIN+I)=I
            Q(JIN+I+NIN)=0
   40    CONTINUE
         JV = LQ(JVOMOT-NIN-1)
         CALL MZPUSH(IXCONS,JV,0,-2,'I')
      ENDIF
  999 RETURN
      END
