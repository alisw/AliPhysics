*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GPDIGI(IUSET,IUDET)
C.
C.    **************************************************************************
C.    *                                                                        *
C.    *       Print DIGIts in detector IUDET of set IUSET                      *
C.    *       (in case IUSET/IUDET = *, take all sets/detectors)               *
C.    *                                                                        *
C.    *       JDI=LQ(JDIGI-ISET)                                               *
C.    *       JDID=LQ(JDI-IDET)                                                *
C.    *       IQ(JDI+IDET)= pointer to LAST USED word in JDID                  *
C.    *                                                                        *
C.    *          Each digit is packed in JDID in the following format          *
C.    *          -- Track  numbers packed                                      *
C.    *          -- Volume numbers packed                                      *
C.    *          -- Digits packed                                              *
C.    *                                                                        *
C.    *    ==>Called by : <USER>, GPRINT                                       *
C.    *       Author    W.Gebel  *********                                     *
C.    *                                                                        *
C.    **************************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      PARAMETER (NDEMX=100,NVMAX=20)
      DIMENSION KDIGI(NDEMX),NUMBV(NVMAX),KWS(130),LTR(3)
      EQUIVALENCE (WS(1),NUMBV(1)),(WS(101),KDIGI(1))
      CHARACTER*4 IUSET,IUDET
C.    ------------------------------------------------------------------
C.
      IF(JDIGI.LE.0)GO TO 999
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
C             Loop on all selected sets
C
      DO 230 ISET=NS1,NS2
      JS=LQ(JSET-ISET)
      JDI=LQ(JDIGI-ISET)
      IF(JDI.LE.0)GO TO 230
      NDET=IQ(JS-1)
      ND1=1
      ND2=NDET
      IF(IUDET(1:1).NE.'*')THEN
         CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
         IF(IDET.EQ.0)GO TO 230
         ND1=IDET
         ND2=IDET
      ENDIF
C
C             Loop on selected detectors for this set
C
      DO 220 IDET=ND1,ND2
      JD=LQ(JS-IDET)
      JDID=LQ(JDI-IDET)
      IF(JDID.LE.0)GO TO 220
      JDDI=LQ(JD-2)
C
      WRITE(CHMAIL,1000)IQ(JS+IDET),IQ(JSET+ISET)
      CALL GMAIL(0,0)
C
C             Get volumes / digitisings names and print header line
C
      ILAST=IQ(JDI+IDET)
      IF(ILAST.EQ.0)GO TO 220
      NV=IQ(JD+2)
      ND=IQ(JD+6)
C
      CALL VBLANK(KWS,130)
      K=3
      IF(NV.GT.0)THEN
C        Number of printed elements limited to 15
         NVM=MIN(NV,15)
         DO 22 I=1,NVM
            CALL UBLOW(IQ(JD+2*I+9),KWS(K),4)
  22     K=K+5
         K=K+5
      ENDIF
      IF(ND.GT.0)THEN
         DO 26 I=1,ND
            IF(K.LE.101) CALL UBLOW(IQ(JDDI+2*I-1),KWS(K),4)
            K=K+8
  26     CONTINUE
         IF(K.GT.104)K=104
      ENDIF
      WRITE(CHMAIL,2000)(KWS(I),I=1,K)
      CALL GMAIL(0,1)
C
C             Now loop on all digits
C             to get track numbers, volume numbers and digits
C
      IDIG=0
      I=0
      NWDI=0
C
   30 CONTINUE
      I=I+NWDI
      IF(I.GE.ILAST)GO TO 220
      NK=2
      IDIG=IDIG+1
      CALL VZERO (LTR(1),3)
C
C             Get unpacked (first 3) tracks producing this digit
C             (2 tracks packed in 1 word; 1st half of 1st word: NTRA-1)
C
      NWDI=IQ(JDID+I+1)
      NTRM1=IBITS(IQ(JDID+I+NK), 0,16)
      NTRA=NTRM1+1
      IF(NTRA.GE.1)LTR(1)=IBITS(IQ(JDID+I+NK),16,16)
      NK=NK+1
      IF(NTRA.GE.2)LTR(2)=IBITS(IQ(JDID+I+NK), 0,16)
      IF(NTRA.GE.3)LTR(3)=IBITS(IQ(JDID+I+NK),16,16)
      NWTR=NTRA/2+1
      NK=NWTR+2
C
C             Get unpacked volume numbers
C
      IF(NV.GT.0)THEN
         K=1
         DO 50 IV=1,NV
            NB=IQ(JD+2*IV+10)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   NK=NK+1
               ENDIF
               NUMBV(IV)=IQ(JDID+I+NK)
               K=1
               NK=NK+1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               NUMBV(IV)=IBITS(IQ(JDID+I+NK),K-1,NB)
               K=K+NB
            ENDIF
  50     CONTINUE
         IF(K.NE.1)NK=NK+1
      ENDIF
C
C             Get unpacked digits
C
      IF(ND.GT.0)THEN
         K=1
         DO 90 ID=1,ND
            NB=IQ(JDDI+2*ID)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   NK=NK+1
               ENDIF
               IF(ID.LE.NDEMX) KDIGI(ID)=IQ(JDID+I+NK)
               K=1
               NK=NK+1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               IF(ID.LE.NDEMX) KDIGI(ID)=IBITS(IQ(JDID+I+NK),K-1,NB)
               K=K+NB
            ENDIF
  90     CONTINUE
      ENDIF
C
C             Do the printout
C             (fitting in 1 line of 128 characters per each digit)
C
      IF(NV.EQ.0)GO TO 119
      IF(NV.GT.15)NV=15
      GO TO (101,102,103,104,105,106,107,108,109,110
     +,      111,112,113,114,115), NV
C
  101 NDP=MIN(ND,12)
      WRITE(CHMAIL,3001)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  102 NDP=MIN(ND,11)
      WRITE(CHMAIL,3002)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  103 NDP=MIN(ND,10)
      WRITE(CHMAIL,3003)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  104 NDP=MIN(ND,10)
      WRITE(CHMAIL,3004)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  105 NDP=MIN(ND, 9)
      WRITE(CHMAIL,3005)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  106 NDP=MIN(ND, 8)
      WRITE(CHMAIL,3006)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  107 NDP=MIN(ND, 8)
      WRITE(CHMAIL,3007)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  108 NDP=MIN(ND, 7)
      WRITE(CHMAIL,3008)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  109 NDP=MIN(ND, 7)
      WRITE(CHMAIL,3009)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  110 NDP=MIN(ND, 6)
      WRITE(CHMAIL,3010)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  111 NDP=MIN(ND, 5)
      WRITE(CHMAIL,3011)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  112 NDP=MIN(ND, 5)
      WRITE(CHMAIL,3012)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  113 NDP=MIN(ND, 4)
      WRITE(CHMAIL,3013)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  114 NDP=MIN(ND, 4)
      WRITE(CHMAIL,3014)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  115 NDP=MIN(ND, 3)
      WRITE(CHMAIL,3015)IDIG,(LTR(J),J=1,3),(NUMBV(L),L=1,NV)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
  119 NDP=MIN(ND,12)
      WRITE(CHMAIL,3000)IDIG,(LTR(J),J=1,3)
     +,               (KDIGI(L),L=1,NDP)
      CALL GMAIL(0,0)
      GO TO 30
C
 220  CONTINUE
 230  CONTINUE
C
C
 1000 FORMAT(' =====>DIGITS OF DETECTOR ** ',A4,
     +' ** OF SET ** ',A4,' **')
 2000 FORMAT(' DIGIT    TR1  TR2  TR3  ',104A1)
 3000 FORMAT(1X,I5,2X,3I5,3X,             12(1X,I7))
 3001 FORMAT(1X,I5,2X,3I5,3X, 1(1X,I4),2X,12(1X,I7))
 3002 FORMAT(1X,I5,2X,3I5,3X, 2(1X,I4),2X,11(1X,I7))
 3003 FORMAT(1X,I5,2X,3I5,3X, 3(1X,I4),2X,10(1X,I7))
 3004 FORMAT(1X,I5,2X,3I5,3X, 4(1X,I4),2X,10(1X,I7))
 3005 FORMAT(1X,I5,2X,3I5,3X, 5(1X,I4),2X, 9(1X,I7))
 3006 FORMAT(1X,I5,2X,3I5,3X, 6(1X,I4),2X, 8(1X,I7))
 3007 FORMAT(1X,I5,2X,3I5,3X, 7(1X,I4),2X, 8(1X,I7))
 3008 FORMAT(1X,I5,2X,3I5,3X, 8(1X,I4),2X, 7(1X,I7))
 3009 FORMAT(1X,I5,2X,3I5,3X, 9(1X,I4),2X, 7(1X,I7))
 3010 FORMAT(1X,I5,2X,3I5,3X,10(1X,I4),2X, 6(1X,I7))
 3011 FORMAT(1X,I5,2X,3I5,3X,11(1X,I4),2X, 5(1X,I7))
 3012 FORMAT(1X,I5,2X,3I5,3X,12(1X,I4),2X, 5(1X,I7))
 3013 FORMAT(1X,I5,2X,3I5,3X,13(1X,I4),2X, 4(1X,I7))
 3014 FORMAT(1X,I5,2X,3I5,2X,14(1X,I4),2X, 4(1X,I7))
 3015 FORMAT(1X,I5,2X,3I5,3X,15(1X,I4),2X, 3(1X,I7))
 999  RETURN
      END
