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
      SUBROUTINE GPHITS(IUSET,IUDET)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints HITS into detector IUDET of set IUSET             *
C.    *                                                                *
C.    *       If IPKHIT in /CGDRAW/ is greater than zero,              *
C.    *       only hit nr. IPKHIT is printed.                          *
C.    *                                                                *
C.    *       JH=LQ(JHITS-ISET)                                        *
C.    *       JHD=LQ(JH-IDET)                                          *
C.    *       IQ(JH+IDET)= pointer to first free word in JHD           *
C.    *         Each hit is packed into JHD in the following format    *
C.    *        --Track number ITRA not packed                          *
C.    *        --Volume numbers packed                                 *
C.    *        --Hits transformed and packed                           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT ,<GXINT>, GINC4, GKHITS       *
C.    *       Authors : R.Brun      **********                         *
C.    *                 P.Zanarini  *                                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
      PARAMETER (NHEMX=100,NVMAX=20)
      DIMENSION HITS(NHEMX),KWS(120),NUMBV(NVMAX)
      EQUIVALENCE (WS(1),NUMBV(1)),(WS(101),HITS(1)),(WS(201),KWS(1))
      CHARACTER*4 IUSET,IUDET
C.
C.    ------------------------------------------------------------------
C.
      IF(JHITS.LE.0)GO TO 999
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
      JH=LQ(JHITS-ISET)
      IF(JS.LE.0)GO TO 130
      IF(JH.LE.0)GO TO 130
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
      JHD=LQ(JH-IDET)
      IF(JHD.LE.0)GO TO 120
      JDH=LQ(JD-1)
C
      WRITE(CHMAIL,1000)IQ(JS+IDET),IQ(JSET+ISET)
      CALL GMAIL(1,1)
C
      ILAST=IQ(JH+IDET)
      IF(ILAST.EQ.0)GO TO 120
      NV=IQ(JD+2)
      NH=IQ(JD+4)
      NW=IQ(JD+1)+IQ(JD+3)+1
C
      CALL VBLANK(KWS,120)
      K=0
      IF(NV.GT.0) THEN
C        Number of printed elements limited to 9
         NVM=MIN(NV,NVMAX)
         DO 22 I=1,NVM
            CALL UBLOW(IQ(JD+2*I+9),KWS(K+1),4)
            K=K+5
  22     CONTINUE
      ENDIF
      IF(NH.GT.0)THEN
         K=K+5
         IP0 = K
         NHM=MIN(NH,NHEMX)
         DO 26 I=1,NHM
            IF(K.GT.110)GO TO 26
            NH1=I
            CALL UBLOW(IQ(JDH+4*I-3),KWS(K+1),4)
            K=K+10
  26     CONTINUE
      ENDIF
      WRITE(CHMAIL,2000)(KWS(I),I=1,K)
      CALL GMAIL(0,0)
C
      IF(NH.GT.0)THEN
         IF(NH1.LT.NHM)THEN
            CALL VBLANK(KWS,IP0)
            DO 40 NHP1=NH1+1,NHM,NH1
               NHP2=NHP1+NH1-1
               NHP2=MIN(NHP2,NHM)
               K=IP0
               DO 30 I=NHP1,NHP2
                  CALL UBLOW(IQ(JDH+4*I-3),KWS(K+1),4)
                  K=K+10
  30           CONTINUE
               WRITE(CHMAIL,2100)(KWS(I),I=1,K)
               CALL GMAIL(0,0)
  40        CONTINUE
         ENDIF
      ENDIF
C
C             Now loop on all hits to get volume numbers and hits
C
      IHIT=0
      DO 110 I=1,ILAST,NW
      IHIT=IHIT+1
C
C             IPKHIT in /CGDRAW/ is used to print just that hit
C
      IF (IPKHIT.GT.0.AND.IHIT.NE.IPKHIT) GO TO 110
C
      ITRA=IQ(JHD+I)
C
C             Get unpacked volume numbers
C
      NK=1
      IF(NV.GT.0)THEN
         K=1
         DO 50 IV=1,NVM
            NB=IQ(JD+2*IV+10)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   K=1
                   NK=NK+1
               ENDIF
               NUMBV(IV)=IQ(JHD+I+NK)
               IF(IV.NE.NV)NK=NK+1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               NUMBV(IV)=IBITS(IQ(JHD+I+NK),K-1,NB)
               K=K+NB
            ENDIF
  50  CONTINUE
      NK=NK+1
      ENDIF
C
C             Get unpacked hits
C             Hits origin is shifted . Division by scale factor
C
      IF(NH.GT.0)THEN
         K=1
         DO 90 IH=1,NHM
            NB=IQ(JDH+4*IH-2)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   NK=NK+1
               ENDIF
               KHIT=IQ(JHD+I+NK)
               K=1
               NK=NK+1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               KHIT=IBITS(IQ(JHD+I+NK),K-1,NB)
               K=K+NB
            ENDIF
            HITS(IH)=FLOAT(KHIT)/Q(JDH+4*IH) - Q(JDH+4*IH-1)
  90  CONTINUE
      ENDIF
C
      WRITE(CHMAIL,3000)IHIT,ITRA,(NUMBV(L),L=1,NV)
      IF(NH.GT.0)THEN
         IP0=NV*5+15
         DO 100 NHP1=1,NHM,NH1
            NHP2=NHP1+NH1-1
            NHP2=MIN(NHP2,NHM)
            WRITE(CHMAIL(IP0+1:),4000)(HITS(L),L=NHP1,NHP2)
            CALL GMAIL(0,0)
            CHMAIL(1:IP0)=' '
 100     CONTINUE
      ENDIF
C
      IF(IHIT.EQ.IPKHIT)GO TO 999
C
 110  CONTINUE
 120  CONTINUE
 130  CONTINUE
C
 1000 FORMAT(' =====>HITS IN DETECTOR ** ',A4,
     +' ** OF SET ** ',A4,' **')
 2000 FORMAT('  HIT TRACK ',120A1)
 2100 FORMAT(12X,120A1)
 3000 FORMAT(1X,I4,I6,20I5)
 4000 FORMAT(12(G10.3))
 999  RETURN
      END
