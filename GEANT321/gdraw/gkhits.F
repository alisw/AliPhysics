*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GKHITS (IUSET,IUDET,EPSILO)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints the hit picked (using the graphics cursor)        *
C.    *       from all the hits of all tracks in detector IUDET        *
C.    *       of set IUSET, in case that the spatial coordinates       *
C.    *       of the hits had been stored as the first three words     *
C.    *       in the argument vector of GSAHIT.                        *
C.    *       If the character typed to pick is 'K' or 'k' then the    *
C.    *       kinematics of the corresponding track is also printed.   *
C.    *       The picking is repeated until the character typed is     *
C.    *       'Q' or 'q'.                                              *
C.    *       EPSILO is the pick aperture, i.e. the first hit that     *
C.    *       falls into the square centered at U0,V0 (the pick point) *
C.    *       with half length EPSILO is kept.                         *
C.    *       If EPSILO is negative then its absolute value is taken   *
C.    *       and in addition the pick square is drawn.                *
C.    *       If EPSILO is zero there is no pick aperture and over all *
C.    *       the hits the one nearest to the pick point is taken.     *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Authors : R.Brun, W.Gebel, P.Zanarini   *********        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
      DIMENSION HITS(3,10)
      DIMENSION U(10),V(10)
      DIMENSION UP(5),VP(5)
      CHARACTER*4 IUSET,IUDET,KIUSET,KIUDET
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH('OFF ',IOFF,4,4)
C
C            Set IOBJ to HIT
C
      IOBJ=3
C
   10 CONTINUE
C
      CALL GDCURS(U0,V0,KCHAR)
      IF (KCHAR.EQ.0) GO TO 160
C
      IF (JHITS.LE.0)           GOTO 160
      NSET = IQ(JSET-1)
      NS1        = 1
      NS2        = NSET
      IF (IUSET(1:1).EQ.'*')   GOTO  20
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF (ISET .LE.0)           GOTO 160
      NS1        = ISET
      NS2        = ISET
C
   20 CONTINUE
C
      EPSIL=ABS(EPSILO)
      EPSU=PLTRNX*EPSIL/20.
      EPSV=PLTRNY*EPSIL/20.
      IHITNR=0
      DISTPK=1000000.
C
C             Draw the pick box
C
      IF (EPSILO.LT.0) THEN
         UP(1)=U0-EPSU
         VP(1)=V0-EPSV
         UP(2)=U0+EPSU
         VP(2)=VP(1)
         UP(3)=UP(2)
         VP(3)=V0+EPSV
         UP(4)=UP(1)
         VP(4)=VP(3)
         UP(5)=UP(1)
         VP(5)=VP(1)
         CALL GDRAWV(UP,VP,5)
      ENDIF
C
C            loop on all selected sets
C
      DO 140   ISET=NS1,NS2
C
         JS = LQ(JSET-ISET)
         JH = LQ(JHITS-ISET)
         IF(JS.LE.0)GO TO 140
         IF(JH.LE.0)GO TO 140
         NDET = IQ(JS-1)
         ND1 = 1
         ND2 = NDET
         IF (IUDET(1:1).EQ.'*') GOTO 30
         CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
         IF (IDET .EQ.0) GOTO 140
         ND1 = IDET
         ND2 = IDET
   30    CONTINUE
C
C            loop on selected detectors for this set
C
         DO 130 IDET=ND1,ND2
C
            JD = LQ(JS-IDET)
            JHD = LQ(JH-IDET)
            IF (JHD .LE.0) GOTO 130
            JDH = LQ(JD-1)
            ILAST = IQ(JH+IDET)
            IF (ILAST.EQ.0) GOTO 130
            NV = IQ(JD+2)
            NH = IQ(JD+4)
            NW = IQ(JD+1) + IQ(JD+3) + 1
C
C            Now loop on all hits to get their space coordinates
C
            MH = 0
            IHIT = 0
            DO 120 I=1,ILAST,NW
               IHIT = IHIT + 1
               ITRA = IQ(JHD+I)
C
C            Is THRZ option set on (R-Z projection)
C            and the cut lateral (vertical or horizontal) ?
C
               IF (ITHRZ.NE.IOFF) THEN
                  IF (ICUT.EQ.1.OR.ICUT.EQ.2) ITR3D=-ITRA
               ELSE
                  ITR3D=0
               ENDIF
C
C            Skip words with volume numbers
C
               NK = 1
               K = 1
               IF (NV.LE.0) GOTO 50
               DO 40 IV=1,NV
                  NB = IQ(JD+2*IV+10)
                  IF (NB.EQ.0) NB=32
                  K = NB + K
                  IF (K.LE.33) GOTO 40
                  NK = NK + 1
                  K = NB + 1
   40          CONTINUE
               NK = NK + 1
C
C            Get unpacked (first three !) HITS.
C            Hits origin is shifted. Division by scale factor.
C
   50          CONTINUE
               IF (NH.LE.2) GOTO 120
               MH = MH + 1
               K=1
               DO 90 IH=1,3
                  NB = IQ(JDH+4*IH-2)
                  IF (NB.GT.0) GOTO 60
                  IF (K.GT.1) THEN
                     NK=NK+1
                  ENDIF
                  KHIT = IQ(JHD+I+NK)
                  K = 1
                  NK = NK + 1
                  GOTO 80
   60             IF (K+NB.LE.33) GOTO 70
                  K = 1
                  NK = NK + 1
   70             KHIT = IBITS(IQ(JHD+I+NK),K-1,NB)
                  K = K + NB
   80             CONTINUE
                  HITS(IH,MH)= FLOAT(KHIT)/Q(JDH+4*IH) - Q(JDH+4*IH-1)
   90          CONTINUE
C
C            Scan the hits (usually 10 at a time);
C            flush HITS vector if
C             --loop on this detector volume is to finish,
C             --loop on current track is to finish,
C             --HITS vector is full.
C
               IF ((ILAST-I).LT.NW) GOTO 100
               ITRN = IQ(JHD+I+NW)
               IF (ITRA.NE.ITRN) GOTO 100
               IF (MH .NE. 10) GOTO 120
C
  100          CONTINUE
C
               CALL GDFR3D(HITS(1,1),MH,U,V)
C
               DO 110 IP=1,MH
C
                  IF (EPSILO.NE.0) THEN
                     IF (U(IP).LT.U0-EPSU) GO TO 110
                     IF (U(IP).GT.U0+EPSU) GO TO 110
                     IF (V(IP).LT.V0-EPSV) GO TO 110
                     IF (V(IP).GT.V0+EPSV) GO TO 110
                     II=IP
                     IHITNR=IHIT-MH+II
                     JS = LQ(JSET-ISET)
                     CALL UHTOC(IQ(JSET+ISET),4,KIUSET,4)
                     CALL UHTOC(IQ(JS+IDET)  ,4,KIUDET,4)
                     KITRA=ITRA
                     GO TO 150
                  ELSE
                     DIST=(U(IP)-U0)*(U(IP)-U0)+(V(IP)-V0)*(V(IP)-V0)
                     IF (DIST.LT.DISTPK) THEN
                        DISTPK=DIST
                        II=IP
                        IHITNR=IHIT-MH+II
                        JS = LQ(JSET-ISET)
                        CALL UHTOC(IQ(JSET+ISET),4,KIUSET,4)
                        CALL UHTOC(IQ(JS+IDET)  ,4,KIUDET,4)
                        KITRA=ITRA
                     ENDIF
                  ENDIF
C
  110          CONTINUE
C
               MH=0
C
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
      IF (EPSILO.EQ.0) GO TO 150
C
C             Hit not picked
C
      WRITE (CHMAIL,1000)
      CALL GMAIL(0,0)
C
C             Restart for another pick operation
C
      GO TO 10
C
C             Hit picked
C
  150 CONTINUE
C
C             Highlight the picked hit
C
      UDELTA=0.05*PLTRNX/20.
      VDELTA=0.05*PLTRNY/20.
      UP(1)=U(II)-UDELTA
      VP(1)=V(II)-VDELTA
      UP(2)=U(II)+UDELTA
      VP(2)=VP(1)
      UP(3)=UP(2)
      VP(3)=V(II)+VDELTA
      UP(4)=UP(1)
      VP(4)=VP(3)
      UP(5)=UP(1)
      VP(5)=VP(1)
      CALL GDRAWV(UP,VP,5)
C
C             Print parameters of the picked hit
C
      IPKHIT=IHITNR
      IF (IPKHIT.GT.0) CALL GPHITS(KIUSET,KIUDET)
      IPKHIT=0
      CALL GPKINE(KITRA)
C
C             Restart for another pick operation
C
      GO TO 10
C
  160 CONTINUE
C
C             Reset ITR3D to 'standard projection'
C
      ITR3D=0
      IOBJ=0
C
 1000 FORMAT (' GKHITS: No hit picked')
      RETURN
      END
