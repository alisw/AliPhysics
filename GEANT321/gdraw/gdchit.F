*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCHIT (IUSET,IUDET,ITRS,ISYMB,SIZMAX,KDHIT,
     +                   HITMIN,HITMAX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw Calorimetry Hits of track ITRS in detector IUDET    *
C.    *       of set IUSET, in case that the spatial coordinates       *
C.    *       of the HITS had been stored as the first three words     *
C.    *       in the argument vector of GSAHIT.                        *
C.    *                                                                *
C.    *       The following characters may be chosen :                 *
C.    *          ISYMB = -1   hardware character    (fast),            *
C.    *                =  0   crosses X by TVDRAW   (default),         *
C.    *                       (size = SSYMB in cm, default = 0.1),     *
C.    *                = 840  characters by HPLSOF                     *
C.    *                   .   according to the table p. 25             *
C.    *                   .   of the HPLOT manual                      *
C.    *                  853  (size = SSYMB in cm)                     *
C.    *             (844,854  white, black asterisk * in addition).    *
C.    *                                                                *
C.    *       SIZMAX is the maximum character size                     *
C.    *                                                                *
C.    *       KDHIT is the hit number selected : the character size    *
C.    *       is a function of HITS(KDHIT) and HITMIN and HITMAX are   *
C.    *       its limits                                               *
C.    *                                                                *
C.    *                                                                *
C.    *       JH =LQ(JHITS-ISET)                                       *
C.    *       JHD=LQ(JH-IDET)                                          *
C.    *       IQ(JH+IDET)=pointer to the LAST USED word in JHD         *
C.    *         Each hit was packed into JHD in the following format:  *
C.    *        --Track number ITRA not packed;                         *
C.    *        --Volume numbers packed;                                *
C.    *        --Hits transformed and packed.                          *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Authors : R.Brun, P.Zanarini   *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
      DIMENSION HITS(3,10),ASIZE(10)
      CHARACTER*4 IUSET,IUDET
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH('OFF ',IOFF,4,4)
C
C            Set IOBJ to HIT
C
      IOBJ=3
C
      IF (JHITS.LE.0)           GOTO 140
      NSET = IQ(JSET-1)
      NS1        = 1
      NS2        = NSET
      IF (IUSET(1:1).EQ.'*')   GOTO  10
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF (ISET .LE.0)           GOTO 140
      NS1        = ISET
      NS2        = ISET
   10 CONTINUE
C
C            loop on all selected sets
C
      DO 130   ISET=NS1,NS2
      JS = LQ(JSET-ISET)
      JH = LQ(JHITS-ISET)
      IF(JS.LE.0)GO TO 130
      IF(JH.LE.0)GO TO 130
      NDET = IQ(JS-1)
      ND1        = 1
      ND2        = NDET
      IF (IUDET(1:1).EQ.'*')   GOTO  20
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF (IDET .EQ.0)   GOTO 130
      ND1        = IDET
      ND2        = IDET
   20 CONTINUE
C
C            loop on selected detectors for this set
C
      DO 120   IDET=ND1,ND2
      JD = LQ(JS-IDET)
      JHD = LQ(JH-IDET)
      IF (JHD  .LE.0)   GOTO 120
      JDH = LQ(JD-1)
      ILAST = IQ(JH+IDET)
      IF (ILAST.EQ.0)   GOTO 120
      NV = IQ(JD+2)
      NH = IQ(JD+4)
      NW = IQ(JD+1) + IQ(JD+3) + 1
C
C            Now loop on all hits to get their space coordinates
C
      MH         = 0
      DO 110   I=1,ILAST,NW
      ITRA = IQ(JHD+I)
C
C            Is this the right track ?
C
      IF (ITRS.NE.0 .AND. ITRS.NE.ITRA) GOTO 110
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
      NK         = 1
      K          = 1
      IF (NV.LE.0)      GOTO  40
      DO 30   IV=1,NV
      NB = IQ(JD+2*IV+10)
      IF (NB.EQ.0)  NB=32
      K          = NB + K
      IF (K.LE.33)      GOTO  30
      NK         = NK + 1
      K          = NB + 1
   30 CONTINUE
      NK         = NK + 1
C
C            Get unpacked (first three !) HITS.
C            Hits origin is shifted. Division by scale factor.
C
   40 CONTINUE
      IF (NH.LE.2)      GOTO 110
      MH         = MH + 1
C
      IF (KDHIT.GT.3.AND.KDHIT.LE.10) THEN
        KDH=KDHIT
      ELSE
        KDH=3
        HITKH=0.
      ENDIF
C
      K=1
      DO 80   IH=1,KDH
      NB = IQ(JDH+4*IH-2)
      IF (NB.GT.0)      GOTO  50
      IF (K.GT.1) THEN
        NK=NK+1
      ENDIF
      KHIT = IQ(JHD+I+NK)
      K          = 1
      NK         = NK + 1
      GOTO  70
   50 IF (K+NB.LE.33)   GOTO  60
      K          = 1
      NK         = NK + 1
   60 KHIT = IBITS(IQ(JHD+I+NK),K-1,NB)
      K          = K + NB
   70 CONTINUE
      IF (IH.LE.3) THEN
        HITS(IH,MH)= FLOAT(KHIT)/Q(JDH+4*IH) - Q(JDH+4*IH-1)
      ELSE
        HITKH=FLOAT(KHIT)/Q(JDH+4*IH) - Q(JDH+4*IH-1)
      ENDIF
   80 CONTINUE
C
C             Compute symbol size
C
      IF (HITMAX.LE.0) THEN
        SSYMB=SIZMAX
      ELSE
        SSYMB=SIZMAX*(HITKH-HITMIN)/HITMAX
        IF (SSYMB.GT.SIZMAX) SSYMB=SIZMAX
C
C             HPLSOF accepts a minimum SIZE of 0.001,
C             otherwise set SIZE=0.25 : protect here this case
C
        IF (SSYMB.LT.0.0011) SSYMB=0.0011
      ENDIF
      ASIZE(MH)=SSYMB
C
C            Draw the hits (usually 10 at a time);
C            flush HITS vector if
C             --loop on this detector volume is to finish,
C             --loop on current track is to finish,
C             --HITS vector is full.
C
      IF ((ILAST-I).LT.NW)   GOTO 90
      ITRN = IQ(JHD+I+NW)
      IF (ITRA.NE.ITRN)      GOTO 90
      IF (MH  .NE.  10)      GOTO 110
C
   90 CONTINUE
C
      DO 100 II=1,MH
        CALL GDAHIT(HITS(1,II),HITS(2,II),HITS(3,II),ISYMB,ASIZE(II))
  100 CONTINUE
      MH=0
C
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
C
C             Reset ITR3D to 'standard projection'
C
      ITR3D=0
      IOBJ=0
C
      END
