*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GFHITS(IUSET,IUDET,NVDIM,NHDIM,NHMAX,ITRS,NUMVS
     +,                 ITRA,NUMBV,HITS,NHITS)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *   Returns the  hits produced  by track ITRS  (if 0,   by all   *
C.    * tracks)  in the physical volume  specified by the list NUMVS   *
C.    * with generic volume name IUDET belonging to set IUSET.         *
C.    * IUSET     user set identifier                                  *
C.    * IUDET     user    detector   identifier    (name   of    the   *
C.    *           corresponding sensitive volume)                      *
C.    * NVDIM     1st dimension of NUMBV and NUMVS (usually =NV, the   *
C.    *           number of volume descriptors which permit to identify*
C.    *           a given detector, eventually smaller than  NV)       *
C.    * NHDIM     1st  dimension  of  array  HITS  (argument  NH  of   *
C.    *           GSDETH)                                              *
C.    * NHMAX     maximum number of hits to be returned                *
C.    * ITRS      number of  the selected  track.   If  ITRS=0,  all   *
C.    *           tracks are taken                                     *
C.    * NUMVS     is a  1-Dim array that  must contain on  input the   *
C.    *           geometric  path  of  the  detector  volume  to  be   *
C.    *           selected.    All 0  interpreted  as 'all  physical   *
C.    *           volumes with generic names IUDET'                    *
C.    * ITRA      is a 1-Dim  array that will contain  on output for   *
C.    *           each  hit  the  number  of  the  track  which  has   *
C.    *           produced it                                          *
C.    * NUMBV     2-Dim array that  will contain on output  for each   *
C.    *           hit the list of volume numbers which identify each   *
C.    *           physical volume.   Zeroed when no more volumes are   *
C.    *           stored                                               *
C.    * HITS      2-Dim array that will contain the NHITS hits         *
C.    * NHITS     returns the number of selected  hits.  In case the   *
C.    *           total number of hits is greater than NHMAX,  NHITS   *
C.    *           is set to NHMAX+1 and only NHMAX hits are returned   *
C.    *   - HITS(1,I)  is the element 1 for hit number I               *
C.    *   - NUMBV(1,I) is the volume number 1 for hit number I         *
C.    *   - ITRA(I)    is  the  track number  corresponding  to  hit   *
C.    * number I                                                       *
C.    *   In the calling routine the arrays NUMVS,  NUMBV,  HITS and   *
C.    * ITRA must be dimensioned to:                                   *
C.    *     NUMVS(NVDIM)                                               *
C.    *     NUMBV(NVDIM,NHMAX)                                         *
C.    *     HITS(NHDIM,NHMAX)                                          *
C.    *     ITRA(NHMAX)                                                *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUDIGI                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      PARAMETER (NVMAX=20)
      DIMENSION NUMVT(NVMAX),NUMVS(NVDIM),NUMBV(NVDIM,1),ITRA(1)
      DIMENSION HITS(NHDIM,1)
      EQUIVALENCE (WS(1),NUMVT(1))
      CHARACTER*4 IUSET,IUDET
C.
C.    ------------------------------------------------------------------
C.
C             Find if selected set, detector exists
C
      NHITS=0
      IF(JHITS.LE.0)GO TO 999
      NSET=IQ(JSET-1)
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF(ISET.LE.0)GO TO 999
C
      JS=LQ(JSET-ISET)
      JH=LQ(JHITS-ISET)
      IF(JS.LE.0)GO TO 999
      IF(JH.LE.0)GO TO 999
      NDET=IQ(JS-1)
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF(IDET.EQ.0)GO TO 999
C
      JD=LQ(JS-IDET)
      JHD=LQ(JH-IDET)
      IF(JHD.LE.0)GO TO 999
      JDH=LQ(JD-1)
C
      ILAST=IQ(JH+IDET)
      IF(ILAST.EQ.0)GO TO 999
      NV=IQ(JD+2)
      NH=IQ(JD+4)
      NW=IQ(JD+1)+IQ(JD+3)+1
C
C
C             Loop on all hits
C
C
      IHIT=0
      DO 100 I=1,ILAST,NW
C
C             Find the selected track
C
      ITRT=IQ(JHD+I)
      IF(ITRS.NE.0 .AND. ITRS.NE.ITRT)GO TO 100
C
C             Find the selected volume
C             (if NO volumes exist take ALL hits)
C
      NK=1
      IF(NV.GT.0)THEN
         K=1
         DO 40 IV=1,NV
            NB=IQ(JD+2*IV+10)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   K=1
                   NK=NK+1
               ENDIF
               IF(IV.LE.NVMAX)NUMVT(IV)=IQ(JHD+I+NK)
               IF(IV.NE.NV)NK=NK+1
            ELSE
               IF(K+NB.GT.33)THEN
                 K=1
                 NK=NK+1
               ENDIF
               IF(IV.LE.NVMAX)NUMVT(IV)=IBITS(IQ(JHD+I+NK),K-1,NB)
               K=K+NB
            ENDIF
            IF(IV.LE.NVDIM)THEN
               IF(NUMVS(IV).NE.0)THEN
                  IF(NUMVS(IV).NE.NUMVT(IV))GO TO 100
               ENDIF
            ENDIF
   40    CONTINUE
         NK=NK+1
      ENDIF
C
C
C
C ========>   Now store track number and volume numbers and fetch hits
C
      IHIT=IHIT+1
      IF(IHIT.GT.NHMAX)GO TO 110
C
      ITRA(IHIT)=ITRT
      NVMIN=MIN(NV,NVDIM)
      DO 55 J=1,NVDIM
  55  NUMBV(J,IHIT)=0
      DO 57 J=1,NVMIN
  57  NUMBV(J,IHIT)=NUMVT(J)
C
C             Get unpacked hits
C             Hits origin is shifted . Division by scale factor
C
      IF(NH.GT.0)THEN
      K=1
      DO 90 IH=1,NH
         NB=IQ(JDH+4*IH-2)
         IF(NB.LE.0)THEN
            IF(K.GT.1)THEN
                K=1
                NK=NK+1
            ENDIF
            KHIT=IQ(JHD+I+NK)
            NK=NK+1
         ELSE
            IF(K+NB.GT.33)THEN
               K=1
               NK=NK+1
            ENDIF
            KHIT=IBITS(IQ(JHD+I+NK),K-1,NB)
            K=K+NB
         ENDIF
         IF(IH.LE.NHDIM)THEN
            HITS(IH,IHIT)=FLOAT(KHIT)/Q(JDH+4*IH) - Q(JDH+4*IH-1)
         ENDIF
  90  CONTINUE
      ENDIF
 100  CONTINUE
C
 110  NHITS=IHIT
C
 999  RETURN
      END
