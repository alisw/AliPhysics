*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:20  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GDCGCL(ISHAPE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This Subroutine allows the clipping of a CG object       *
C.    *       built with the Hidden Line Removal by means of any       *
C.    *       kind of shape (moreover it's possible to clip the        *
C.    *       same object more than once and by different shapes)      *
C.    *       defined by 'MCVOL' Command.                              *
C.    *                                                                *
C.    *    ==>Called by :  GDCGHI                                      *
C.    *                                                                *
C.    *       Authors :  J.Salt ; S.Giani     *********                *
C.    *                                                                *
C.    ******************************************************************
C
#include "geant321/gcbank.inc"
#include "geant321/gcgobj.inc"
#include "geant321/gchiln.inc"
#include "geant321/gcspee.inc"
#include "geant321/gcmutr.inc"
*
      DIMENSION VMIN(3),VMAX(3)
*
*
*     Volume substraction. The algorithm is the following :
*
*     Check if the Clipping volume is inside Volume 'I' (First Check)
*
*      a) If Yes , Volume 'I' is Seen (IVFUN=1)
*      b) If Not , Check the following 3 cases (Second Check):
*
*          1) C. Vol. intersects volume 'I', but the volume does not include it
*              (IVFUN=2).
*          2) Volume 'I' is inside C. Vol., Then Volume is Unseen (IVFUN=0)
*          3) Volume 'I' is outside C. Vol., Then Volume is Seen (IVFUN=1)
*
*SG
      IA=JCGOBJ+1
      IB=JCGOBJ+8000
*SG
      CALL CGMNMX(Q(IA),VMIN,VMAX)
*
*      First Check
*
***SG
*
      IF(NAIN.EQ.1)THEN
        ISUB=JCGOBJ+20000
        IF(ISHAPE.EQ.29.OR.ISHAPE.EQ.99)THEN
          CALL CGCOPY(Q(IB),8000,Q(ISUB))
        ELSE
          CALL CGCOPY(Q(IA),8000,Q(ISUB))
        ENDIF
        RETURN
      ENDIF
*
*      Do it for all the volumes cutting 'I'
*
      DO 11 IJ=1,JJJ
      DO 10 K=1,3
         IF(VMIN(K).LE.BMIN(K+3*IJ-3).AND.BMIN(K+3*IJ-3).LE.VMAX(K).AND.
     +VMIN(K).LE.BMAX(K+3*IJ-3).AND.BMAX(K+3*IJ-3).LE.VMAX(K))THEN
***SG
            GOTO 10
         ELSE
            GOTO 20
         ENDIF
   10 CONTINUE
      IVFUN=1
      GOTO 50
*
*      C. Vol. is not inside 'I'  Volume. Second Check:
***SG
*
   20 IKON=0
      IDISJ=0
      DO 30 J=1,3
******         IDISJ=0
*      Do it for all the volumes cutting 'I'
*
         IF(BMIN(J+3*IJ-3).LE.VMIN(J).AND.VMIN(J).LE.BMAX(J+3*IJ-3))THEN
            IKON=IKON+1
         ELSE
            IDISJ=IDISJ+1
         ENDIF
         IF(BMIN(J+3*IJ-3).LE.VMAX(J).AND.VMAX(J).LE.BMAX(J+3*IJ-3))THEN
***SG
            IKON=IKON+1
         ELSE
            IDISJ=IDISJ+1
         ENDIF
         IF(IDISJ.EQ.6)GOTO 40
   30 CONTINUE
   40 IF(IDISJ.EQ.6)THEN
         IVFUN=1
      ELSE
         IF(IKON.EQ.6)THEN
          IF(ICUBE.EQ.JJJ)THEN
**            IVFUN=0
             IVFUN=2
          ELSE
            IVFUN=2
          ENDIF
         ELSE
            IVFUN=2
         ENDIF
      ENDIF
   50 CONTINUE
**      IF(IVFUN.EQ.0)GOTO 11
*
*****SG
*
*    Multiple clipping: you can clip, as a sequence, the same
*    volume by two different shapes
*
      IF(JJJ.EQ.2)THEN
        ISUB1=JCGOBJ+12000
        ISUB =JCGOBJ+20000
        IF(IJ.EQ.1)THEN
         IF(ISHAPE.EQ.29.OR.ISHAPE.EQ.99)THEN
           IF(IVFUN.EQ.2) CALL CGSUB(Q(IB),Q(ICLIP1),8000,Q(ISUB1))
           IF(IVFUN.EQ.1) CALL CGCOPY(Q(IB),8000,Q(ISUB1))
         ELSE
           IF(IVFUN.EQ.2) CALL CGSUB(Q(IA),Q(ICLIP1),8000,Q(ISUB1))
           IF(IVFUN.EQ.1) CALL CGCOPY(Q(IA),8000,Q(ISUB1))
         ENDIF
        ENDIF
        IF(IJ.EQ.2)THEN
          IF(IVFUN.EQ.2) CALL CGSUB(Q(ISUB1),Q(ICLIP2),8000,Q(ISUB))
          IF(IVFUN.EQ.1) CALL CGCOPY(Q(ISUB1),8000,Q(ISUB))
        ENDIF
      ELSE
        ISUB=JCGOBJ+20000
        IF(ISHAPE.EQ.29.OR.ISHAPE.EQ.99)THEN
          IF(IVFUN.EQ.2) CALL CGSUB(Q(IB),Q(ICLIP1),8000,Q(ISUB))
          IF(IVFUN.EQ.1) CALL CGCOPY(Q(IB),8000,Q(ISUB))
        ELSE
          IF(IVFUN.EQ.2) CALL CGSUB(Q(IA),Q(ICLIP1),8000,Q(ISUB))
          IF(IVFUN.EQ.1) CALL CGCOPY(Q(IA),8000,Q(ISUB))
        ENDIF
      ENDIF
   11 CONTINUE
*
*****SG
*
      END
