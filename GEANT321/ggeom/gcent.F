*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GCENT(IAX,X0)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       SUBROUTINE TO CENTRE A SLICE OF THE NLEVEL TH            *
C     *       SHAPE, WHEN THE CENTRE OF THE SLICE ALONG AXIS           *
C     *       NUMBER IAX AFTER ROTATION BY MATRIX IROT IS AT           *
C     *       X0(IAX). THE CENTER IN THE UNROTATED FRAME IS            *
C     *       RETURNED.                                                *
C     *                                                                *
C.    *    ==>Called by : GFTRAC, GINVOL, GLVOLU, GMEDIV, GNEXT,       *
C.    *                   GTMEDI, GTNEXT                               *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
      DIMENSION X0(3),PAR(100),ATT(20)
C.
C.    ------------------------------------------------------------------
C.
C
      IVO=LVOLUM(NLEVEL)
      JVO=LQ(JVOLUM-IVO)
      ISH=Q(JVO+2)
C
      IF(ISH.EQ.4) GO TO 20
C
C          NOT A GENERAL TRAPEZOID.
C
      IF(ISH.EQ.10.AND.IAX.NE.1) GO TO 20
C
C         NOT A PARALLELEPIPED WITH NON X AXIS SLICING.
C
      DO 10 I=1,3
      IF(I.NE.IAX) X0(I)=0.0
   10 CONTINUE
C
      GO TO 999
   20 CONTINUE
C
C
      JIN=0
      IF(NLEVEL.LE.1) GO TO 30
C
      IN   = LINDEX(NLEVEL)
      IVOO = LVOLUM(NLEVEL-1)
      JVOO = LQ(JVOLUM-IVOO)
      IF (Q(JVOO+3).LT.0.) THEN
         JIN = LQ(JVOO-1)
      ELSE
         JIN = LQ(JVOO-IN)
      ENDIF
C
   30 CONTINUE
C
C
      IF (NLEVEL.LT.NLDEV(NLEVEL)) THEN
         JPAR = 0
      ELSE
*       (case with structure JVOLUM locally developed)
         JPAR = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         IF (NLEVEL.EQ.NLDEV(NLEVEL)) GO TO 120
         DO 110 ILEV = NLDEV(NLEVEL), NLEVEL-1
            IF (IQ(JPAR+1).EQ.0) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
               IF (JPAR.EQ.0) GO TO 120
            ELSE IF (IQ(JPAR-3).GT.1) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
            ELSE
               JPAR = LQ(JPAR-1)
            ENDIF
            IF (ILEV.EQ.NLEVEL-1) THEN
               JPAR = JPAR + 5
               NPAR = IQ(JPAR)
               CALL UCOPY (Q(JPAR+1), PAR, NPAR)
            ENDIF
  110    CONTINUE
         GO TO 130
      ENDIF
*      (normal case)
  120 CONTINUE
*
C
      CALL GFIPAR(JVO,JIN,IN,NPAR,NATT,PAR,ATT)
C
  130 CONTINUE
      IF(ISH.EQ.10) GO TO 40
C
C          A GENERAL TRAPEZOID
C
      IF(IAX.NE.3) GO TO 41
C
      X0(1)=PAR(2)*X0(IAX)
      X0(2)=PAR(3)*X0(IAX)
C
      GO TO 999
C
   41 IF(IAX.NE.2) GO TO 900
C
      X0(1) = X0(IAX)*PAR(7)
      GO TO 999
C
   40 CONTINUE
C
C          PARALLELEPIPED WITH NON X AXIS SLICING.
C
      X0(1)=PAR(IAX+2)*X0(IAX)
C
      IF(IAX.EQ.3) GO TO 50
C
      X0(3)=0.0
C
      GO TO 999
   50 CONTINUE
C
      X0(2)=PAR(6)*X0(3)
C
      GO TO 999
C
  900 CONTINUE
C
C          GENERAL TRAPEZOID SLICED OTHER THAN IN Z
C          DOES NOT MAKE SENSE.
C
      WRITE(CHMAIL,1000) ISH,IAX
      CALL GMAIL(0,0)
      WRITE(CHMAIL,2000)
      CALL GMAIL(0,0)
 1000 FORMAT(' ISH,IAX =',2I5,' DOES NOT MAKE SENSE.')
 2000 FORMAT(' IT IS A GENERAL TRAPEZOID SLICED OTHER THAN IN Z')
C
C
  999 CONTINUE
      END
