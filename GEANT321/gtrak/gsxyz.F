*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.24  by  S.Giani
*-- Author :
      SUBROUTINE GSXYZ
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Store space point VECT in banks JXYZ                     *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Authors    R.Brun, F.Carena  *********                   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcpush.inc"
#include "geant321/gckine.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF(JXYZ.EQ.0)THEN
C
C               Create mother bank JXYZ
C
         CALL MZBOOK(IXDIV,JXYZ,JXYZ,1,'JXYZ',NCKINE,NCKINE,0,3,0)
         IQ(JXYZ-5)=0
      ENDIF
C
      NXYZ=IQ(JXYZ-2)
      IF(ITRA.GT.NXYZ)THEN
         NPUSH=MAX(10,ITRA-NXYZ)
         CALL MZPUSH(IXDIV,JXYZ,NPUSH,0,'I')
      ENDIF
C
C              Create bank to store points for primary track ITRA
C
      JX=LQ(JXYZ-ITRA)
      IF(JX.EQ.0)THEN
         CALL MZBOOK(IXDIV,JX,JXYZ,-ITRA,'JXYZ',0,0,NCJXYZ,3,0)
         Q(JX+1)=NCJXYZ-2.
         Q(JX+2)=3.
         IQ(JX)=NTMULT
C
C              Check if there is enough space to store current point
C
      ELSE
         NFREE=Q(JX+1)
         IF(NFREE.LE.6)THEN
            CALL MZPUSH(IXDIV,JX,0,NPJXYZ,'I')
            JX=LQ(JXYZ-ITRA)
            Q(JX+1)=Q(JX+1)+NPJXYZ
         ENDIF
C
         IF(NTMULT.NE.IQ(JX))THEN
            IQ(JX)=NTMULT
            IPOINT=JX+Q(JX+2)
            IF(Q(IPOINT).NE.0.) THEN
               Q(JX+2)=Q(JX+2)+3.*Q(IPOINT)+2.
               Q(JX+1)=Q(JX+1)-2.
            ENDIF
         ENDIF
      ENDIF
C
      IPOINT=JX+Q(JX+2)
      JP=IPOINT+3.*Q(IPOINT)+1.
      Q(IPOINT)=Q(IPOINT)+1.
      Q(IPOINT+1)=IPART
      DO 10 I=1,3
  10  Q(JP+I)=VECT(I)
      Q(JX+1)=Q(JX+1)-3.
C
      END
