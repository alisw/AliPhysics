*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 02/07/94  18.24.47  by  S.Giani
*-- Author :
      SUBROUTINE GSROTM(NMAT,THETA1,PHI1,THETA2,PHI2,THETA3,PHI3)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       STORE ROTATION MATRICES                                  *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *         Author  R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcnum.inc"
      DIMENSION ANGLES(6),IP(3)
      SAVE SINMIN
#if defined(CERNLIB_SINGLE)
      DIMENSION ROTMAT(9)
      PARAMETER (ONE=1.0,ZERO=0.0)
      DATA SINMIN/1.00E-5/
#endif
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION ROTMAT(9),ONE,ZERO
      DOUBLE PRECISION PROD1,PROD2,PROD3,HMOD,SINMIN
      PARAMETER (ONE=1.D0,ZERO=0.D0)
      DATA SINMIN/1.00D-5/
#endif
C.
C.    ------------------------------------------------------------------
C.
      IF(NMAT.LE.0)GO TO 999
      IF(JROTM.LE.0)CALL MZBOOK(IXCONS,JROTM,JROTM,1,'ROTM',NROTM,NROTM,
     +0,3,0)
      IF(NMAT.GT.IQ(JROTM-2)) THEN
         NPUSH=NMAT-IQ(JROTM-2)+50
         CALL MZPUSH(IXCONS,JROTM,NPUSH,0,'I')
         NROTM=IQ(JROTM-2)
         JR1=0
      ELSE
         JR1=LQ(JROTM-NMAT)
         IF(JR1.GT.0)THEN
            WRITE(CHMAIL,1000)
            CALL GMAIL(1,0)
            CALL GPROTM(NMAT)
            CALL MZDROP(IXCONS,LQ(JROTM-NMAT),' ')
         ENDIF
      ENDIF
      CALL MZBOOK(IXCONS,JR,JROTM,-NMAT,'ROTM',0,0,16,3,0)
C
      Q(JR + 11) = THETA1
      Q(JR + 12) = PHI1
      Q(JR + 13) = THETA2
      Q(JR + 14) = PHI2
      Q(JR + 15) = THETA3
      Q(JR + 16) = PHI3
C
      DO  10 N = 1,3
        THERAD = Q(JR+ 9+2*N)*DEGRAD
        PHIRAD = Q(JR+10+2*N)*DEGRAD
        SINTHE = SIN(THERAD)
        Q(JR+3*N-2) = SINTHE * COS(PHIRAD)
        Q(JR+3*N-1) = SINTHE * SIN(PHIRAD)
        Q(JR+3*N  ) = COS(THERAD)
        CALL VUNIT (Q(JR+3*N-2),Q(JR+3*N-2),3)
  10  CONTINUE
C.
C.---       Test orthonormality
      DO 20 J=1,9
        ROTMAT(J)=Q(JR+J)
  20  CONTINUE
      PROD2=ZERO
C.
C.               X - Y
      PROD1=
     +ROTMAT(1)*ROTMAT(4)+ROTMAT(2)*ROTMAT(5)+ROTMAT(3)*ROTMAT(6)
      IF(ABS(PROD1).GT.SINMIN) GO TO 30
C.
C.               X - Z
      PROD2=
     +ROTMAT(1)*ROTMAT(7)+ROTMAT(2)*ROTMAT(8)+ROTMAT(3)*ROTMAT(9)
      IF(ABS(PROD2).GT.SINMIN) GO TO 30
C.
C.               Y - Z
      PROD3=
     +ROTMAT(7)*ROTMAT(4)+ROTMAT(8)*ROTMAT(5)+ROTMAT(9)*ROTMAT(6)
      IF(ABS(PROD3).LE.SINMIN) GO TO 110
  30  CONTINUE
C.
C.---       Orthonormalization needed
C.
C.          Assume X correct
      HMOD=ZERO
      DO 40 J=4,6
        ROTMAT(J)=ROTMAT(J)-ROTMAT(J-3)*PROD1
        HMOD=HMOD+ROTMAT(J)*ROTMAT(J)
  40  CONTINUE
      HMOD=ONE/SQRT(HMOD)
      DO 50 J=4,6
        ROTMAT(J)=ROTMAT(J)*HMOD
  50  CONTINUE
C.
C.          Y done, do Z
C.
*      IF(PROD2.EQ.ZERO) THEN
*        PROD2=
*     +  ROTMAT(1)*ROTMAT(7)+ROTMAT(2)*ROTMAT(8)+ROTMAT(3)*ROTMAT(9)
*      ENDIF
*      PROD3=
*     +ROTMAT(4)*ROTMAT(7)+ROTMAT(5)*ROTMAT(8)+ROTMAT(6)*ROTMAT(9)
*      HMOD = ZERO
*      DO 60 J=1,3
*        ROTMAT(J+6)=ROTMAT(J+6)-ROTMAT(J+3)*PROD3-ROTMAT(J)*PROD2
*        HMOD = HMOD+ROTMAT(J)*ROTMAT(J)
*  60  CONTINUE
* == AV ==
      ROTMAT(7) = ROTMAT(2)*ROTMAT(6) - ROTMAT(3)*ROTMAT(5)
      ROTMAT(8) = ROTMAT(3)*ROTMAT(4) - ROTMAT(1)*ROTMAT(6)
      ROTMAT(9) = ROTMAT(1)*ROTMAT(5) - ROTMAT(2)*ROTMAT(4)
      HMOD = ZERO + ROTMAT(7)*ROTMAT(7) + ROTMAT(8)*ROTMAT(8)
     &            + ROTMAT(9)*ROTMAT(9)
* == AV ==
      HMOD = ONE/SQRT(HMOD)
      DO 70 J=7,9
        ROTMAT(J) = ROTMAT(J)*HMOD
  70  CONTINUE
C.
C.          Put back the matrix in place
      DO 80 J=1,9
        Q(JR+J)=ROTMAT(J)
  80  CONTINUE
C.
C.          Now recompute the angles
      DO 90 J=1,3
        ANGLES(J*2-1) = ACOS(MAX(-ONE,MIN(ONE,ROTMAT(J*3))))*RADDEG
        ANGLES(J*2)   = ZERO
        IF(ROTMAT(J*3-1).NE.ZERO) THEN
          ANGLES(J*2) = ATAN2(ROTMAT(J*3-1),ROTMAT(J*3-2))*RADDEG
          IF(ANGLES(2*J).LT.0.0) ANGLES(2*J) = ANGLES(2*J)+360.0
        ENDIF
  90  CONTINUE
      WRITE(CHMAIL,2000) NMAT
      CALL GMAIL(1,2)
      WRITE(CHMAIL,2001) (Q(JR+10+J),J=1,6)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,2002) ANGLES
      CALL GMAIL(0,1)
C.
C.          Put back the angles in place
      DO 100 J=1,6
        Q(JR+10+J) = ANGLES(J)
 100  CONTINUE
C.
C.---       Orthonormalization ended
 110  CONTINUE
C
      DO 130 J = 1,3
        IP(J)  = 3
        JJR=JR+J*3-3
C
        DO 120 I = 1,3
          IF(ABS(Q(JJR+I)).LT.0.99999) GO TO 120
C
          IP(J)  = I + 3
          IF(Q(JJR+I).GE.0.) GO TO 130
C
          IP(J)  = 3 - I
          GO TO 130
C
 120    CONTINUE
 130  CONTINUE
C
      Q(JR + 10) = IP(1) + 10* IP(2) + 100* IP(3)
C
      IF(JR1.GT.0) THEN
         CALL GPROTM(-NMAT)
      ENDIF
C
1000  FORMAT(' *** GSROTM ***: Warning, rotation matrix redefinition:')
2000  FORMAT(' *** GSROTM ***: ',
     +       'Parameters of matrix no. ',I4,' changed:')
2001  FORMAT(' Old values: ',6(F14.7,3X))
2002  FORMAT(' New values: ',6(F14.7,3X))
 999  RETURN
      END
