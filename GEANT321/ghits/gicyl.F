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
      SUBROUTINE GICYL(R,X1,X2,S1,S2,IC,XINT,SINT,PZINT,IFLAG)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Intersection of a Track with a Cylinder or a Plane         *
C.    *     --------------------------------------------------         *
C.    *                                                                *
C.    *   Calculates intersection of track (x1,x2)  with cylindrical   *
C.    * detector of radius R.  The track  is approximated by a cubic   *
C.    * in the track length.  To  improve stability,  the coordinate   *
C.    * system is shifted.                                             *
C.    * R         radius of cylinder in cm                             *
C.    * X1        x,y,z,xp,yp,zp of 1st point                          *
C.    * X2        x,y,z,xp,yp,zp of 2nd point                          *
C.    * S1(2)     S at 1st (2nd) point                                 *
C.    * IC        =1 straight line defined by x+xp                     *
C.    * IC        =2 straight line defined by x1+x2                    *
C.    * IC        =3 cubic model                                       *
C.    * XINT      x,y,z,xp,yp,zp at intersection point                 *
C.    * SINT      S at intersection point                              *
C.    * PZINT     phi,z,dphi/dr,dz/dr                                  *
C.    * IFLAG     =1 if track intersects cylinder, =0 if not           *
C.    *   Calculates  intersection  of  track  (x1,x2)   with  plane   *
C.    * parallel to (X-Z).   The track is approximated by a cubic in   *
C.    * the  track length.   To improve  stability,  the  coordinate   *
C.    * system is shifted.                                             *
C.    * YC        Y coordinate of plane                                *
C.    * X1,...    as for GICYL                                         *
C.    * IFLAG     =1 if track intersects plane,                        *
C.    *           =0 if not                                            *
C.    * Warning:  the default accuracy is  10 microns.  The value of   *
C.    *         EPSI  (internal variable)   must  be  changed for  a   *
C.    *         better precision.                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUDIGI                               *
C.    *                                                                *
C.    *       AUTHORS:R.BRUN/JJ.DUMONT FROM AN ORIGINAL ROUTINE BY     *
C.    *       H. BOERNER  KEK  OCTOBER 1982                            *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X1(6),X2(6),XINT(6),PZINT(4),A(4),B(4),C(4)
C
      DATA MAXHIT/100/
      DATA EPSI2/0.000001/
C.
C.
C.    ------------------------------------------------------------------
C.
C.
      IFLAG  = 1
      R12    = X1(1) * X1(1) + X1(2) * X1(2)
      R22    = X2(1) * X2(1) + X2(2) * X2(2)
      R2     = R * R
      DR2  = R22-R2
C
C             TRACK CROSSING THE CYLINDER FROM INSIDE OR OUTSIDE ?
C
      IF (R22.LT.R12)                            GO TO 5
      IF (R2.LT.R12)                             GO TO 90
      IF (R2.GT.R22)                             GO TO 90
      DRCTN  = 1.
      IF (IC.EQ.3) GO TO 7
C
      IF(IC.EQ.2) GOTO 30
      S=S1
      DXDS=X1(4)
      DYDS=X1(5)
      DZDS=X1(6)
      BX=X1(1)-DXDS*S
      BY=X1(2)-DYDS*S
      BZ=X1(3)-DZDS*S
      GO TO 40
C
   5  IF (R2.LT.R22)                             GO TO 90
      IF (R2.GT.R12)                             GO TO 90
      DRCTN  = - 1.
C
      IF(IC.EQ.3) GOTO 7
      IF(IC.EQ.2) GOTO 30
      S=S2
      DXDS=X2(4)
      DYDS=X2(5)
      DZDS=X2(6)
      BX=X2(1)-DXDS*S
      BY=X2(2)-DYDS*S
      BZ=X2(3)-DZDS*S
      GOTO 40
C
   30 DX=X2(1)-X1(1)
      DY=X2(2)-X1(2)
      DZ=X2(3)-X1(3)
      DS=SQRT(DX*DX+DY*DY+DZ*DZ)
      S=S1
      DXDS=DX/DS
      DYDS=DY/DS
      DZDS=DZ/DS
      BX=X1(1)-DXDS*S
      BY=X1(2)-DYDS*S
      BZ=X1(3)-DZDS*S
C
   40 AE=DYDS*DYDS+DXDS*DXDS
      IF(AE.EQ.0.) GO TO 30
      BE=DXDS*BX+DYDS*BY
      CE=BY*BY+BX*BX-R2
      SG=SIGN(1.,DR2)
      XX=BE*BE-AE*CE
      IF(XX.LE.0.) GOTO 30
      TRLEN=(SG*SQRT(ABS(XX))-BE)/AE
      XINT(1)=DXDS*TRLEN+BX
      XINT(2)=DYDS*TRLEN+BY
      XINT(3)=DZDS*TRLEN+BZ
      SINT=TRLEN
      XINT(4)=DXDS
      XINT(5)=DYDS
      XINT(6)=DZDS
      GO TO 200
C
C             SHIFT COORDINATE SYSTEM SUCH THAT CENTER OF GRAVITY=0
C
   7  IF(R.LE.0.)                                GO TO 90
      SHIFTX = (X1(1) + X2(1)) * 0.5
      SHIFTY = (X1(2) + X2(2)) * 0.5
      SHIFTZ = (X1(3) + X2(3)) * 0.5
      SHIFTS = (S1 + S2) * 0.5
C
C             ONLY ONE VALUE NECESSARY SINCE X1= -X2 ETC.
C
      XSHFT  = X1(1) - SHIFTX
      YSHFT  = X1(2) - SHIFTY
      ZSHFT  = X1(3) - SHIFTZ
      SSHFT  = S1 - SHIFTS
C
      PABS1  = SQRT(X1(4)**2 + X1(5)**2 + X1(6)**2)
      PABS2  = SQRT(X2(4)**2 + X2(5)**2 + X2(6)**2)
      IF (PABS1.EQ.0..OR.PABS2.EQ.0.)            GO TO 90
C
C             PARAMETRIZE THE TRACK BY A CUBIC THROUGH X1,X2
C
      CALL GCUBS(SSHFT,XSHFT,X1(4)/PABS1,X2(4)/PABS2,A)
      CALL GCUBS(SSHFT,YSHFT,X1(5)/PABS1,X2(5)/PABS2,B)
      CALL GCUBS(SSHFT,ZSHFT,X1(6)/PABS1,X2(6)/PABS2,C)
C
C             ITERATE TO FIND THE TRACK LENGTH CORRESPONDING TO
C             THE INTERSECTION OF TRACK AND CYLINDER.
C             START AT S=0. MIDDLE OF THE SHIFTED INTERVAL.
C
      DINTER = ABS(S2 - S1) * 0.5
      S      = 0.
C
      DO 10 I = 1,MAXHIT
         X = SHIFTX + A(1) + S * (A(2) + S * (A(3) + S * A(4)))
         Y = SHIFTY + B(1) + S * (B(2) + S * (B(3) + S * B(4)))
         RN2    = X * X + Y * Y
         DR2    = (R2 - RN2) * DRCTN
         IF (ABS(DR2).LT.EPSI2)                     GO TO 20
         DINTER = DINTER * 0.5
         IF (DR2.LT.0.)S = S - DINTER
         IF (DR2.GE.0.)S = S + DINTER
  10  CONTINUE
C
C             COMPUTE INTERSECTION IN ORIGINAL COORDINATES
C
  20  CONTINUE
      XINT(1) = SHIFTX + A(1) + S * (A(2) + S * (A(3) + S * A(4)))
      XINT(2) = SHIFTY + B(1) + S * (B(2) + S * (B(3) + S * B(4)))
      XINT(3) = SHIFTZ + C(1) + S * (C(2) + S * (C(3) + S * C(4)))
      XINT(4) = A(2) + S * (2. * A(3) + 3. * S * A(4))
      XINT(5) = B(2) + S * (2. * B(3) + 3. * S * B(4))
      XINT(6) = C(2) + S * (2. * C(3) + 3. * S * C(4))
C
C             COMPUTE PHIHIT,ZHIT AND CORRESPONDING DERIVATIVES
C
      SINT   = S + SHIFTS
  200 TERM   = 1. / (XINT(4) * XINT(1) + XINT(5) * XINT(2))
      PZINT(1) = ATAN2(XINT(2),XINT(1))
      PZINT(2) = XINT(3)
      PZINT(3) = (XINT(1) * XINT(5) - XINT(2) * XINT(4)) * TERM / R
      PZINT(4) = TERM * XINT(6) * R
      RETURN
C
  90  IFLAG  = 0
      END
