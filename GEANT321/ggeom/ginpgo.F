*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GINPGO(X,P,IYES)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE TO PERFORM GINME FUNCTION FOR A POLYGON.            *
C.    *        X  -  A VECTOR DIMENSION 3 GIVING THE POINT.            *
C.    *        P  -  THE VECTOR OF PARAMETERS, (1) = PHIMIN,           *
C.    *                   (2) = DPHI, (3) = NDIV, (4) = NZ, (5) = Z(1) *
C.    *                   (6) = RMIN(1), (7) = RMAX(1), (8) = Z(2),    *
C.    *                   (9) = RMIN(2), (10) = RMAX(2), ETC.          *
C.    *                   RMIN IS THE DISTANCE FROM THE ORIGIN TO THE  *
C.    *                   INSIDE EDGE PROJECTED ALONG THE CENTRE LINE  *
C.    *                   OF THE DIVISION. RMAX IS THE CORRESPONDING   *
C.    *                   PROJECTED DISTANCE FOR THE OUTER EDGE.       *
C.    *     IYES  -  RETURNED AS 1 IF POINT IS INSIDE, OTHERWISE 0.    *
C.    *                                                                *
C.    *    ==>Called by : GINME                                        *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcpoly.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION XX1,XX2,PHI,DELPHI,PH0,R,DZRAT,DIVPHI
      DOUBLE PRECISION RMN,RMX,HNDIV
#endif
      DIMENSION P(10),X(3)
C.
C.                ---------------------------------
C.
      IYES=0
C
C                 FIRST CHECK TOTAL Z RANGE.
C
      NZ=P(4)
      IF(X(3).LT.P(5).OR.X(3).GT.P(NZ*3+2)) GO TO 999
C
C             Find smallest RMIN and largest RMAX
C
      RMIN=P(6)
      RMAX=P(7)
      DO 10 I=2,NZ
         IP=I*3+3
         IF(P(IP).LT.RMIN)RMIN=P(IP)
         IF(P(IP+1).GT.RMAX) RMAX=P(IP+1)
   10 CONTINUE
*
* ***         Check R against RMIN
*
      XX1 = X(1)
      XX2 = X(2)
      R   = SQRT(XX1**2+XX2**2)
      IF(R.LT.RMIN)GO TO 999
      HNDIV=P(3)
      DELPHI=P(2)
      DIVPHI=DELPHI/HNDIV
*
* ***     Find the radius of the outscribed circle
*
      RMAX = RMAX/COS(0.5*DIVPHI*DEGRAD)
*
* ***         Check R against RMAX
*
      IF(R.GT.RMAX)GO TO 999
C
C               If R sufficiently small ignore phi structure.
C
      IPSEC=1
C
      IF(R.LT.1.0E-5) GO TO 20
C
C                 NOW CHECK PHI.
C
      PHI=ATAN2(XX2,XX1)*RADDEG
      IF(PHI.LT.P(1)) PHI=PHI+360.0
      IF(PHI.LT.P(1).OR.PHI-P(1).GT.P(2)) GO TO 999
C
C                 NOW FIND PHI DIVISION.
C
      IPSEC=MIN((PHI-P(1))/DIVPHI+1.,HNDIV)
C
      PH0=(P(1)+DIVPHI*(IPSEC-0.5))*DEGRAD
C
C                 NOW CHECK PROJECTED DISTANCE. FIRST COMPUTE IT.
C
      R=XX1*COS(PH0)+XX2*SIN(PH0)
C
   20 CONTINUE
C
C                 NOW FIND WHICH Z SECTION THE POINT IS IN.
C
      IZL=1
      IZH=NZ
      IPL=IZL*3+2
      IPH=IZH*3+2
C
   30 IF(IZH-IZL.GT.1) THEN
C
         IZT=(IZL+IZH)/2
         IPT=IZT*3+2
         IF(X(3).LE.P(IPT)) THEN
C
            IZH=IZT
            IPH=IPT
         ELSE
C
            IZL=IZT
            IPL=IPT
         ENDIF
C
         GO TO 30
      ENDIF
C
C          THE POINT IS IN THE SECTION BOUNDED BY THE IZL TH AND THE
C          IZH TH Z VALUES.
C
C          NOW COMPUTE RMIN AND RMAX AND TEST THE VALUE OF R .
C
      IZSEC=IZL
      DZRAT=X(3)-P(IPL)
      DZRAT=DZRAT/(P(IPH)-P(IPL))
      RMN=P(IPL+1)+DZRAT*(P(IPH+1)-P(IPL+1))
*
* ***     Is the point inside the 'hole' at the centre of the volume?
*
      IF(R.LT.RMN) GO TO 999
*
      RMX=P(IPL+2)+DZRAT*(P(IPH+2)-P(IPL+2))
*
* ***     Is the point outside the volume?
*
      IF(R.GT.RMX) GO TO 999
C
C           POINT IS IN POLYGON.
C
      IYES=1
C
  999 END
