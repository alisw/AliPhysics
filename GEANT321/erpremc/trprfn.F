*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:35  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE TRPRFN(X1,P1,H1,X2,P2,H2,CH,XL,R,MVAR,IFLAG,ITRAN,IERR)
************************************************************************
*
*
*     SUBR. TRPRFN(X1,P1,H1,X2,P2,H2,CH,XL,*R*,MVAR,IFLAG,ITRAN,IERR*)
*
*     Origin W.Wittek    EMCSW/81/18
*
*     Finite step length case coded by V.Innocente ( Feb. 88 )
*
*     code improved:                   V.Innocente ( April. 90 )
*                   inline code replaces external function
*     code improved:                   V.Innocente ( January 91 )
*                   effect of energy loss added
*
*_______________________________________________________________________
*
* *** ERROR PROPAGATION ALONG A PARTICLE TRAJECTORY IN A MAGNETIC FIELD
*     ROUTINE ASSUMES THAT IN THE INTERVAL (X1,X2) THE QUANTITIES 1/P
*     AND (HX,HY,HZ) ARE CONSTANT.
*
* *** IFLAG  =  -1   INITIALIZATION, TRANSFORMATION OF ERROR MATRIX FROM
*                    EXTERNAL TO SC VARIABLES
*            =   0   ERROR PROPAGATION FROM X1 TO X2
*            =   1   TRANSFORMATION OF ERROR MATRIX FROM SC TO
*                    EXTERNAL VARIABLES
*
*     ITRAN          USED FOR IFLAG = 0 OR 1 ONLY
*            =   0   TRANSFORMATION MATRIX IS UPDATED ,BUT ERROR MATRIX IS NOT
*                    TRANSFORMED
*           =    1   TRANSF. MATRIX IS UPDATED  AND ERROR MATRIX IS TRANSFORMED
*
*     MVAR           SPECIFIES TYPE OF EXTERNAL VARIABLES
*            =   0   ( 1/P,LAMBDA,PHI,YT, ZT ;   SC   )
*            =   1   ( 1/P,  Y',  Z',  Y,  Z ; SPLINE )
*
* *** X1, P1, H1     X,Y,Z COMPONENTS OF POSITION, MOMENTUM AND MAGNETIC   INPUT
*                    FIELD VECTOR/GRADIENT AT STARTING POINT OF INTERVAL
*     X2, P2, H2     ......  AT END POINT OF INTERVAL                      INPUT
*     CH             CHARGE OF PARTICLE                                    INPUT
*     XL             PATHLENGTH FROM X1 TO X2   ( NEGATIVE IF OPPOSITE
*                    TO ACTUAL MOVEMENT OF PARTICLE )                      INPUT
*     R              ERROR MATRIX  (TRIANGLE)                       INPUT/OUTPUT
*     B              5 * 5 TRANSFORMATION MATRIX FOR ERRORS IN
*                    SC VARIABLES                                         OUTPUT
*
* *** IERR   =  1    ILLEGAL VALUE OF MVAR                                OUTPUT
*               2    MOMENTUM IS ZERO
*               3    H*ALFA/P AT X1 AND X2 DIFFER TOO MUCH
*               4    PARTICLE MOVES IN Z - DIRECTION
*
************************************************************************
*
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL  X1,P1,H1,X2,P2,H2,R,CH,PS,PC,XL,SPX
#endif
#include "geant321/trcom3.inc"
#include "geant321/gcunit.inc"
      DIMENSION X1(3),P1(3),H1(9),X2(3),P2(3),H2(9)
      DIMENSION R(15),PS(3),PC(3)
*
      DIMENSION T1(3),T2(3),U1(3),U2(3),V1(3),V2(3),HN(9)
      DIMENSION AN1(3),AN2(3),DX(3)
      DIMENSION HV1(3),HU1(3)
*
      SAVE INIT,DELHP6,CFACT8
*
      DATA INIT/0/
#if !defined(CERNLIB_SINGLE)
      DATA DELHP6/300.D0/
*
      DATA CFACT8 / 2.997925 D-4 /
#endif
#if defined(CERNLIB_SINGLE)
      DATA DELHP6/300./
*
      DATA CFACT8 / 2.997925 E-4 /
#endif
*
*____________________________________________________________________
*
      IERR=0
      IF(IFLAG) 10, 20, 80
*
* *** TRANSFORM ERROR MATRIX FROM EXTERNAL TO INTERNAL VARIABLES;
*
   10 NEW=1
      IF(MVAR.NE.1) GO TO 11
      PA1=SQRT(P1(1)**2+P1(2)**2+P1(3)**2)
      IF(PA1.EQ.0.) GO TO 902
      PS(1)=1./PA1
      IF(P1(1).EQ.0.) GO TO 904
      PS(2)=P1(2)/P1(1)
      PS(3)=P1(3)/P1(1)
      SPX=1.
      IF(P1(1).LT.0.) SPX=-1.
      CALL TRSPSC(PS,R,PC,R,H1,CH,IERR,SPX)
      GO TO 19
*
   11 IF(MVAR.NE.0) GO TO 901
   19 GO TO 900
*
* *** ERROR PROPAGATION ON A HELIX ASSUMING SC VARIABLES
 
*
   20 PA1=SQRT(P1(1)**2+P1(2)**2+P1(3)**2)
      PA2=SQRT(P2(1)**2+P2(2)**2+P2(3)**2)
      IF(PA1*PA2.EQ.0.) GO TO 902
C      DPA = PA2 - PA1
      PM1=1./PA1
      PM2=1./PA2
      DPM = PM2 - PM1
*
      DO 201 I=1,3
        T1(I) = P1(I)*PM1
        T2(I) = P2(I)*PM2
201   CONTINUE
*
      SINL=T2(3)
      SINL0=T1(3)
*
      COSL=SQRT(ABS(1.-SINL**2))
      IF(COSL.EQ.0.) GO TO 904
      COSL1=1./COSL
      COSL0=SQRT(ABS(1.-SINL0**2))
*
* *** DEFINE TRANSFORMATION MATRIX BETWEEN X1 AND X2 FOR
* *** NEUTRAL PARTICLE OR FIELDFREE REGION
*
      DO 26 I=1,5
         DO 15 K=1,5
            A(I,K)=0.
   15    CONTINUE
         A(I,I)=1.
   26 CONTINUE
      A(4,3)=XL*COSL
      A(5,2)=XL
*
      IF(CH.EQ.0.) GO TO 45
      HA1=SQRT(H1(1)**2+H1(2)**2+H1(3)**2)
      HA2=SQRT(H2(1)**2+H2(2)**2+H2(3)**2)
      HAM1=HA1*PM1
      HAM2=HA2*PM2
      HAMX=MAX(HAM1,HAM2)
      IF(HAMX.EQ.0.) GO TO 45
*
*
*
* *** CHECK WHETHER H*ALFA/P IS TOO DIFFERENT AT X1 AND X2
*
*
      IF(HA2.NE.0.) THEN
         GAM=(H2(1)*T2(1)+H2(2)*T2(2)+H2(3)*T2(3))/HA2
      ELSE
         GAM=(H1(1)*T1(1)+H1(2)*T1(2)+H1(3)*T1(3))/HA1
      ENDIF
*
      ALFA2=1.-GAM**2
*
      DH2=(H1(1)*PM1-H2(1)*PM2)**2+
     1    (H1(2)*PM1-H2(2)*PM2)**2+
     1    (H1(3)*PM1-H2(3)*PM2)**2
      IF(DH2*ALFA2.GT.DELHP6**2) GO TO 903
*
* *** DEFINE AVERAGE MAGNETIC FIELD AND GRADIENT
*
      PM12=(PM1+PM2)*0.5
      P12=1./(2.*PM12)
      HN(1)=(H1(1)*PM1+H2(1)*PM2)*P12*CH*CFACT8
      HN(2)=(H1(2)*PM1+H2(2)*PM2)*P12*CH*CFACT8
      HN(3)=(H1(3)*PM1+H2(3)*PM2)*P12*CH*CFACT8
CC    HN(4)=(H1(4)*PM1+H2(4)*PM2)*P12*CH*CFACT8
CC    HN(5)=(H1(5)*PM1+H2(5)*PM2)*P12*CH*CFACT8
CC    HN(6)=(H1(6)*PM1+H2(6)*PM2)*P12*CH*CFACT8
CC    HN(7)=(H1(7)*PM1+H2(7)*PM2)*P12*CH*CFACT8
CC    HN(8)=(H1(8)*PM1+H2(8)*PM2)*P12*CH*CFACT8
CC    HN(9)=(H1(9)*PM1+H2(9)*PM2)*P12*CH*CFACT8
*
      HM = SQRT(HN(1)**2+HN(2)**2+HN(3)**2)
      OVER = 1./HM
      HN(1) = OVER*HN(1)
      HN(2) = OVER*HN(2)
      HN(3) = OVER*HN(3)
      PAV = .5*(PA1+PA2)
      Q = - HM/PAV
      THETA = Q*XL
      SINT = SIN(THETA)
      COST = COS(THETA)
      GAMMA=HN(1)*T2(1)+HN(2)*T2(2)+HN(3)*T2(3)
      AN2(1) = HN(2)*T2(3)-HN(3)*T2(2)
      AN2(2) = HN(3)*T2(1)-HN(1)*T2(3)
      AN2(3) = HN(1)*T2(2)-HN(2)*T2(1)
*
      AU = 1./SQRT(T1(1)**2+T1(2)**2)
      U1(1) = -AU*T1(2)
      U1(2) =  AU*T1(1)
      U1(3) =  0.D0
      V1(1) = -T1(3)*U1(2)
      V1(2) =  T1(3)*U1(1)
      V1(3) =  T1(1)*U1(2)-T1(2)*U1(1)
*
      AU = 1./SQRT(T2(1)**2+T2(2)**2)
      U2(1) = -AU*T2(2)
      U2(2) =  AU*T2(1)
      U2(3) =  0.D0
      V2(1) = -T2(3)*U2(2)
      V2(2) =  T2(3)*U2(1)
      V2(3) =  T2(1)*U2(2)-T2(2)*U2(1)
*
      DX(1) = X1(1) - X2(1)
      DX(2) = X1(2) - X2(2)
      DX(3) = X1(3) - X2(3)
*
*
* *** COMPLETE TRANSFORMATION MATRIX BETWEEN ERRORS AT X1 AND X2
* *** FIELD GRADIENT PERPENDICULAR TO TRACK IS PRESENTLY NOT
* *** TAKEN INTO ACCOUNT
*
   30 CONTINUE
      QP  = Q*PAV
      ANV = -(HN(1)*U2(1)+HN(2)*U2(2)            )
      ANU =  (HN(1)*V2(1)+HN(2)*V2(2)+HN(3)*V2(3))
      OMCOST = 1.-COST
      TMSINT = THETA-SINT
*
      HU1(1) =            -HN(3)*U1(2)
      HU1(2) = HN(3)*U1(1)
      HU1(3) = HN(1)*U1(2)-HN(2)*U1(1)
*
      HV1(1) = HN(2)*V1(3)-HN(3)*V1(2)
      HV1(2) = HN(3)*V1(1)-HN(1)*V1(3)
      HV1(3) = HN(1)*V1(2)-HN(2)*V1(1)
*
***   1/P
*
      A(1,1) = 1.-DPM*PAV*(1.+(T2(1)*DX(1)+T2(2)*DX(2)+T2(3)*DX(3))/XL)
     +           +2.*DPM*PAV
*
      A(1,2) =  -DPM/THETA*
     1           ( TMSINT*GAMMA*(HN(1)*V1(1)+HN(2)*V1(2)+HN(3)*V1(3)) +
     2             SINT*(V1(1)*T2(1)+V1(2)*T2(2)+V1(3)*T2(3)) +
     3             OMCOST*(HV1(1)*T2(1)+HV1(2)*T2(2)+HV1(3)*T2(3)) )
*
      A(1,3) =  -COSL0*DPM/THETA*
     1           ( TMSINT*GAMMA*(HN(1)*U1(1)+HN(2)*U1(2)            ) +
     2             SINT*(U1(1)*T2(1)+U1(2)*T2(2)            ) +
     3             OMCOST*(HU1(1)*T2(1)+HU1(2)*T2(2)+HU1(3)*T2(3)) )
*
      A(1,4) =  -DPM/XL*(U1(1)*T2(1)+U1(2)*T2(2)            )
*
      A(1,5) =  -DPM/XL*(V1(1)*T2(1)+V1(2)*T2(2)+V1(3)*T2(3))
*
***   Lambda
*
      A(2,1) = -QP*ANV*(T2(1)*DX(1)+T2(2)*DX(2)+T2(3)*DX(3))
     +         *(1.+DPM*PAV)
*
      A(2,2) = COST*(V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)) +
     +         SINT*(HV1(1)*V2(1)+HV1(2)*V2(2)+HV1(3)*V2(3)) +
     1         OMCOST*(HN(1)*V1(1)+HN(2)*V1(2)+HN(3)*V1(3))*
     A                (HN(1)*V2(1)+HN(2)*V2(2)+HN(3)*V2(3)) +
     2         ANV*( -SINT*(V1(1)*T2(1)+V1(2)*T2(2)+V1(3)*T2(3)) +
     3         OMCOST*(V1(1)*AN2(1)+V1(2)*AN2(2)+V1(3)*AN2(3)) -
     4         TMSINT*GAMMA*(HN(1)*V1(1)+HN(2)*V1(2)+HN(3)*V1(3)) )
*
      A(2,3) = COST*(U1(1)*V2(1)+U1(2)*V2(2)            ) +
     +         SINT*(HU1(1)*V2(1)+HU1(2)*V2(2)+HU1(3)*V2(3)) +
     1         OMCOST*(HN(1)*U1(1)+HN(2)*U1(2)            )*
     A                (HN(1)*V2(1)+HN(2)*V2(2)+HN(3)*V2(3)) +
     2         ANV*( -SINT*(U1(1)*T2(1)+U1(2)*T2(2)            ) +
     3         OMCOST*(U1(1)*AN2(1)+U1(2)*AN2(2)             ) -
     4         TMSINT*GAMMA*(HN(1)*U1(1)+HN(2)*U1(2)            ) )
      A(2,3) = COSL0*A(2,3)
*
      A(2,4) = -Q*ANV*(U1(1)*T2(1)+U1(2)*T2(2)            )
*
      A(2,5) = -Q*ANV*(V1(1)*T2(1)+V1(2)*T2(2)+V1(3)*T2(3))
*
***   Phi
*
      A(3,1) = -QP*ANU*(T2(1)*DX(1)+T2(2)*DX(2)+T2(3)*DX(3))*COSL1
     +         *(1.+DPM*PAV)
*
      A(3,2) = COST*(V1(1)*U2(1)+V1(2)*U2(2)            ) +
     +         SINT*(HV1(1)*U2(1)+HV1(2)*U2(2)             ) +
     1         OMCOST*(HN(1)*V1(1)+HN(2)*V1(2)+HN(3)*V1(3))*
     A                (HN(1)*U2(1)+HN(2)*U2(2)            ) +
     2         ANU*( -SINT*(V1(1)*T2(1)+V1(2)*T2(2)+V1(3)*T2(3)) +
     3         OMCOST*(V1(1)*AN2(1)+V1(2)*AN2(2)+V1(3)*AN2(3)) -
     4         TMSINT*GAMMA*(HN(1)*V1(1)+HN(2)*V1(2)+HN(3)*V1(3)) )
      A(3,2) = COSL1*A(3,2)
*
      A(3,3) = COST*(U1(1)*U2(1)+U1(2)*U2(2)            ) +
     +         SINT*(HU1(1)*U2(1)+HU1(2)*U2(2)             ) +
     1         OMCOST*(HN(1)*U1(1)+HN(2)*U1(2)            )*
     A                (HN(1)*U2(1)+HN(2)*U2(2)            ) +
     2         ANU*( -SINT*(U1(1)*T2(1)+U1(2)*T2(2)            ) +
     3         OMCOST*(U1(1)*AN2(1)+U1(2)*AN2(2)             ) -
     4         TMSINT*GAMMA*(HN(1)*U1(1)+HN(2)*U1(2)            ) )
      A(3,3) = COSL1*COSL0*A(3,3)
*
      A(3,4) = -Q*ANU*(U1(1)*T2(1)+U1(2)*T2(2)            )*COSL1
*
      A(3,5) = -Q*ANU*(V1(1)*T2(1)+V1(2)*T2(2)+V1(3)*T2(3))*COSL1
*
***   Yt
*
      A(4,1) = PAV*(U2(1)*DX(1)+U2(2)*DX(2)            )
     +         *(1.+DPM*PAV)
*
      A(4,2) = (   SINT*(V1(1)*U2(1)+V1(2)*U2(2)            ) +
     1           OMCOST*(HV1(1)*U2(1)+HV1(2)*U2(2)             ) +
     2           TMSINT*(HN(1)*U2(1)+HN(2)*U2(2)            )*
     3                  (HN(1)*V1(1)+HN(2)*V1(2)+HN(3)*V1(3)) )/Q
*
      A(4,3) = (   SINT*(U1(1)*U2(1)+U1(2)*U2(2)            ) +
     1           OMCOST*(HU1(1)*U2(1)+HU1(2)*U2(2)             ) +
     2           TMSINT*(HN(1)*U2(1)+HN(2)*U2(2)            )*
     3                  (HN(1)*U1(1)+HN(2)*U1(2)            ) )*COSL0/Q
*
      A(4,4) = (U1(1)*U2(1)+U1(2)*U2(2)            )
*
      A(4,5) = (V1(1)*U2(1)+V1(2)*U2(2)            )
*
***   Zt
*
      A(5,1) = PAV*(V2(1)*DX(1)+V2(2)*DX(2)+V2(3)*DX(3))
     +         *(1.+DPM*PAV)
*
      A(5,2) = (   SINT*(V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)) +
     1           OMCOST*(HV1(1)*V2(1)+HV1(2)*V2(2)+HV1(3)*V2(3)) +
     2           TMSINT*(HN(1)*V2(1)+HN(2)*V2(2)+HN(3)*V2(3))*
     3                  (HN(1)*V1(1)+HN(2)*V1(2)+HN(3)*V1(3)) )/Q
*
      A(5,3) = (   SINT*(U1(1)*V2(1)+U1(2)*V2(2)            ) +
     1           OMCOST*(HU1(1)*V2(1)+HU1(2)*V2(2)+HU1(3)*V2(3)) +
     2           TMSINT*(HN(1)*V2(1)+HN(2)*V2(2)+HN(3)*V2(3))*
     3                  (HN(1)*U1(1)+HN(2)*U1(2)            ) )*COSL0/Q
*
      A(5,4) = (U1(1)*V2(1)+U1(2)*V2(2)            )
*
      A(5,5) = (V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))
   45 CONTINUE
*
* *** NEW = 0  TRANSFORMATION MATRIX IS UPDATED
*           1  TRANSFORMATION MATRIX IS INITIALIZED
*
      IF(NEW.EQ.0) GO TO 23
      NEW=0
      DO 25 I=1,5
         DO 24 K=1,5
            B(I,K)=A(I,K)
   24    CONTINUE
   25 CONTINUE
      GO TO 27
   23 CONTINUE
*
      CALL XMM55(A,B,B)
*
   27 CONTINUE
   80 IF(ITRAN.EQ.0) GO TO 90
*
*
      J=0
      DO 22 I=1,5
         DO 21 K=I,5
            J=J+1
            S(J)=R(J)
   21    CONTINUE
   22 CONTINUE
*
*
* *** TRANSFORM ERROR MATRIX
*
      CALL SSMT5T(B,S,S)
*
      NEW=1
 
      J=0
      DO 41 I=1,5
         DO 40 K=I,5
            J=J+1
            R(J)=S(J)
   40    CONTINUE
   41 CONTINUE
*
   90 IF(IFLAG.LE.0) GO TO 900
*
*
* *** TRANSFORM ERROR MATRIX FROM INTERNAL TO EXTERNAL VARIABLES;
*
*
      NEW=1
      IF(MVAR.NE.1) GO TO 91
      PC(1)=PM2
      PC(2)=ASIN(P2(3)*PC(1))
      IF (ABS (P2(1)) .LT. 1.E-30) P2(1) = 1.E-30
      PC(3)=ATAN2(P2(2),P2(1))
      CALL TRSCSP(PC,R,PS,R,H2,CH,IERR,SPX)
      GO TO 900
*
   91 IF(MVAR.NE.0) GO TO 901
      GO TO 900
*
* *** ERROR EXITS
*
  901 IERR=1
      GO TO 999
  902 IERR=2
      GO TO 999
  903 IERR=3
C     IF(INIT.NE.0) GO TO 30
*     WRITE (LOUT, 998) DH2,ALFA2,XL
  998 FORMAT('0',' *** S/R TRPROP   DELTA(H*ALFA/P)',5X
     1,'EXCEEDS TOLERANCE    '/'0',3E12.5//' **********    ',///)
      INIT=1
      GO TO 30
  904 IERR=4
  999 WRITE (LOUT, 1000) IERR
 1000 FORMAT(1H ,' *** S/R ERPROP   IERR =',I5)
*
  900 CONTINUE
      END
 
