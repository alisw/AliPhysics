*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:35  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_NEVER)
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GSYNIT
C
C  Preperation of synchrotron radiation tables.
C  This routine can be used to compute the YT and FTI used
C  in GSYNGE.
C  Author: Igor Gavrilenko
C
      COMMON /SYNTAB/ YT(54),FT(54),FTI(53)
      EXTERNAL BK53,BKS53
      DOUBLE PRECISION YM,YMA,EPS,XX
C...
      DATA EPS/.0001/
C...
      L = 0
      YS = .000001
      DO 20 I = 1,6
         YS = YS*10.
         DO 10 J = 1,9,1
            L = L+1
            YM = J*YS
            YT(L) = YM
            YMA = YM+YS
            IF(L.EQ.54) YMA = 20.
            FT(L) = DGAUSS(BK53,YM,YMA,EPS)
   10    CONTINUE
   20 CONTINUE
C...
      DO 30 I = 53,1,-1
         FT(I) = FT(I)+FT(I+1)
   30 CONTINUE
C...
      DO 40 I = 53,1,-1
         YM = YT(I)
         YMA = YT(I+1)
         IF(I.LT.53) FTI(I) = DGAUSS(BKS53,YM,YMA,EPS)+FTI(I+1)
         IF(I.EQ.53) FTI(I) = DGAUSS(BKS53,YM,YMA,EPS)
   40 CONTINUE
C...
      E N D
C-------------------------------------------------------
      DOUBLE PRECISION FUNCTION B K 5 3 (X)
C
      DOUBLE PRECISION X
C...
      BK53 = 4./(3.*X)*DBSKR3(X,2)+DBSKR3(X,1)
      RETURN
      E N D
C------------------------------------------------------
      DOUBLE PRECISION FUNCTION B K S 5 3 (X)
C
      COMMON /SYNTAB/ YT(54),FT(54),FTI(53)
      DOUBLE PRECISION X
C...
      IF(X.GT.YT(1)) GO TO 10
      BKS53 = FT(1)
      GO TO 60
C...
   10 IF(X.LT.YT(54)) GO TO 20
      BKS53 = 0.
      RETURN
C...
   20 K1 = 1
      K2 = 54
   30 K = (K1+K2)/2
      IF(X.GT.YT(K)) K1 = K
      IF(X.LT.YT(K)) K2 = K
      IF(X.NE.YT(K)) GO TO 40
      BKS53 = FT(K)
      GO TO 60
C...
   40 IF(K2-K1.GT.1) GO TO 30
      F1 = FT(K1)
      F2 = FT(K2)
      X1 = YT(K1)
      X2 = YT(K2)
C...
      IF(X.LE.1.) THEN
         X1 = LOG(X1)
         X2 = LOG(X2)
         B = F1+(F2-F1)/(X2-X1)*(LOG(X)-X1)
         GO TO 50
      ENDIF
C...
      IF(X.GT.1.) THEN
         F1 = LOG(F1)
         F2 = LOG(F2)
         B = F1+(F2-F1)/(X2-X1)*(X-X1)
         B = EXP(B)
      ENDIF
C...
   50 BKS53 = B
   60 BKS53 = .1909859*BKS53
      RETURN
      E N D
#endif
