*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:25  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GLANDG(YRAN)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Copy of the CERN library routine GENLAN                       *
C.    *  Generation of LANDAU-distributed random numbers by 4-point    *
C.    *  interpolation in the previously-tabulated inverse cumulative  *
C.    *  distribution                                                  *
C.    *                                                                *
C.    *    ==>Called by : GLANDO                                       *
C.    *                                                                *
C.    ******************************************************************
C.
      PARAMETER (RANGE1=0.807069,RANGE2=0.994869)
      PARAMETER (X1BOT =0.0     ,X2BOT =0.791240)
      PARAMETER (GP1INV=118.836 ,GP2INV=486.178 )
C
      DIMENSION XCUM1(100) , XCUM2(100)
      DIMENSION RNDM(1)
C
C====>  1ST TABLE OF INVERSE CUMULATIVE LANDAU POINTS. ( 0.0<P<0.833 )
      DATA XCUM1(  1) / -2.57000 /   ,   XCUM1(  2) / -2.15412 /
      DATA XCUM1(  3) / -1.94167 /   ,   XCUM1(  4) / -1.79583 /
      DATA XCUM1(  5) / -1.67975 /   ,   XCUM1(  6) / -1.58090 /
      DATA XCUM1(  7) / -1.49341 /   ,   XCUM1(  8) / -1.41397 /
      DATA XCUM1(  9) / -1.34057 /   ,   XCUM1( 10) / -1.27185 /
      DATA XCUM1( 11) / -1.20686 /   ,   XCUM1( 12) / -1.14490 /
      DATA XCUM1( 13) / -1.08545 /   ,   XCUM1( 14) / -1.02809 /
      DATA XCUM1( 15) /  -.97249 /   ,   XCUM1( 16) /  -.91839 /
      DATA XCUM1( 17) /  -.86557 /   ,   XCUM1( 18) /  -.81386 /
      DATA XCUM1( 19) /  -.76308 /   ,   XCUM1( 20) /  -.71311 /
      DATA XCUM1( 21) /  -.66384 /   ,   XCUM1( 22) /  -.61515 /
      DATA XCUM1( 23) /  -.56696 /   ,   XCUM1( 24) /  -.51919 /
      DATA XCUM1( 25) /  -.47176 /   ,   XCUM1( 26) /  -.42462 /
      DATA XCUM1( 27) /  -.37769 /   ,   XCUM1( 28) /  -.33092 /
      DATA XCUM1( 29) /  -.28426 /   ,   XCUM1( 30) /  -.23767 /
      DATA XCUM1( 31) /  -.19109 /   ,   XCUM1( 32) /  -.14447 /
      DATA XCUM1( 33) /  -.09779 /   ,   XCUM1( 34) /  -.05100 /
      DATA XCUM1( 35) /  -.00405 /   ,   XCUM1( 36) /   .04308 /
      DATA XCUM1( 37) /   .09044 /   ,   XCUM1( 38) /   .13806 /
      DATA XCUM1( 39) /   .18598 /   ,   XCUM1( 40) /   .23423 /
      DATA XCUM1( 41) /   .28285 /   ,   XCUM1( 42) /   .33187 /
      DATA XCUM1( 43) /   .38134 /   ,   XCUM1( 44) /   .43128 /
      DATA XCUM1( 45) /   .48174 /   ,   XCUM1( 46) /   .53275 /
      DATA XCUM1( 47) /   .58435 /   ,   XCUM1( 48) /   .63658 /
      DATA XCUM1( 49) /   .68948 /   ,   XCUM1( 50) /   .74309 /
      DATA XCUM1( 51) /   .79746 /   ,   XCUM1( 52) /   .85263 /
      DATA XCUM1( 53) /   .90865 /   ,   XCUM1( 54) /   .96556 /
      DATA XCUM1( 55) /  1.02342 /   ,   XCUM1( 56) /  1.08228 /
      DATA XCUM1( 57) /  1.14219 /   ,   XCUM1( 58) /  1.20322 /
      DATA XCUM1( 59) /  1.26542 /   ,   XCUM1( 60) /  1.32887 /
      DATA XCUM1( 61) /  1.39362 /   ,   XCUM1( 62) /  1.45976 /
      DATA XCUM1( 63) /  1.52736 /   ,   XCUM1( 64) /  1.59650 /
      DATA XCUM1( 65) /  1.66727 /   ,   XCUM1( 66) /  1.73976 /
      DATA XCUM1( 67) /  1.81407 /   ,   XCUM1( 68) /  1.89032 /
      DATA XCUM1( 69) /  1.96860 /   ,   XCUM1( 70) /  2.04905 /
      DATA XCUM1( 71) /  2.13180 /   ,   XCUM1( 72) /  2.21699 /
      DATA XCUM1( 73) /  2.30477 /   ,   XCUM1( 74) /  2.39531 /
      DATA XCUM1( 75) /  2.48878 /   ,   XCUM1( 76) /  2.58540 /
      DATA XCUM1( 77) /  2.68535 /   ,   XCUM1( 78) /  2.78889 /
      DATA XCUM1( 79) /  2.89626 /   ,   XCUM1( 80) /  3.00775 /
      DATA XCUM1( 81) /  3.12364 /   ,   XCUM1( 82) /  3.24429 /
      DATA XCUM1( 83) /  3.37005 /   ,   XCUM1( 84) /  3.50136 /
      DATA XCUM1( 85) /  3.63866 /   ,   XCUM1( 86) /  3.78246 /
      DATA XCUM1( 87) /  3.93334 /   ,   XCUM1( 88) /  4.09194 /
      DATA XCUM1( 89) /  4.25900 /   ,   XCUM1( 90) /  4.43533 /
      DATA XCUM1( 91) /  4.62186 /   ,   XCUM1( 92) /  4.81960 /
      DATA XCUM1( 93) /  5.02974 /   ,   XCUM1( 94) /  5.25368 /
      DATA XCUM1( 95) /  5.49312 /   ,   XCUM1( 96) /  5.74987 /
      DATA XCUM1( 97) /  6.02605 /   ,   XCUM1( 98) /  6.32428 /
      DATA XCUM1( 99) /  6.64773 /   ,   XCUM1(100) /  7.00000 /
C
C====>  2ND TABLE OF INVERSE CUMULATIVE LANDAU POINTS. ( 0.791<P<0.995 )
      DATA XCUM2(  1) /   5.50000 /   ,   XCUM2(  2) /   5.56120 /
      DATA XCUM2(  3) /   5.62347 /   ,   XCUM2(  4) /   5.68684 /
      DATA XCUM2(  5) /   5.75133 /   ,   XCUM2(  6) /   5.81699 /
      DATA XCUM2(  7) /   5.88383 /   ,   XCUM2(  8) /   5.95191 /
      DATA XCUM2(  9) /   6.02125 /   ,   XCUM2( 10) /   6.09190 /
      DATA XCUM2( 11) /   6.16391 /   ,   XCUM2( 12) /   6.23732 /
      DATA XCUM2( 13) /   6.31219 /   ,   XCUM2( 14) /   6.38855 /
      DATA XCUM2( 15) /   6.46646 /   ,   XCUM2( 16) /   6.54597 /
      DATA XCUM2( 17) /   6.62712 /   ,   XCUM2( 18) /   6.70998 /
      DATA XCUM2( 19) /   6.79460 /   ,   XCUM2( 20) /   6.88103 /
      DATA XCUM2( 21) /   6.96935 /   ,   XCUM2( 22) /   7.05962 /
      DATA XCUM2( 23) /   7.15192 /   ,   XCUM2( 24) /   7.24633 /
      DATA XCUM2( 25) /   7.34294 /   ,   XCUM2( 26) /   7.44182 /
      DATA XCUM2( 27) /   7.54306 /   ,   XCUM2( 28) /   7.64676 /
      DATA XCUM2( 29) /   7.75300 /   ,   XCUM2( 30) /   7.86188 /
      DATA XCUM2( 31) /   7.97351 /   ,   XCUM2( 32) /   8.08800 /
      DATA XCUM2( 33) /   8.20548 /   ,   XCUM2( 34) /   8.32610 /
      DATA XCUM2( 35) /   8.44997 /   ,   XCUM2( 36) /   8.57725 /
      DATA XCUM2( 37) /   8.70808 /   ,   XCUM2( 38) /   8.84262 /
      DATA XCUM2( 39) /   8.98103 /   ,   XCUM2( 40) /   9.12349 /
      DATA XCUM2( 41) /   9.27021 /   ,   XCUM2( 42) /   9.42141 /
      DATA XCUM2( 43) /   9.57730 /   ,   XCUM2( 44) /   9.73812 /
      DATA XCUM2( 45) /   9.90410 /   ,   XCUM2( 46) /  10.07552 /
      DATA XCUM2( 47) /  10.25265 /   ,   XCUM2( 48) /  10.43584 /
      DATA XCUM2( 49) /  10.62540 /   ,   XCUM2( 50) /  10.82169 /
      DATA XCUM2( 51) /  11.02508 /   ,   XCUM2( 52) /  11.23598 /
      DATA XCUM2( 53) /  11.45487 /   ,   XCUM2( 54) /  11.68222 /
      DATA XCUM2( 55) /  11.91855 /   ,   XCUM2( 56) /  12.16441 /
      DATA XCUM2( 57) /  12.42045 /   ,   XCUM2( 58) /  12.68734 /
      DATA XCUM2( 59) /  12.96580 /   ,   XCUM2( 60) /  13.25663 /
      DATA XCUM2( 61) /  13.56075 /   ,   XCUM2( 62) /  13.87912 /
      DATA XCUM2( 63) /  14.21280 /   ,   XCUM2( 64) /  14.56296 /
      DATA XCUM2( 65) /  14.93088 /   ,   XCUM2( 66) /  15.31799 /
      DATA XCUM2( 67) /  15.72593 /   ,   XCUM2( 68) /  16.15650 /
      DATA XCUM2( 69) /  16.61170 /   ,   XCUM2( 70) /  17.09378 /
      DATA XCUM2( 71) /  17.60523 /   ,   XCUM2( 72) /  18.14884 /
      DATA XCUM2( 73) /  18.72792 /   ,   XCUM2( 74) /  19.34622 /
      DATA XCUM2( 75) /  20.00796 /   ,   XCUM2( 76) /  20.71792 /
      DATA XCUM2( 77) /  21.48176 /   ,   XCUM2( 78) /  22.30618 /
      DATA XCUM2( 79) /  23.19862 /   ,   XCUM2( 80) /  24.16809 /
      DATA XCUM2( 81) /  25.22547 /   ,   XCUM2( 82) /  26.38320 /
      DATA XCUM2( 83) /  27.65658 /   ,   XCUM2( 84) /  29.06490 /
      DATA XCUM2( 85) /  30.63048 /   ,   XCUM2( 86) /  32.38211 /
      DATA XCUM2( 87) /  34.35555 /   ,   XCUM2( 88) /  36.59607 /
      DATA XCUM2( 89) /  39.16212 /   ,   XCUM2( 90) /  42.13195 /
      DATA XCUM2( 91) /  45.61204 /   ,   XCUM2( 92) /  49.74758 /
      DATA XCUM2( 93) /  54.74189 /   ,   XCUM2( 94) /  60.89935 /
      DATA XCUM2( 95) /  68.67370 /   ,   XCUM2( 96) /  78.81480 /
      DATA XCUM2( 97) /  92.61047 /   ,   XCUM2( 98) / 112.50807 /
      DATA XCUM2( 99) / 143.78539 /   ,   XCUM2(100) / 200.00000 /
C.
C.    ------------------------------------------------------------------
C.
      CALL GRNDM(RNDM,1)
      X = RNDM(1)
      IF(X .LT. 0.004) THEN
*
* *** Extreme left-hand tail
         YRAN = -SQRT(ABS(LOG(X)))
*
      ELSEIF(X.LE.RANGE1) THEN
*
* *** 4-point interpolation in the first cumulative table
         TABPO1 = (X-X1BOT)*GP1INV
         J = TABPO1 + 1
         J = MAX(J,2)
         J = MIN(J,98)
         P = TABPO1 - (J-1)
         A = (P+1.0) * XCUM1(J+2) - (P-2.0) * XCUM1(J-1)
         B = (P-1.0) * XCUM1(J)   - P       * XCUM1(J+1)
         YRAN = A*P*(P-1.0)*0.16666667+B*(P+1.0)*(P-2.0)*0.5
*
      ELSEIF(X.LE.RANGE2) THEN
*
* *** 4-point interpolation in the first cumulative table
         TABPO2 = (X-X2BOT)*GP2INV
         J = TABPO2 + 1
         J = MAX(J,2)
         J = MIN(J,98)
         P = TABPO2 - (J-1)
         A = (P+1.0) * XCUM2(J+2) - (P-2.0) * XCUM2(J-1)
         B = (P-1.0) * XCUM2(J)   - P       * XCUM2(J+1)
         YRAN = A*P*(P-1.0)*0.16666667+B*(P+1.0)*(P-2.0)*0.5
*
      ELSE
*
* *** 1/x**2 sampling for extreme Landau tail
         CALL GRNDM(RNDM,1)
         YRAN=200./RNDM(1)
*
      ENDIF
      END
