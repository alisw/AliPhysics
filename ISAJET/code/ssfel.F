#include "isajet/pilot.h"
      FUNCTION SSFEL(X,INIT)
C***********************************************************************
C* Computes the electron spectrum as a convolution of the beam- and    *
C* bremsstrahlung-spectra, including leading-log summation for the lat-*
C* ter (in one-loop order), and Chen's approximate expression for the  *
C* former. X is the e energy in units of the nominal beam energy, and  *
C* BETA is 2 alpha_em / pi (log s/me^2 - 1). If more than 99.5% of all *
C* electrons are in the delta-peak, beamstrahlung is ignored. Other-   *
C* wise, beamstrahlung is included. In the latter case, the complete   *
C* spectrum is computed at the first call (with INIT=1), and fitted in *
C* a cubic spline; in later calls (with INIT=0), only the spline is    *
C* used. This reduces the necessary amount of CPU time considerably.   *
C* This subroutine needs the programs BEAMEL, SIMAU8, and SPLINE.      *
C***********************************************************************
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/eepar.inc"
#include "isajet/brembm.inc"
C
      REAL X,SSFEL
      INTEGER INIT
      REAL Y,XLMM,XL,GAM,RE,XKAPPA,NUCL,NUGAM,NGAM,DC,
     $DX,TAU(100),C(4,100),XM,Z,RES,SSXINT,Y2,H,S,ESTRUC,Y1
      INTEGER I
      SAVE DC,NGAM,C,TAU
      EXTERNAL FBRBM
C
      IF(INIT.NE.0) THEN
C       Compute delta function contribution
        Y=UPSLON
        XLMM=SIGZ
        XLMM = 2*SQRT(3.)*XLMM
        XL = XLMM*1.E12/.197327
        GAM = EB/5.11E-4
        RE = 1./(137.*5.11E-4)
        XKAPPA = 2./(3.*Y)
        NUCL = 2.5*Y/(SQRT(3.)*137.**2*GAM*RE)
        NUGAM = NUCL/SQRT(1.+Y**.6666666)
        NGAM=.5*NUGAM*XL
        DC = (1.-EXP(-NGAM))/NGAM
        SSFEL=0.
C       No initialization needed if >.995 included in delta peak
        IF(DC.GT..995) RETURN

C  ***  Computation of 'knots'   ***

        DX = .05
        DO 100 I = 1, 19                              
 100    TAU(I) = FLOAT(I-1)*DX     
        DO 110 I = 1, 9                               
 110    TAU(19+I) = .9 + FLOAT(I)*1.E-2
        DO 120 I = 1, 5                               
 120    TAU(28+I) = .99 + FLOAT(I)*1.E-3
        DO 121 I = 1, 12                              
 121    TAU(33+I) = .995 + FLOAT(I)*2.5E-4
        DO 130 I = 1, 20                              
 130    TAU(45+I) = .998 + FLOAT(I)*1.E-4

C   ***  Computation of corresponding y-values (electron densities)  ***

        XM = TAU(65)
        DO 140 I = 1,65
        Z = TAU(I)
        XMIN = Z
        RES=SSXINT(Z,FBRBM,XM)
140     C(1,I) = RES +DC*ESTRUC(Z,QSQBM)

C   ***  Computation of derivative at zero   ***

        Z = 1.E-5
        XMIN = Z
        RES=SSXINT(Z,FBRBM,XM)
        Y1 = RES + DC*ESTRUC(Z,QSQBM)
        Z = 1.E-4
        XMIN = Z
        RES=SSXINT(Z,FBRBM,XM)
        Y2 = RES + DC*ESTRUC(Z,QSQBM)
        C(1,2) = (Y2-Y1)/(1.E-4 - 1.E-5)
147     CALL SPLINE(TAU,C,65,1,0)
        RETURN
      ENDIF
      IF(X.GT..999999) THEN
        Z = .999999
      ELSE
        Z = X
      ENDIF
      DC = (1.-EXP(-NGAM))/NGAM
      IF(DC.GT..995) THEN
        SSFEL = DC*ESTRUC(Z,QSQBM)
        RETURN
      ENDIF

      DO 2 I = 1, 64
    2 IF(Z.LT.TAU(I+1)) GOTO 3
    3 H = Z - TAU(I)
      S = C(1,I) + H * ( C(2,I) + H*(C(3,I)+H*C(4,I)) )
      SSFEL = S
      RETURN
      END
