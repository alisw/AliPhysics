#include "isajet/pilot.h"
        SUBROUTINE SSWZBF
C-----------------------------------------------------------------------
C       This subroutine gives chargino (wi) and neutralino (zi)
C       decays according to Baer, et al.
C       Valid for all scalar masses (functions in double precision)
C       Includes Higgs sector radiative corrections (Aug. 31)
C
C       Auxiliary functions are called SSWxyi, SSZxyi, where normally
C       x indicates the SUSY particle, y the SM particle(s), and i is
C       a counter.
C
C       Ver. 7.07: For w_i -> z_j or z_i -> w_j decays, require that
C                  decay be allowed by a factor FUDGE = 1.0
C       Ver. 7.28: Calculate full large tan(beta) decays.
C                  Calculate tau polarizations.
C       Ver. 7.33: Add GMSB modes.
C       Ver. 7.41: Add decay matrix elements.
C                  Split real work into sszibf and sswfbf.
C
C       Baer's GAUGBF
C
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/ssmode.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstype.inc"
#include "isajet/sstmp.inc"
#include "isajet/sspols.inc"
C
      INTEGER IZ
C
C          Initialize polarizations to zero
C
      DO 150 IZ=1,4
        PTAU1(IZ)=0
        PTAU2(IZ)=0
150   CONTINUE
      PTAUZZ=0
      PTAUWZ=0
C-----------------------------------------------------------------------
C          Generate Neutralino zi Branching Fractions
C-----------------------------------------------------------------------
      CALL SSZIBF
C-----------------------------------------------------------------------
C          Generate Chargino Branching Fractions
C-----------------------------------------------------------------------
      CALL SSWIBF
C
C          Set more neutralino polarizations
C
      IF (ABS(AMZISS(2)).GT.(AMTAU+AML1SS)) PTAUZ2(1)=PTAU1(2)
      IF (ABS(AMZISS(2)).GT.(AMTAU+AML2SS)) PTAUZ2(2)=PTAU2(2)
      IF (ABS(AMZISS(3)).GT.(AMTAU+AML1SS)) PTAUZ3(1)=PTAU1(3)
      IF (ABS(AMZISS(3)).GT.(AMTAU+AML2SS)) PTAUZ3(2)=PTAU2(3)
      IF (ABS(AMZISS(4)).GT.(AMTAU+AML1SS)) PTAUZ4(1)=PTAU1(4)
      IF (ABS(AMZISS(4)).GT.(AMTAU+AML2SS)) PTAUZ4(2)=PTAU2(4)
C
      RETURN
      END
