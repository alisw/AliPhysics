#include "isajet/pilot.h"
      SUBROUTINE SSHIBF
C-----------------------------------------------------------------------
C
C     This subroutine calculates the decay widths for decays of the 
C     Higgs scalars present in the minimal SUSY model.
C
C     NOTE: Decays into sfermions are not yet incorporated.
C
C     Standard model parameters are hard wired in  SSMSSM. To get
C     the 1987-8 values corresponding to the Gunion et al. papers
C     (Intl. J. Mod. Phys. 2(4):1035; Nucl. Phys. B307:445) you must
C     change
C          ALFA3 = 0.12  --> 0.136
C          AMW   = 80.0  --> 81.3
C          AMZ   = 91.17 --> 92.7
C
C     2/9/91:
C     I've modified the program slightly.  The ALPHA3 = 0.12 value
C     above is the recent empirical value from LEP. Using the equation
C     from page 220 in Barger and Phillips yields ALPHA3 = 0.136.
C
C     10/1/92:
C     Now includes vertex corrections for triple Higgs couplings.
C     (See Kunszt and Zwirner, CERN-TH.6150/91 for all but HH-HC-HC
C     correction which is in our Higgs --> SUSY paper: Baer et al. 
C     FSU-HEP-920630 or UH-511-749-92)
C
C     Bisset's HIGSBF
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sstype.inc"
C
C          Decays into fermions
      CALL SSHFF
C          Loop decays into photons and gluons
      CALL SSHGM
      CALL SSHGL
C          Decays into WW(*), ZZ(*)
      CALL SSHWW
C          Decays into neutralinos and charginos
      CALL SSHNN
      CALL SSHCC
C          Decays into other Higgs bosons
      CALL SSHHX
C          Decays to sfermions
      CALL SSHSF
C          Normalize branching ratios
C
      CALL SSNORM(ISHL)
      CALL SSNORM(ISHH)
      CALL SSNORM(ISHA)
      CALL SSNORM(ISHC)
C
      RETURN
      END
