#ifndef CERNLIB_ISAJET_SSPOLS_INC
#define CERNLIB_ISAJET_SSPOLS_INC
*
*
* sspols.inc
*
C          Polarizations in SUSY decays
C          PTAUj(i)  = P_tau for tauj -> ziss tau
C          PTAUZi(j) = P_tau for ziss -> tauj tau
C          PTAUZZ    = P_tau for z2ss -> z1ss tau tau
C          PTAUWZ    = P_tau for w1ss -> z1ss tau nutau
      COMMON/SSPOLS/PTAU1(4),PTAU2(4),PTAUZ2(2),PTAUZ3(2),PTAUZ4(2),
     $PTAUZZ,PTAUWZ
      SAVE /SSPOLS/
      REAL PTAU1,PTAU2,PTAUZ2,PTAUZ3,PTAUZ4,PTAUZZ,PTAUWZ

#endif
