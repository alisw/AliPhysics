#include "isajet/pilot.h"
      SUBROUTINE DECPS2(IP,NADD,PGEN,PREST,BETA,GAMMA)
C
C          Carry out decays using masses from DECPS1 or special matrix
C          elements.
C          Auxiliary routine for DECAY.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
C
#include "isajet/itapes.inc"
#include "isajet/partcl.inc"
#include "isajet/const.inc"
C
      INTEGER IP,NADD
      REAL PGEN(5,5),PREST(4,6)
      REAL PCM,DBLPCM,RANF
      REAL U(3),BETA(3)
      REAL QCM,PHI,A,B,C,GAMMA,BP
      INTEGER I,J,NADD1,II,K,K1
C
C          Function definitions.
C
#if defined(CERNLIB_SINGLE)
      PCM(A,B,C)=SQRT((A-B-C)*(A+B+C)*(A-B+C)*(A+B-C))/(2.*A)
#endif
#if defined(CERNLIB_DOUBLE)
      PCM(A,B,C)=DBLPCM(A,B,C)
#endif
C
C          Carry out two-body decays in PGEN frames.
C
      NADD1=NADD-1
100   CONTINUE
      DO 110 I=1,NADD1
        QCM=PCM(PGEN(5,I),PGEN(5,I+1),PPTCL(5,NPTCL+I))
        U(3)=2.*RANF()-1.
        PHI=2.*PI*RANF()
        U(1)=SQRT(1.-U(3)**2)*COS(PHI)
        U(2)=SQRT(1.-U(3)**2)*SIN(PHI)
        DO 120 J=1,3
          PPTCL(J,NPTCL+I)=QCM*U(J)
          PGEN(J,I+1)=-PPTCL(J,NPTCL+I)
120     CONTINUE
        PPTCL(4,NPTCL+I)=SQRT(QCM**2+PPTCL(5,NPTCL+I)**2)
        PGEN(4,I+1)=SQRT(QCM**2+PGEN(5,I+1)**2)
110   CONTINUE
C
      DO 130 J=1,4
        PPTCL(J,NPTCL+NADD)=PGEN(J,NADD)
130   CONTINUE
C
C          Boost PGEN frames to lab frame, saving momenta in rest frame.
C
      DO 200 II=1,NADD1
        I=NADD-II
        DO 210 J=1,3
          BETA(J)=PGEN(J,I)/PGEN(4,I)
210     CONTINUE
        GAMMA=PGEN(4,I)/PGEN(5,I)
        DO 220 K=I,NADD
          K1=NPTCL+K
          BP=BETA(1)*PPTCL(1,K1)+BETA(2)*PPTCL(2,K1)+BETA(3)*PPTCL(3,K1)
          DO 230 J=1,3
            PREST(J,K)=PPTCL(J,K1)
            PPTCL(J,K1)=PPTCL(J,K1)+GAMMA*BETA(J)*(PPTCL(4,K1)
     $      +BP*GAMMA/(GAMMA+1.))
230       CONTINUE
          PREST(4,K)=PPTCL(4,K1)
          PPTCL(4,K1)=GAMMA*(PPTCL(4,K1)+BP)
220     CONTINUE
200   CONTINUE
C
      RETURN
      END
