#include "isajet/pilot.h"
      SUBROUTINE DECPS1(IP,NADD,PGEN)
C
C          Generate masses for uniform NADD-body phase space in DECPS2.
C          Auxiliary routine for DECAY.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
C
#include "isajet/itapes.inc"
#include "isajet/partcl.inc"
C
      INTEGER IP,NADD
      REAL PGEN(5,5)
      REAL REDUCE(5),RND(5)
      REAL RANF,PCM,DBLPCM
      REAL WTMAX,SUM1,SUM2,SUM,RNEW,WT,A,B,C
      INTEGER I,NADD1,J,I1,JJ1,JSAVE
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
      DATA REDUCE/1.,1.,2.,5.,15./
C
C          Calculate maximum phase-space weight.
C
      IF(NADD.LE.2) RETURN
      NADD1=NADD-1
      WTMAX=1./REDUCE(NADD)
      SUM=0
      DO 100 I=1,NADD
        SUM=SUM+PPTCL(5,NPTCL+I)
100   CONTINUE
      SUM1=PGEN(5,1)
      SUM2=SUM-PPTCL(5,NPTCL+1)
      DO 110 I=1,NADD1
        WTMAX=WTMAX*PCM(SUM1,SUM2,PPTCL(5,NPTCL+I))
        SUM1=SUM1-PPTCL(5,NPTCL+I)
        SUM2=SUM2-PPTCL(5,NPTCL+I+1)
110   CONTINUE
C
C          Generate masses for uniform NADD-body phase space.
C
200   CONTINUE
      RND(1)=1.
      DO 210 I=2,NADD1
        RNEW=RANF()
        I1=I-1
        DO 220 JJ1=1,I1
          J=I-JJ1
          JSAVE=J+1
          IF(RNEW.LE.RND(J)) GO TO 210
          RND(JSAVE)=RND(J)
220     CONTINUE
210   RND(JSAVE)=RNEW
      RND(NADD)=0.
      WT=1.
      SUM1=SUM
      DO 230 I=2,NADD
        SUM1=SUM1-PPTCL(5,NPTCL+I-1)
        PGEN(5,I)=SUM1+RND(I)*(PGEN(5,1)-SUM)
        IF(PGEN(5,I-1).LE.PGEN(5,I)+PPTCL(5,NPTCL+I-1)) GO TO 200
        WT=WT*PCM(PGEN(5,I-1),PGEN(5,I),PPTCL(5,NPTCL+I-1))
230   CONTINUE
      IF(WT.LT.RANF()*WTMAX) GO TO 200
C
      RETURN
      END
