*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE RMMAQ(ISEED,ISEQ,CHOPT)
C
C CERN PROGLIB# V113    RMMAQ           .VERSION KERNFOR  1.0
C ORIG. 01/03/89 FCA + FJ
C
      COMMON/RANMA2/JSEQ,IU(103)
      DIMENSION   I97(0:1),J97(0:1),C(0:1),NTOT(0:1)
      DIMENSION   NTOT2(0:1),IJKL(0:1),U(97)
      EQUIVALENCE (IJKL(0),IU(1)),(NTOT(0),IU(2)),(NTOT2(0),IU(3))
      EQUIVALENCE (U(1),IU(4)),(C(0),IU(101)),(I97(0),IU(102))
      EQUIVALENCE (J97(0),IU(103))
      DIMENSION   ISEED(*),UU(97)
      CHARACTER   CHOPT*(*), CCHOPT*12
      PARAMETER (TWOM24=2.**(-24),TWOM48=2.**(-48))
      PARAMETER (CD=7654321.*TWOM24,CM=16777213.*TWOM24)
      PARAMETER (CINT=362436.*TWOM24,MODCNS=1000000000)
      SAVE /RANMA2/
C
      CCHOPT = CHOPT
      IF(CCHOPT.EQ.' ') THEN
        ISEED(1) = 54217137
        ISEED(2) = 0
        ISEED(3) = 0
        CCHOPT = 'S'
        JSEQ = 1
      ENDIF
C
      IF(INDEX(CCHOPT,'S').NE.0) THEN
        IF(ISEQ.GT.0) JSEQ=ISEQ
        IBASE = (JSEQ-1)*103
        IF(INDEX(CCHOPT,'V').NE.0) THEN
          DO 10 JJ=1, 103
            IU(IBASE+JJ) = ISEED(JJ)
  10      CONTINUE
        ELSE
          IJKL(IBASE)  = ISEED(1)
          NTOT(IBASE)  = ISEED(2)
          NTOT2(IBASE) = ISEED(3)
          IJ = IJKL(IBASE)/30082
          KL = IJKL(IBASE) - 30082*IJ
          I = MOD(IJ/177, 177) + 2
          J = MOD(IJ, 177)     + 2
          K = MOD(KL/169, 178) + 1
          L = MOD(KL, 169)
          DO 30 II= 1, 97
            S = 0.
            T = .5
            DO 20 JJ= 1, 24
              M = MOD(MOD(I*J,179)*K, 179)
              I = J
              J = K
              K = M
              L = MOD(53*L+1, 169)
              IF (MOD(L*M,64) .GE. 32)  S = S+T
              T = 0.5*T
  20        CONTINUE
            UU(II) = S
  30      CONTINUE
          CC   = CINT
          II97 = 97
          IJ97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
          NITER = MODCNS
          DO 50 LOOP2= 1, NTOT2(IBASE)+1
            IF(LOOP2.GT.NTOT2(IBASE)) NITER=NTOT(IBASE)
            DO 40 IDUM = 1, NITER
              UNI = UU(II97)-UU(IJ97)
              IF (UNI .LT. 0.) UNI=UNI+1.
              UU(II97) = UNI
              II97 = II97-1
              IF (II97 .EQ. 0)  II97=97
              IJ97 = IJ97-1
              IF (IJ97 .EQ. 0)  IJ97=97
              CC = CC - CD
              IF (CC .LT. 0.)   CC=CC+CM
   40       CONTINUE
   50     CONTINUE
          I97(IBASE) = II97
          J97(IBASE) = IJ97
          C(IBASE)   = CC
          DO 60 JJ=1, 97
            U(IBASE+JJ) = UU(JJ)
  60      CONTINUE
        ENDIF
      ELSE IF(INDEX(CCHOPT,'R').NE.0) THEN
        IF(ISEQ.GT.0) THEN
          JSEQ=ISEQ
        ELSE
          ISEQ=JSEQ
        ENDIF
        IBASE = (JSEQ-1)*103
        IF(INDEX(CCHOPT,'V').NE.0) THEN
          NCOPY = 103
        ELSE
          NCOPY = 3
        ENDIF
        DO 70 JJ=1,NCOPY
          ISEED(JJ) = IU(IBASE+JJ)
  70    CONTINUE
      ENDIF
C
      END
