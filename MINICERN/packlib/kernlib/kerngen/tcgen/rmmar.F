*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE RMMAR(RVEC,LENV,ISEQ)
#if defined(CERNLIB_QMCRY)
CDIR$ STACK
#endif
C
C CERN PROGLIB# V113    RMMAR           .VERSION KERNFOR  1.0
C ORIG. 01/03/89 FCA + FJ
C
      DIMENSION RVEC(*)
C
      COMMON/RANMA2/JSEQ,IU(103)
      DIMENSION I97(0:1),J97(0:1),C(0:1),NTOT(0:1)
      DIMENSION NTOT2(0:1),IJKL(0:1),U(97)
      EQUIVALENCE (IJKL(0),IU(1)),(NTOT(0),IU(2)),(NTOT2(0),IU(3))
      EQUIVALENCE (U(1),IU(4)),(C(0),IU(101)),(I97(0),IU(102))
      EQUIVALENCE (J97(0),IU(103))
      PARAMETER (TWOM24=2.**(-24),TWOM48=2.**(-48))
      PARAMETER (CD=7654321.*TWOM24,CM=16777213.*TWOM24)
      PARAMETER (MODCNS=1000000000)
      SAVE /RANMA2/
C
      IF(ISEQ.GT.0) JSEQ = ISEQ
      IBASE = (JSEQ-1)*103
C
      DO 100 IVEC= 1, LENV
        UNI = U(IBASE+I97(IBASE))-U(IBASE+J97(IBASE))
        IF (UNI .LT. 0.) UNI=UNI+1.
        U(IBASE+I97(IBASE)) = UNI
        I97(IBASE) = I97(IBASE)-1
        IF (I97(IBASE) .EQ. 0)  I97(IBASE)=97
        J97(IBASE) = J97(IBASE)-1
        IF (J97(IBASE) .EQ. 0)  J97(IBASE)=97
        C(IBASE) = C(IBASE) - CD
        IF (C(IBASE) .LT. 0.)   C(IBASE)=C(IBASE)+CM
        UNI = UNI-C(IBASE)
        IF (UNI .LT. 0.) UNI=UNI+1.
C
C   Replace exact zeroes by uniform distr. *2**-24
C
        IF (UNI .EQ. 0.)  THEN
          UNI = TWOM24*U(2)
C
C   An exact zero here is very unlikely, but let's be safe.
C
          IF (UNI .EQ. 0.) UNI= TWOM48
        ENDIF
        RVEC(IVEC) = UNI
  100 CONTINUE
C
      NTOT(IBASE) = NTOT(IBASE) + LENV
      IF (NTOT(IBASE) .GE. MODCNS)  THEN
        NTOT2(IBASE) = NTOT2(IBASE) + 1
        NTOT(IBASE)  = NTOT(IBASE) - MODCNS
      ENDIF
      END
