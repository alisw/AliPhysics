*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE RANMAR(RVEC,LENV)
#if defined(CERNLIB_QMCRY)
CDIR$ STACK
#endif
C
C CERN PROGLIB# V113    RANMAR          .VERSION KERNFOR  4.21  890323
C ORIG. 01/03/89 FCA + FJ
C
      DIMENSION RVEC(*)
C
      COMMON/RANMA1/IJKL,NTOT,NTOT2,I97,J97,C,U(97)
      LOGICAL FIRST
      PARAMETER (TWOM24=2.**(-24),TWOM48=2.**(-48))
      PARAMETER (CD=7654321.*TWOM24,CM=16777213.*TWOM24)
      PARAMETER (CINT=362436.*TWOM24,MODCNS=1000000000)
      SAVE /RANMA1/, FIRST
      DATA FIRST/.TRUE./
C
      IF(FIRST) THEN
        IJKL = 54217137
        NTOT = 0
        NTOT2 = 0
        GO TO 70
      ENDIF
C
   80 CONTINUE
      DO 100 IVEC= 1, LENV
        UNI = U(I97)-U(J97)
        IF (UNI .LT. 0.) UNI=UNI+1.
        U(I97) = UNI
        I97 = I97-1
        IF (I97 .EQ. 0)  I97=97
        J97 = J97-1
        IF (J97 .EQ. 0)  J97=97
        C = C - CD
        IF (C .LT. 0.)   C=C+CM
        UNI = UNI-C
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
      NTOT = NTOT + LENV
      IF (NTOT .GE. MODCNS)  THEN
        NTOT2 = NTOT2 + 1
        NTOT  = NTOT - MODCNS
      ENDIF
      RETURN
      ENTRY RMARIN(IJKLIN,NTOTIN,NTO2IN)
C
      FIRST = .FALSE.
      IJKL  = IJKLIN
      NTOT  = NTOTIN
      NTOT2 = NTO2IN
C
   70 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
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
  20    CONTINUE
        U(II) = S
  30  CONTINUE
      C   = CINT
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      NITER = MODCNS
      DO 50 LOOP2= 1, NTOT2+1
        IF(LOOP2.GT.NTOT2) NITER=NTOT
        DO 40 IDUM = 1, NITER
          UNI = U(I97)-U(J97)
          IF (UNI .LT. 0.) UNI=UNI+1.
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. 0.)   C=C+CM
   40   CONTINUE
   50 CONTINUE
      NTOT  = 0
      NTOT2 = 0
      IF(FIRST) THEN
        FIRST = .FALSE.
        GO TO 80
      ENDIF
      RETURN
      ENTRY RMARUT(IJKLUT,NTOTUT,NTO2UT)
C
      NTOTUT = NTOT
      NTO2UT = NTOT2
      IJKLUT = IJKL
C
      END
