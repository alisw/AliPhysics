#include "isajet/pilot.h"
C----------------------------------------------------------------------
      FUNCTION SUALFS(QSQ,ALAM4,TMASS,LOOP)
C----------------------------------------------------------------------
C
C       This function returns the 1, 2, or 3-loop value of alpha_s
C       Input:
C               QSQ     = Q**2 (real)
C               ALAM4   = Lambda for 4 active quark flavors (real)
C               TMASS   = top quark mass to determine lambda-6 (real)
C               LOOP    = number of loops for alpha_s (= 1, 2 or 3)
C       Parametrization of the strong coupling constant according to
C       LOOP = 1, 2 : from the book;
C       LOOP = 3:     W. J. Marciano, Phys. Rev. D 29 (1984) 580.
C       Note : threshold at 2*Mq
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
      REAL SUALFS, QSQ, ALAM4, TMASS
      INTEGER LOOP
      REAL PI, BMASS
      REAL ANF, ALAM, ALAMSQ, ALAM5, T, TT, B0, B1, B2, X, ALPHAS
      DATA PI/3.1415927/, BMASS/5.0/
C
      IF (QSQ .LT. 4.0*BMASS**2) THEN
        ANF   = 4.0
        ALAM  = ALAM4
      ELSE IF (QSQ .LT. 4.0*TMASS**2) THEN
        ANF   = 5.0
        ALAM  = ALAM4*(ALAM4/(2.0*BMASS))**(2.0/23.0)
     1           *(ALOG(4.0*BMASS**2/ALAM4**2))**(-963.0/13225.0)
      ELSE
        ANF   = 6.0
        ALAM5 = ALAM4*(ALAM4/(2.0*BMASS))**(2.0/23.0)
     1           *(ALOG(4.0*BMASS**2/ALAM4**2))**(-963.0/13225.0)
        ALAM  = ALAM5*(ALAM5/(2.0*TMASS))**(2.0/21.0)
     1           *(ALOG(4.0*TMASS**2/ALAM5**2))**(-107.0/1127.0)
      END IF
      B0       = 11.0-2.0/3.0*ANF
      ALAMSQ   = ALAM**2
      T        = ALOG(QSQ/ALAMSQ)
      IF (T .LE. 1.0) T = ALOG(4.0/ALAMSQ)
      ALPHAS   = 4*PI/B0/T
      IF (LOOP .EQ. 1) THEN
        SUALFS = ALPHAS
      ELSE IF (LOOP .EQ. 2) THEN
        B1 = 102.0-38.0/3.0*ANF
        X  = B1/(B0**2*T)
        TT = ALOG(T)
        SUALFS = ALPHAS*(1.0-X*TT)
      ELSE IF (LOOP .EQ. 3) THEN
        B1 = 102.0-38.0/3.0*ANF
        B2 = 0.5*(2857.0-5033.0/9.0*ANF+325.0/27.0*ANF**2)
        X  = B1/(B0**2*T)
        TT = ALOG(T)
        SUALFS = ALPHAS*(1.0-X*TT+X**2*((TT-0.5)**2
     $      +B2*B0/B1**2-1.25))
      ELSE
        WRITE(LOUT,*) ' WRONG LOOP NUMBER IN ALPHA-S EVALUATION!'
        STOP 99
      END IF
C
      RETURN
      END
