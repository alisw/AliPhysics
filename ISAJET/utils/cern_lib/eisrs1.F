*
* $Id$
*
* $Log$
* Revision 1.9  2000/07/25 14:53:05  mclareni
* Version 7.51 from author
*
*
*#include "sys/CERNLIB_machine.h"
#include "isajet/pilot.h"
#if defined(CERNLIB_NOCERN)
      SUBROUTINE EISRS1(NM,N,AR,WR,ZR,IERR,WORK)
C     ALL EIGENVALUES AND CORRESPONDING EIGENVECTORS OF A REAL
C     SYMMETRIC MATRIX
C     FROM CERN PROGRAM LIBRARY
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      INTEGER NM,N,IERR
      REAL AR(NM,NM),WR(N),ZR(NM,NM),WORK(1)
      CALL TRED2(NM,N,AR,WR,WORK,ZR)
      CALL TQL2(NM,N,WR,WORK,ZR,IERR)
      RETURN
      END
#endif
