*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE NAMES(NXARG,NYARG)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PURPOSE                                                            C
C     TO DEFINE WHICH LABELS SHALL BE PRINTED AT THE AXES BY           C
C     SUBSEQUENT CALLS TO THE SUBROUTINE CONT.                         C
C   USAGE                                                              C
C     IF YOU FOR EXAMPLE WANT THE WORD KSI TO BE PRINTED AT THE        C
C     HORIZONTAL AXIS AND THE WORD ETA TO BE PRINTED AT THE VERTICAL   C
C     AXIS, EXECUTE THE STATEMENT                                      C
C     CALL NAMES(NXARG,NYARG)                                          C
C     BEFORE THE CALL TO CONT.WHERE -                                  C
C     DATA NXARG,NYARG/6HKSI   ,6HETA   /
C   PARAMETERS                                                         C
C     NXARG  -  THE NAME TO BE PRINTED AT THE HORIZONTAL AXIS.         C
C     NYARG  -  THE NAME TO BE PRINTED AT THE VERTICAL AXIS.           C
C   REMARKS                                                            C
C     THE DEFAULT OPTION OF BLANK NAMES IS AUTOMATICALLY PROVIDED.     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON /J509C1/XBL,YBL,XUR,YUR,NX(6),NY(6),LETT1,LETT2,LINE1,LINE2
C
      CALL J509BD
      CALL UBLOW(NXARG,NX,6)
      CALL UBLOW(NYARG,NY,6)
      RETURN
      END
