*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:20  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
*
      SUBROUTINE GDBOMB(ICPOIN,ISHAPE)
*
********************************************************************************
*
*  This subroutine allows exploded view of detectors !
*
*  GBOOM > 0 => EXPLOSION
*
*  GBOOM = 0 => RESET
*
*  GBOOM < 0 => IMPLOSION
*
********************************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcgobj.inc"
*SG
#include "geant321/gcmutr.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gchiln.inc"
#include "geant321/gcspee.inc"
      DIMENSION RRMIN(3),RRMAX(3)
      SAVE SEXP1, SEXP2, SEXP3
C.    ------------------------------------------------------------------
*
*
      IF(GBOOM.NE.0)THEN
         IF(ISUBLI.EQ.1) THEN
            CALL CGSHIF(SEXP1,SEXP2,SEXP3,Q(ICPOIN))
            ITSTCU=1
         ELSE
            CALL CGMNMX(Q(ICPOIN),RRMIN,RRMAX)
            EXP1=(RRMAX(1)+RRMIN(1))*.5
            EXP2=(RRMAX(2)+RRMIN(2))*.5
            EXP3=(RRMAX(3)+RRMIN(3))*.5
            AEXP1=ABS(EXP1)
            AEXP2=ABS(EXP2)
            AEXP3=ABS(EXP3)
            AA1=MAX(ABS(RRMAX(1)),ABS(RRMIN(1)))
            AA2=MAX(ABS(RRMAX(2)),ABS(RRMIN(2)))
            AA3=MAX(ABS(RRMAX(3)),ABS(RRMIN(3)))
            IF((AEXP1/AA1).LE.0.05.AND.(AEXP2/AA2).LE.0.05.AND.
     +                 (AEXP3/AA3).LE.0.05)THEN
               ITSTCU=0
               IF(IPORLI.EQ.1)LPASS=1
            ELSE
               ITSTCU=1
               IF(ISHAPE.EQ.11.OR.ISHAPE.EQ.12)LPASS=1
               SEXP1=GBOOM*EXP1
               SEXP2=GBOOM*EXP2
               SEXP3=GBOOM*EXP3
               CALL CGSHIF(SEXP1,SEXP2,SEXP3,Q(ICPOIN))
            ENDIF
         ENDIF
      ENDIF
*
 999  END
 
