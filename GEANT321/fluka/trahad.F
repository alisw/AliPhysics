*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE TRAHAD.FOR
*COPY TRAHAD
*
*=== trahad ===========================================================*
*
      SUBROUTINE TRAHAD(N,NAUX,GAM,BGX,BGY,BGZ)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C
C     LORENTZ TRANSFORMATION OF  N PARTICLES IN  FINPAR  TO BE
C     STORED IN  AUXPAR  STARTING AT NAUX
C
C********************************************************************
C
#include "geant321/hadpar.inc"
#include "geant321/auxpar.inc"
      DO 1  I=1,N
         J = NAUX + I - 1
         AMH(J)   = AMA(I)
         ICHH(J)  = ICHA(I)
         IBARH(J) = IBARA(I)
         ANH(J)   = ANA(I)
         NREH(J)  = NREA(I)
         CALL ALTRA(GAM,BGX,BGY,BGZ,PXA(I),PYA(I),PZA(I),HEPA(I),
     &              PPA,PXH(J),PYH(J),PZH(J),HEPH(J))
   1  CONTINUE
      RETURN
      END
