*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE DECAUX.FOR
*COPY DECAUX
*
*=== decaux ===========================================================*
*
      SUBROUTINE DECAUX(NHAD,NUMU)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*  Decaux89: slight revision by A. Ferrari                             *
*----------------------------------------------------------------------*
*
#include "geant321/finpar2.inc"
#include "geant321/auxpar.inc"
C      STORE VARIABLES FROM HADPAR INTO FINPAR
      DO 1 I=1,NHAD
         PXF(I)  = PXA(I)
         PYF(I)  = PYA(I)
         PZF(I)  = PZA(I)
         HEF(I)  = HEPA(I)
         AMF(I)  = AMA(I)
         ICHF(I) = ICHA(I)
         IBARF(I)= IBARA(I)
         ANF(I)  = ANA(I)
         NREF(I) = NREA(I)
C     PRINT 3358,PXA(I),PYA(I),PZA(I),HEPA(I),AMA(I),ICHA(I),IBARA(I),
C    1ANA(I),NREA(I)
C3358 FORMAT (5E15.5,2I10,A8,I10)
    1 CONTINUE
C
      CALL FKDECA (NHAD,NUMU)
C
      DO 2 I=1,NHAD
         PXA(I)  = PXF(I)
         PYA(I)  = PYF(I)
         PZA(I)  = PZF(I)
         HEPA(I) = HEF(I)
         AMA(I)  = AMF(I)
         ICHA(I) = ICHF(I)
         IBARA(I)= IBARF(I)
         ANA(I)  = ANF(I)
         NREA(I) = NREF(I)
C     PRINT 3358,PXA(I),PYA(I),PZA(I),HEPA(I),AMA(I),ICHA(I),IBARA(I),
C    1ANA(I),NREA(I)
    2 CONTINUE
      RETURN
      END
