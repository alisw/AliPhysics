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
*$ CREATE PARJET.FOR
*COPY PARJET
*
*=== parjet ===========================================================*
*
      SUBROUTINE PARJET(IHAD,I)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
#include "geant321/finpar.inc"
#include "geant321/part.inc"
#include "geant321/inpdat2.inc"
      IHAD=1
      NREF(1)=I
      PXF(1)=0.D0
      PYF(1)=0.D0
      PZF(1)=0.D0
      HEP(1)=AM(I)
      AMF(1)=AM(I)
      ICHF(1)=ICH(I)
      IBARF(1)=IBAR(I)
      ANF(1)=ANAME(I)
      RETURN
      END
