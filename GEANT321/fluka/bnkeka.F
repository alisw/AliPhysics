*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE BNKEKA.FOR
*COPY BNKEKA
*                                                                      *
*=== bnkeka ===========================================================*
*                                                                      *
      FUNCTION BNKEKA(IT,TO,AMSS,SQAMSS )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C*** CASCADE PARAMETER CALCULATION
      DIMENSION A(2),B(2),C(2),D(2)
      SAVE A,B,C,D
      DATA A/1.D0,1.3D0/
      DATA B/.15D0,.225D0/
      DATA C /1.35D0,2.28D0/
      DATA D/.3D0,.4D0/
      IF(TO-0.1D0) 1,1,2
    1 CONTINUE
      BNKEKA=SQAMSS*B(IT)*0.1D0
      RETURN
    2 CONTINUE
      IF (TO-10.D0) 3,4,4
    3 CONTINUE
      BNKEKA=0.1D0*SQAMSS*(0.5D0+A(IT)*(1.D0+LOG10(TO))**2)*D(IT)
      RETURN
    4 CONTINUE
      BNKEKA=SQAMSS*0.1D0*C(IT)
      RETURN
      END
