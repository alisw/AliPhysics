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
*$ CREATE BKEKA.FOR
*COPY BKEKA
*                                                                      *
*=== bkeka ============================================================*
*                                                                      *
      FUNCTION BKEKA(IT,TO,AMSS)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C*** NUCLEAR CASCADE NUCLEON EMISSION PARAMETERS
      DIMENSION A(2),B(2)
      SAVE A,B
      DATA A/0.11D0,0.1D0/
      DATA B/0.21D0,0.20D0/
      IF (TO-10.D0) 1,2,2
    1 CONTINUE
C    ENERGYADJUST
      ENADJ=1.D0
      BKEKA=ENADJ*(1.D0-0.001D0*AMSS)*(A(IT)+0.01D0*TO)
      RETURN
    2 CONTINUE
      BKEKA=B(IT)*(1.D0-0.001D0*AMSS)
      RETURN
      END
