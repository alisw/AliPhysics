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
*$ CREATE BEKEKA.FOR
*COPY BEKEKA
*                                                                      *
*=== bekeka ===========================================================*
*                                                                      *
      FUNCTION BEKEKA(IX,TO,AMSS,SQAMSS)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C*** CALCULATION OF AVERAGE CASCADE AND EXCITATION ENERGY
C*****IX=1 EEV    2EPK    3 ENK     4 EEX=EPK+EEV    5 EEXT=EEX+ENK
      GO TO (1,2,3,1,1),IX
    1 CONTINUE
      IF(TO-0.1D0) 11,11,12
   11 CONTINUE
      AA=0.001D0*SQAMSS
      GO TO 19
   12 CONTINUE
      APAR=0.035D0
      BPAR=3.D0
      CPAR=0.1D0
      AA=CPAR*SQAMSS*(0.01D0+APAR*(BPAR+LOG10(TO))**2)
   19 CONTINUE
      IF(IX.GT.3) GO TO 2
C    ENERGYADJUST
      ENADJ=1.D0
      AA=AA*ENADJ
      BEKEKA=AA
      RETURN
    2 CONTINUE
      AN=BNKEKA(1,TO,AMSS,SQAMSS)
      A =BKEKA (1,TO,AMSS)
      EXTOA=0.D0
      IF(TO.LT.5.D0*A) EXTOA=EXP(-TO/A)
      TPKAV=A*(1.D0-(TO/A+1.D0)*EXTOA)/(1.D0-EXTOA)
      BB=TPKAV*AN
      IF(IX.EQ.4) GO TO 4
      IF(IX.EQ.5) GO TO 3
      BEKEKA=BB
      RETURN
    4 CONTINUE
      BEKEKA=AA+BB
      RETURN
    3 CONTINUE
      AN=BNKEKA(2,TO,AMSS,SQAMSS)
      A =BKEKA (2,TO,AMSS)
      EXTOA=0.D0
      IF(TO.LT.5.D0*A) EXTOA=EXP(-TO/A)
      TNKAV=A*(1.D0-(TO/A+1.D0)*EXTOA)/(1.D0-EXTOA)
      CC=TNKAV*AN
      IF(IX.EQ.5) GO TO 5
      BEKEKA=CC
      RETURN
    5 CONTINUE
      BEKEKA=AA+BB+CC
      RETURN
      END
