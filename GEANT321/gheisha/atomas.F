*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      REAL FUNCTION ATOMAS(A,Z)
C
C *** DETERMINATION OF THE ATOMIC MASS ***
C *** NVE 19-MAY-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (02-DEC-1986)
C
#include "geant321/s_consts.inc"
C
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION AA,ZZ,MASS
#endif
C
C --- GET ATOMIC (= ELECTRONS INCL.) MASSES (IN MEV) FROM RMASS ARRAY ---
C --- ELECTRON ---
      RMEL=RMASS(4)*1000.
C --- PROTON ---
      RMP=RMASS(14)*1000.
C --- NEUTRON ---
      RMN=RMASS(16)*1000.
C --- DEUTERON ---
      RMD=RMASS(30)*1000.+RMEL
C --- ALPHA ---
      RMA=RMASS(32)*1000.+2.*RMEL
C
      ATOMAS = 0.
      AA = A * 1.D0
      ZZ = Z * 1.D0
      IA = IFIX(A + 0.5)
      IF(IA.LT.1) RETURN
      IZ = IFIX(Z + 0.5)
      IF(IZ.LT.0) RETURN
      IF(IZ.GT.IA) RETURN
      IF(IA.GT.4) GOTO 50
      MASS=0.D0
      GOTO (10,20,50,40),IA
   10 IF(IZ.EQ.0) MASS=RMN
      IF(IZ.EQ.1) MASS=RMP+RMEL
      GOTO 60
   20 IF(IZ.NE.1) GOTO 50
      MASS=RMD
      GOTO 60
   40 IF(IZ.NE.2) GOTO 50
      MASS=RMA
      GOTO 60
   50 MASS=(AA-ZZ)*RMN + ZZ*RMP +ZZ*RMEL - 15.67*AA
     *     + 17.23*(AA**0.6666667) + 93.15*((AA/2.-ZZ)**2)/AA
     *     +0.6984523*ZZ**2/(AA**0.3333333)
      IPP=MOD(IA-IZ,2)
      IZZ=MOD(IZ,2)
      IF(IPP.NE.IZZ) GOTO 60
      MASS = MASS + (IPP+IZZ- 1)*12.00*(AA**(-0.5))
   60 ATOMAS = MASS*0.001
      RETURN
      END
