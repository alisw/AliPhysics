*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      FUNCTION GEFCH2(Z1,Z2,TK)
C.    *******************************************************************
C.    *                                                                 *
C.    *  Calculates  gamma**2=( eff.charge/z1)**2 for heavy ions        *
C.    *   (see e.g. F.Hubert et al. Atomic Data and Nuclear Data        *
C.    *             Tables,Vol.40,1.(1990) )                            *
C.    *                                                                 *
C.    *  The stopping power S for a given combination of stopping       *
C.    *  medium and heavy ion can be calculated according to            *
C.    *  the scaling law                                                *
C.    *                                                                 *
C.    *            S = (gamma*Z1)**2*Sref/Zref**2                       *
C.    *                                                                 *
C.    *    where Z1 the ion charge , gamma*Z1 is the effective charge,  *
C.    *          Sref is the stopping power of the same medium          *
C.    *          for a reference ion of the same velocity and           *
C.    *          of effective charge Zref.                              *
C.    *                                                                 *
C.    *  Input: Z1 atomic number of the ion                             *
C.    *         Z2 atomic number of the medium                          *
C.    *         TK E/A in GeV (ion kin.energy/atomic weight)            *
C.    *                                                                 *
C.    *******************************************************************
C.
      ALZ2=LOG(Z2)
      TKMEV=1000.*TK
C
      X1=1.164+0.2319*EXP(-0.004302*Z2)+1.658*EXP(-0.05170*Z1)
      X2=8.144+0.09876*ALZ2
      X3=0.3140+0.01072*ALZ2
      X4=0.5218+0.02521*ALZ2
      EE=MAX(-50.,-X2*TKMEV**X3/Z1**X4)
C     .........................
      GA=1.-X1*EXP(EE)
C
      GEFCH2=GA**2
      END
