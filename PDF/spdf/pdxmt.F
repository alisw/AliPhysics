*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:02  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:26  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
      FUNCTION PDXMT (ISET, IPARTON, X, Q, IRT)

C             For ISET = 1, 2 .. , returns sets of Parton Distributions
C             (in the proton) with parton label Iparton (6, 5, ...,0, ..
C             for (t, b, c, s, d, u, g, u-bar, ... t-bar), and kinematic
C             variables (X, Q).   IRT is a return error code.
C
C     Iset =  1, 2, 3, 4 corresponds to the S1, B1, B2, and E1 fits of M
C             Tung (Fermilab-Pub-90/74, IIT-90/11) to NLO in the DIS sch
C
C             5 corresponds to a new set with input strange quark distri
C             being 1/2 of the non-strange sea quarks (as prefered by so
C
C             6, 7, 8, 9 corresponds to the S1, B1, B2, and E1 fits of M
C             Tung (Fermilab-Pub-90/74, IIT-90/11) to NLO in the MSbar s
C
C             10 corresponds to a new set with input strange quark distr
C             being 1/2 of the non-strange sea quarks (as prefered by so
C             in the MSbar scheme
C
C             All the above sets assume a SU(3)-symmetric sea.
C
C             11 corresponds to a set of LO distributions suitable to be
C             with LO hard scattering matrix elements.
C
C  The "lambda" parameter (4-flavors) for each parton distribution set c
C  obtained by making the following FUNCTION call:
C             Alam = Vlambd (Iset, Iorder)
C  where Iset is the (input) set #, Iorder is the (output) order of the
C  for set 11, 2 for all the others), and Alam is the value of the effec
C  lambda for 4 flavors.

C             Details about the 1 - 5 distributions are
C             given in the above-mentioned preprint.
C
#include "pdf/impdp.inc"
      REAL PDZXMT,X,Q

      DIMENSION  THRSLD(0:6)

      DATA (THRSLD(I), I=0,6) / 4*0.0, 1.5, 5.0, 90.0 /

      IFL = IPARTON
      JFL = ABS(IFL)
C                                                   Return 0 if below th
      IF (Q .LE. THRSLD(JFL)) THEN
        PDXMT = 0.0
        RETURN
      ENDIF
C
      IF (IFL .LE. 0) THEN
        VL = 0
      ELSEIF (IFL .LE. 2) THEN
        VL = PDZXMT(ISET, IFL, X, Q, IRT)
      ELSE
        VL = 0
      ENDIF
C
      SEA = PDZXMT (ISET, -JFL, X, Q, IRT)

      PDXMT = VL + SEA

      RETURN
C                         *************************
      END
