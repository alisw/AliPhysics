*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:51  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:15  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C     Version 3 CTEQ distribution function in a parametrized form.

C   By: H.L. Lai, J. Botts, J. Huston, J.G. Morfin, J.F. Owens, J. Qiu,
C       W.K. Tung & H. Weerts;  Preprint MSU-HEP/41024, CTEQ 404

C   This file contains three versions of the same CTEQ3 parton distributions:
C
C Two "front-end" subprograms:
C     FUNCTION Ctq3Pf (Iset, Iparton, X, Q, Irt)
C         returns the PROBABILITY density for a GIVEN flavor;
C     SUBROUTINE Ctq3Pds(Iset, Pdf, XX, QQ, Irt)
C         returns an array of MOMENTUM densities for ALL flavors;
C One lower-level subprogram:
C     FUNCTION Ctq3Pd (Iset, Iprtn, XX, QQ, Irt)
C         returns the MOMENTUM density of a GIVEN valence or sea distribution.

C      One supplementary function to return the QCD lambda parameter
C      concerning these distributions is also included (see below).

C     Although DOUBLE PRECISION is used, conversion to SINGLE PRECISION
C     is straightforward by removing the
C     Implicit Double Precision statements.

C     Since this is an initial distribution of version 3, it is
C     useful for the authors to maintain a record of the distribution
C     list in case there are revisions or corrections.
C     In the interest of maintaining the integrity of this package,
C     please do not freely distribute this program package; instead, refer
C     any interested colleagues to direct their request for a copy to:
C     Lai@cteq11.pa.msu.edu or Tung@msupa.pa.msu.edu.

C   If you have detailed questions concerning these CTEQ3 distributions,
C   or if you find problems/bugs using this initial distribution, direct
C   inquires to Hung-Liang Lai or Wu-Ki Tung.

C     -------------------------------------------
C     Detailed instructions follow.

C     Name convention for CTEQ distributions:  CTEQnSx  where
C        n : version number                      (currently n = 3)
C        S : factorization scheme label: = [M L D] for [MS-bar LO DIS]
c               resp.
C        x : special characteristics, if any
C        (e.g. S(F) for singular (flat) small-x, L for "LEP lambda value")
C        (not applicable to CTEQ3 since only three standard sets are given.)

C    Explanation of functional arguments:

C    Iset is the set label; in this version, Iset = 1, 2, 3
C                           correspond to the following CTEQ global fits:

C          cteq3M  : best fit in the MS-bar scheme
C          cteq3L  : best fit in Leading order QCD
C          cteq3D  : best fit in the DIS scheme

C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)
C  *** WARNING: We use the parton label 2 as D-quark, and 1 as U-quark which
C               might be different with your labels.

C   X, Q are the usual x, Q;
C   Irt is a return error code (see individual modules for explanation).
C
C     ---------------------------------------------

C  Since the QCD Lambda value for the various sets are needed more often than
C  the other parameters in most applications, a special function
C     Wlamd3 (Iset, Iorder, Neff)                    is provided
C  which returns the lambda value for Neff = 4,5,6 effective flavors as well as
C  the order these values pertain to.

C     ----------------------------------------------
C     The range of (x, Q) used in this round of global analysis is, approxi-
C     mately,  0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2 for fixed target
C     experiments and 0.0001 < x < 0.1 from HERA data.

C    The range of (x, Q) used in the reparametrization of the QCD evolved
C    parton distributions is 10E-6 < x < 1 ; 1.6 GeV < Q < 10 TeV.  The
C    functional form of this parametrization is:

C      A0 * x^A1 * (1-x)^A2 * (1 + A3 * x^A4) * [log(1+1/x)]^A5

C   with the A'coefficients being smooth functions of Q.  For heavy quarks,
C   a threshold factor is applied to A0 which simulates the proper Q-dependence
C   of the QCD evolution in that region according to the renormalization
C   scheme defined in Collins-Tung, Nucl. Phys. B278, 934 (1986).

C   Since this function is positive definite and smooth, it provides sensible
C   extrapolations of the parton distributions if they are called beyond
C   the original range in an application. There is no artificial boundaries
C   or sharp cutoff's.
C    ------------------------------------------------

      FUNCTION Ctq3Pf (Iset, Iparton, X, Q, Irt)

C   This function returns the CTEQ parton distributions f^Iset_Iprtn/proton
C   --- the PROBABILITY density

C   (Iset, Iparton, X, Q): explained above;

C    Irt : return error code: see module Ctq3Pd for explanation.

C     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C+SEQ, IMPDP.

      Ifl = Iparton
      JFL = ABS(Ifl)
C                                                             Valence
      IF (Ifl.Eq.1 .or. Ifl.Eq.2) THEN
        VL = Ctq3Pd(Iset, Ifl, X, Q, Irt)
      ELSE
        VL = 0.
      ENDIF
C                                                             Sea
      SEA = Ctq3Pd (Iset, -JFL, X, Q, Irt)
C                                              Full (probability) Distribution
      Ctq3Pf = (VL + SEA) / X

      Return
C                         *************************
      END
