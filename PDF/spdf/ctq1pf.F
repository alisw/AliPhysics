*
* $Id$
*
* $Log$
* Revision 1.1.1.5  2000/04/19 13:27:52  mclareni
* Import version 8.04
*
* Revision 1.1.1.1  1996/04/12 15:29:08  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      FUNCTION Ctq1pf (ISET, IPARTON, X, Q, IRT)

C     CTEQ distribution function in a parametrized form.

C     (No data tables are needed.)

C The returned function value is the PROBABILITY density for a given FLAVOR.

C  !! A companion function (next module), which this one depends on,
C  !!        ctq1pd (Iset, Iparton, X, Q, Irt)
C  !! gives the VALENCE and SEA MOMENTUM FRACTION distributions.

C  \\  A parallel (independent) program Ctq1Ps (not included in this file)
C  ||  in Subroutine form is also available.
C  ||  It returns ALL the parton flavors at once in an array form.
C  //  See details in that separate file if you are interested.

C Ref.: "CTEQ Parton Distributions and Flavor Dependence of the Sea Quarks"
C       by: J. Botts, J.G. Morfin, J.F. Owens, J. Qiu, W.K. Tung & H. Weerts
C       MSUHEP-92-27, Fermilab-Pub-92/371, FSU-HEP-92-1225, ISU-NP-92-17

C     Since this is an initial distribution, and there may be updates, it is
C     useful for the authors to maintain a record of the distribution list.
C     Please do not freely distribute this program package; instead, refer any
C     interested colleague to direct their request for a copy to:
C     Botts@msupa.pa.msu.edu  or  Botts@msupa (bitnet)  or  MSUHEP::Botts

C    If you have any questions concerning these distributions, direct inquires
C    to Jim Botts or Wu-Ki Tung (username Tung at same E-mail nodes as above).

C$Header$
C$Log$
CRevision 1.1.1.5  2000/04/19 13:27:52  mclareni
CImport version 8.04
C
CRevision 1.1.1.1  1996/04/12 15:29:08  plothow
CVersion 7.01
C
C Revision 1.2  93/02/26  10:42:43  wkt
C Version with heavy quark threshold factor and faster algorithm.
C
C Revision 1.1  93/02/14  17:30:21  botts
C The new Faster version.
C Revision 1.0  93/02/08  18:35:25  wkt
C Initial revision

C   This function returns the CTEQ parton distributions f^Iset_Iprtn/proton
C     where Iset (= 1, 2, ..., 5) is the set label;

C       Name convention for CTEQ distributions:  CTEQnSx  where
C           n : version number                      (currently n = 1)
C           S : factorization scheme label: = [M D L] for [MS-bar DIS LO]
C               resp.
C           x : special characteristics, if any
C                    (e.g. S for singular gluon, L for "LEP lambda value")

C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)

C X, Q are the usual x, Q; Irt is a return error code (not implemented yet).

C --> Iset = 1, 2, 3, 4, 5 correspond to the following CTEQ global fits:
C     cteq1M, cteq1MS, cteq1ML, cteq1D, cteq1L  respectively.

C --> QCD parameters for parton distribution set Iset can be obtained inside
C         the user's program by:
C     Dum = PrCtq1
C    >        (Iset, Iord, Ischeme, MxFlv,
C    >         Alam4, Alam5, Alam6, Amas4, Amas5, Amas6,
C    >         Xmin, Qini, Qmax, ExpNor)
C     where all but the first argument are output parameters.
C     They should be self-explanary -- see details in next module.

C     The range of (x, Q) used in this round of global analysis is, approxi-
C     mately,  0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2.

C    The range of (x, Q) used in the reparametrization of the QCD evolved
C    parton distributions is 10E-5 < x < 1 ; 2 GeV < Q < 1 TeV.  The
C    functional form of this parametrization is:

C      A0 * x^A1 * (1-x)^A2 * (1 + A3 * x^A4) * [log(1+1/x)]^A5

C   with the A'coefficients being smooth functions of Q.  For heavy quarks,
C   an additional threshold factor is applied which simulate the Q-dependence
C   of the QCD evolution in that region.

C   Since this function is positive definite and smooth, it provides sensible
C   extrapolations of the parton distributions if they are called beyond
C   the original range in an application. There is no artificial boundaries
C   or sharp cutoff's.

C     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C+SEQ, IMPDP.
      REAL Ctq1Pd,X,Q

      IFL = IPARTON
      JFL = ABS(IFL)
C                                                             Valence
      IF (IFL .LE. 0) THEN
        VL = 0.
      ELSEIF (IFL .LE. 2) THEN
        VL = Ctq1Pd (ISET, IFL, X, Q, IRT)
      ELSE
        VL = 0.
      ENDIF
C                                                                         Sea
      SEA = Ctq1Pd (ISET, -JFL, X, Q, IRT)
C                                              Full (probability) Distribution
      Ctq1pf = (VL + SEA) / X

      RETURN
C                         *************************
      END
