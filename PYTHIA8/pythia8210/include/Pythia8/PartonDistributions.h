// PartonDistributions.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for parton densities.
// PDF: base class.
// LHAPDF: derived class for interface to the LHAPDF library.
// GRV94L: derived class for the GRV 94L parton densities.
// CTEQ5L: derived class for the CTEQ 5L parton densities.
// MSTWpdf: derived class for MRST LO*, LO**, MSTW 2008 LO, NLO.
// CTEQ6pdf: derived class for CTEQ 6L, 6L1, 66, CT09 MC1, MC2, (MCS?).
// ProtonPoint: unresolved proton with equivalent photon spectrum.
// GRVpiL: derived class for the GRV LO pion parton densities.
// PomFix: derived class for Q2-independent Pomeron parton densities.
// PomH1FitAB: derived class for the H1 2006 Fit A and Fit B Pomeron PDFs.
// PomH1Jets: derived class for the H1 2007 Jets Pomeron PDFs.
// Lepton: derived class for parton densities inside a lepton.
// LeptonPoint: derived class for unresolved lepton (mainly dummy).
// NNPDF: derived class for the NNPDF2.3 QCD+QED PDF sets.

#ifndef Pythia8_PartonDistributions_H
#define Pythia8_PartonDistributions_H

#include "Pythia8/Basics.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Base class for parton distribution functions.

class PDF {

public:

  // Constructor.
  PDF(int idBeamIn = 2212) {idBeam = idBeamIn; idBeamAbs = abs(idBeam);
    setValenceContent(); idSav = 9; xSav = -1.; Q2Sav = -1.;
    xu = 0.; xd = 0.; xs = 0.; xubar = 0.; xdbar = 0.; xsbar = 0.; xc = 0.;
    xb = 0.; xg = 0.; xlepton = 0.; xgamma = 0.; xuVal = 0.; xuSea = 0.;
    xdVal = 0.; xdSea = 0.; isSet = true; isInit = false;}

  // Destructor.
  virtual ~PDF() {}

  // Confirm that PDF has been set up (important for LHAPDF and H1 Pomeron).
  virtual bool isSetup() {return isSet;}

  // Dynamic choice of meson valence flavours for pi0, K0S, K0L, Pomeron.
  virtual void newValenceContent(int idVal1In, int idVal2In) {
    idVal1 = idVal1In; idVal2 = idVal2In;}

  // Allow extrapolation beyond boundaries. This is optional.
  virtual void setExtrapolate(bool) {}

  // Read out parton density
  virtual double xf(int id, double x, double Q2);

  // Read out valence and sea part of parton densities.
  virtual double xfVal(int id, double x, double Q2);
  virtual double xfSea(int id, double x, double Q2);

  // Check whether x and Q2 values fall inside the fit bounds (LHAPDF6 only).
  virtual bool insideBounds(double, double) {return true;}

  // Access the running alpha_s of a PDF set (LHAPDF6 only).
  virtual double alphaS(double) { return 1.;}

  // Return quark masses used in the PDF fit (LHAPDF6 only).
  virtual double mQuarkPDF(int) { return -1.;}

protected:

  // Allow the LHAPDF class to access these methods.
  friend class LHAPDF;

  // Store relevant quantities.
  int    idBeam, idBeamAbs, idSav, idVal1, idVal2;
  double xSav, Q2Sav;
  double xu, xd, xs, xubar, xdbar, xsbar, xc, xb, xg, xlepton, xgamma,
         xuVal, xuSea, xdVal, xdSea;
  bool   isSet, isInit;

  // Resolve valence content for assumed meson. Possibly modified later.
  void setValenceContent();

  // Update parton densities.
  virtual void xfUpdate(int id, double x, double Q2) = 0;

};

//==========================================================================

// Gives the GRV 94L (leading order) parton distribution function set
// in parametrized form. Authors: M. Glueck, E. Reya and A. Vogt.

class GRV94L : public PDF {

public:

  // Constructor.
  GRV94L(int idBeamIn = 2212) : PDF(idBeamIn) {}

private:

  // Update PDF values.
  void xfUpdate(int , double x, double Q2);

  // Auxiliary routines used during the updating.
  double grvv (double x, double n, double ak, double bk, double a,
    double b, double c, double d);
  double grvw (double x, double s, double al, double be, double ak,
    double bk, double a, double b, double c, double d, double e, double es);
  double grvs (double x, double s, double sth, double al, double be,
    double ak, double ag, double b, double d, double e, double es);

};

//==========================================================================

// Gives the CTEQ 5L (leading order) parton distribution function set
// in parametrized form. Parametrization by J. Pumplin. Authors: CTEQ.

class CTEQ5L : public PDF {

public:

  // Constructor.
  CTEQ5L(int idBeamIn = 2212) : PDF(idBeamIn) {}

private:

  // Update PDF values.
  void xfUpdate(int , double x, double Q2);

};

//==========================================================================

// The MSTWpdf class.
// MRST LO*(*) and MSTW 2008 PDF's, specifically the LO one.
// Original C++ version by Jeppe Andersen.
// Modified by Graeme Watt <watt(at)hep.ucl.ac.uk>.
// Sets available:
// iFit = 1 : MRST LO*  (2007).
// iFit = 2 : MRST LO** (2008).
// iFit = 3 : MSTW 2008 LO, central member.
// iFit = 4 : MSTW 2008 NLO, central member. (Warning!)

class MSTWpdf : public PDF {

public:

  // Constructor.
  MSTWpdf(int idBeamIn = 2212, int iFitIn = 1, 
    string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0) 
    : PDF(idBeamIn) {init( iFitIn,  xmlPath, infoPtr);}

private:

  // Constants: could only be changed in the code itself.
  static const int    np, nx, nq, nqc0, nqb0;
  static const double xmin, xmax, qsqmin, qsqmax, xxInit[65], qqInit[49];

  // Data read in from grid file or set at initialization.
  int    iFit, alphaSorder, alphaSnfmax;
  double mCharm, mBottom, alphaSQ0, alphaSMZ, distance, tolerance,
         xx[65], qq[49], c[13][64][48][5][5];

  // Initialization of data array.
  void init( int iFitIn, string xmlPath, Info* infoPtr);

  // Update PDF values.
  void xfUpdate(int , double x, double Q2);

  // Evaluate PDF of one flavour species.
  double parton(int flavour,double x,double q);
  double parton_interpolate(int flavour,double xxx,double qqq);
  double parton_extrapolate(int flavour,double xxx,double qqq);

  // Auxiliary routines for evaluation.
  int locate(double xx[],int n,double x);
  double polderivative1(double x1, double x2, double x3, double y1,
    double y2, double y3);
  double polderivative2(double x1, double x2, double x3, double y1,
    double y2, double y3);
  double polderivative3(double x1, double x2, double x3, double y1,
    double y2, double y3);

};

//==========================================================================

// The CTEQ6pdf class.
// Sets available:
// iFit = 1 : CTEQ6L
// iFit = 2 : CTEQ6L1
// iFit = 3 : CTEQ66.00 (NLO, central member)
// iFit = 4 : CT09MC1
// iFit = 5 : CT09MC2
// iFit = 6 : CT09MCS (not yet implemented)

class CTEQ6pdf : public PDF {

public:

  // Constructor.
  CTEQ6pdf(int idBeamIn = 2212, int iFitIn = 1, 
    string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0) 
    : PDF(idBeamIn) {init( iFitIn, xmlPath, infoPtr);}

private:

  // Constants: could only be changed in the code itself.
  static const double EPSILON, XPOWER;

  // Data read in from grid file or set at initialization.
  int    iFit, order, nQuark, nfMx, mxVal, nX, nT, nG,
         iGridX, iGridQ, iGridLX, iGridLQ;
  double lambda, mQ[7], qIni, qMax, tv[26], xMin, xv[202], upd[57773],
         xvpow[202], xMinEps, xMaxEps, qMinEps, qMaxEps, fVec[5],
         tConst[9], xConst[9], xLast, qLast;

  // Initialization of data array.
  void init( int iFitIn, string xmlPath, Info* infoPtr);

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

  // Evaluate PDF of one flavour species.
  double parton6(int iParton, double x, double q);

  // Interpolation in grid.
  double polint4F(double xgrid[], double fgrid[], double xin);

};

//==========================================================================

// SA Unresolved proton: equivalent photon spectrum from
// V.M. Budnev, I.F. Ginzburg, G.V. Meledin and V.G. Serbo,
// Phys. Rept. 15 (1974/1975) 181.

class ProtonPoint : public PDF {

public:

  // Constructor.
  ProtonPoint(int idBeamIn = 2212, Info* infoPtrIn = 0) :
              PDF(idBeamIn), m_infoPtr(infoPtrIn) {}

private:

  // Stored value for PDF choice.
  static const double ALPHAEM, Q2MAX, Q20, A, B, C;

  // Update PDF values.
  void xfUpdate(int , double x, double Q2);

  // phi function from Q2 integration.
  double phiFunc(double x, double Q);

  // Info and errors
  Info* m_infoPtr;

};

//==========================================================================

// Gives the GRV 1992 pi+ (leading order) parton distribution function set
// in parametrized form. Authors: Glueck, Reya and Vogt.

class GRVpiL : public PDF {

public:

  // Constructor.
  GRVpiL(int idBeamIn = 221) : PDF(idBeamIn) {}

private:

  // Update PDF values.
  void xfUpdate(int , double x, double Q2);

};

//==========================================================================

// Gives generic Q2-independent Pomeron PDF.

class PomFix : public PDF {

public:

  // Constructor.
  PomFix(int idBeamIn = 990, double PomGluonAIn = 0.,
    double PomGluonBIn = 0., double PomQuarkAIn = 0.,
    double PomQuarkBIn = 0., double PomQuarkFracIn = 0.,
    double PomStrangeSuppIn = 0.) : PDF(idBeamIn),
    PomGluonA(PomGluonAIn), PomGluonB(PomGluonBIn),
    PomQuarkA(PomQuarkAIn), PomQuarkB(PomQuarkBIn),
    PomQuarkFrac(PomQuarkFracIn), PomStrangeSupp(PomStrangeSuppIn)
    {init();}

private:

  // Stored value for PDF choice.
  double PomGluonA, PomGluonB, PomQuarkA, PomQuarkB, PomQuarkFrac,
         PomStrangeSupp, normGluon, normQuark;

  // Initialization of some constants.
  void init();

  // Update PDF values.
  void xfUpdate(int , double x, double);

};

//==========================================================================

// The H1 2006 Fit A and Fit B Pomeron parametrization.
// H1 Collaboration, A. Aktas et al., "Measurement and QCD Analysis of
// the Diffractive Deep-Inelastic Scattering Cross Section at HERA",
// DESY-06-049, Eur. Phys. J. C48 (2006) 715. e-Print: hep-ex/0606004.

class PomH1FitAB : public PDF {

public:

  // Constructor.
 PomH1FitAB(int idBeamIn = 990, int iFit = 1, double rescaleIn = 1.,
   string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0) 
   : PDF(idBeamIn) {rescale = rescaleIn; init( iFit, xmlPath, infoPtr);}

private:

  // Limits for grid in x, in Q2, and data in (x, Q2).
  int    nx, nQ2;
  double rescale, xlow, xupp, dx, Q2low, Q2upp, dQ2;
  double gluonGrid[100][30];
  double quarkGrid[100][30];

  // Initialization of data array.
  void init( int iFit, string xmlPath, Info* infoPtr);

  // Update PDF values.
  void xfUpdate(int , double x, double );

};

//==========================================================================

// The H1 2007 Jets Pomeron parametrization..
// H1 Collaboration, A. Aktas et al., "Dijet Cross Sections and Parton
// Densities in Diffractive DIS at HERA", DESY-07-115, Aug 2007. 33pp.
// Published in JHEP 0710:042,2007. e-Print: arXiv:0708.3217 [hep-ex]

class PomH1Jets : public PDF {

public:

  // Constructor.
  PomH1Jets(int idBeamIn = 990,  double rescaleIn = 1.,
   string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0) 
   : PDF(idBeamIn) {rescale = rescaleIn; init( xmlPath, infoPtr);}

private:

  // Arrays for grid in x, in Q2, and data in (x, Q2).
  double rescale;
  double xGrid[100];
  double Q2Grid[88];
  double gluonGrid[100][88];
  double singletGrid[100][88];
  double charmGrid[100][88];

  // Initialization of data array.
  void init( string xmlPath, Info* infoPtr);

  // Update PDF values.
  void xfUpdate(int id, double x, double );

};

//==========================================================================

// Gives electron (or muon, or tau) parton distribution.

class Lepton : public PDF {

public:

  // Constructor.
  Lepton(int idBeamIn = 11) : PDF(idBeamIn) {}

private:

  // Constants: could only be changed in the code itself.
  static const double ALPHAEM, ME, MMU, MTAU;

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

  // The squared lepton mass, set at initialization.
  double m2Lep;

};

//==========================================================================

// Gives electron (or other lepton) parton distribution when unresolved.

class LeptonPoint : public PDF {

public:

  // Constructor.
  LeptonPoint(int idBeamIn = 11) : PDF(idBeamIn) {}

private:

  // Update PDF values in trivial way.
  void xfUpdate(int , double , double ) {xlepton = 1; xgamma = 0.;}

};

//==========================================================================

// Gives neutrino parton distribution when unresolved (only choice for now).
// Note factor of 2 since only lefthanded implies no spin averaging.

class NeutrinoPoint : public PDF {

public:

  // Constructor.
  NeutrinoPoint(int idBeamIn = 12) : PDF(idBeamIn) {}

private:

  // Update PDF values, with spin factor of 2.
  void xfUpdate(int , double , double ) {xlepton = 2; xgamma = 0.;}

};

//==========================================================================

// The NNPDF class.
// Sets available:
// Leading order QCD+QED Proton PDF sets
// iFit = 1 : NNPDF2.3 QCD+QED LO, alphas(MZ) = 0.130
// iFit = 2 : NNPDF2.3 QCD+QED LO, alphas(MZ) = 0.119
// (Next-to-)Next-to-Leading order QCD+QED Proton PDF sets
// iFit = 3 : NNPDF2.3 QCD+QED NLO, alphas(MZ) = 0.119
// iFit = 4 : NNPDF2.3 QCD+QED NNLO, alphas(MZ) = 0.119
// Code provided by Juan Rojo and Stefano Carrazza.

class NNPDF : public PDF {

public:

  // Constructor.
  NNPDF(int idBeamIn = 2212, int iFitIn = 1, 
    string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0) 
    : PDF(idBeamIn), fPDFGrid(NULL), fXGrid(NULL), fLogXGrid(NULL), 
    fQ2Grid(NULL), fLogQ2Grid(NULL), fRes(NULL) {
    init( iFitIn, xmlPath, infoPtr); };

  // Destructor.
  ~NNPDF() {
    if (fPDFGrid) {
      for (int i = 0; i < fNFL; i++) {
        for (int j = 0; j < fNX; j++)
          if (fPDFGrid[i][j]) delete[] fPDFGrid[i][j];
        if (fPDFGrid[i]) delete[] fPDFGrid[i];
      }
      delete[] fPDFGrid;
    }
    if (fXGrid) delete[] fXGrid;
    if (fLogXGrid) delete[] fLogXGrid;
    if (fQ2Grid) delete[] fQ2Grid;
    if (fLogQ2Grid) delete[] fLogQ2Grid;
    if (fRes) delete[] fRes;
  };

private:

  // Constants: could only be changed in the code itself.
  static const double fXMINGRID;

  // Number of flavors (including photon) and interpolation parameters.
  static const int fNFL = 14;
  static const int fM = 4;
  static const int fN = 2;

  // Variables to be set during code initialization.
  int iFit, fNX, fNQ2;
  double ***fPDFGrid;
  double *fXGrid;
  double *fLogXGrid;
  double *fQ2Grid;
  double *fLogQ2Grid;
  double *fRes;

  // Initialization of data array.
  void init( int iFitIn, string xmlPath, Info* infoPtr);

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

  // Interpolation in the grid for a given PDF flavour.
  void xfxevolve(double x, double Q2);

  // 1D and 2D polynomial interpolation.
  void polint(double xa[], double ya[], int n, double x,
    double& y, double& dy);
  void polin2(double x1a[], double x2a[], double ya[][fN],
    double x1, double x2, double& y, double& dy);

};

//==========================================================================

// LHAPDF plugin interface class.

class LHAPDF : public PDF {

public:

  // Constructor and destructor.
  LHAPDF(int idIn, string pSet, Info* infoPtrIn);
  ~LHAPDF();

  // Confirm that PDF has been set up.
  bool isSetup() {if (pdfPtr) return pdfPtr->isSetup(); return false;}

  // Dynamic choice of meson valence flavours for pi0, K0S, K0L, Pomeron.
  void newValenceContent(int idVal1In, int idVal2In) {
    if (pdfPtr) pdfPtr->newValenceContent(idVal1In, idVal2In);}

  // Allow extrapolation beyond boundaries.
  void setExtrapolate(bool extrapolate) {
    if (pdfPtr) pdfPtr->setExtrapolate(extrapolate);}

  // Read out parton density
  double xf(int id, double x, double Q2) {
    if (pdfPtr) return pdfPtr->xf(id, x, Q2); else return 0;}

  // Read out valence and sea part of parton densities.
  double xfVal(int id, double x, double Q2) {
    if (pdfPtr) return pdfPtr->xfVal(id, x, Q2); else return 0;}
  double xfSea(int id, double x, double Q2) {
    if (pdfPtr) return pdfPtr->xfSea(id, x, Q2); else return 0;}

  // Check whether x and Q2 values fall inside the fit bounds (LHAPDF6 only).
  bool insideBounds(double x, double Q2) {
    if(pdfPtr) return pdfPtr->insideBounds(x, Q2); else return true;}

  // Access the running alpha_s of a PDF set (LHAPDF6 only).
  double alphaS(double Q2) {
    if(pdfPtr) return pdfPtr->alphaS(Q2); else return 1.;}

  // Return quark masses used in the PDF fit (LHAPDF6 only).
  double mQuarkPDF(int idIn) {
    if(pdfPtr) return pdfPtr->mQuarkPDF(idIn); else return -1.;}

private:

  // Resolve valence content for assumed meson.
  void setValenceContent() {if (pdfPtr) pdfPtr->setValenceContent();}

  // Update parton densities.
  void xfUpdate(int id, double x, double Q2) {
    if (pdfPtr) pdfPtr->xfUpdate(id, x, Q2);}

  // Typedefs of the hooks used to access the plugin.
  typedef PDF* NewLHAPDF(int, string, int, Info*);
  typedef void DeleteLHAPDF(PDF*);
  typedef void (*Symbol)();

  // Acccess a plugin library symbol.
  Symbol symbol(string symName);

  // The loaded LHAPDF object, info pointer, and plugin library name.
  PDF   *pdfPtr;
  Info  *infoPtr;
  string libName;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PartonDistributions_H
