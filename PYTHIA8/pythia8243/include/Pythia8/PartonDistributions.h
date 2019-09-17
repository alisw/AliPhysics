// PartonDistributions.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
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
// CJKL: derived class for CJKL parton densities for photons.

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
  PDF(int idBeamIn = 2212) : idBeam(idBeamIn), idBeamAbs(abs(idBeam)),
    idVal1(), idVal2(), xsVal(), xcVal(), xbVal(), xsSea(), xcSea(), xbSea()
    { setValenceContent();
    idSav = 9; xSav = -1.; Q2Sav = -1.;
    xu = 0.; xd = 0.; xs = 0.; xubar = 0.; xdbar = 0.; xsbar = 0.; xc = 0.;
    xb = 0.; xg = 0.; xlepton = 0.; xgamma = 0.; xuVal = 0.; xuSea = 0.;
    xdVal = 0.; xdSea = 0.; isSet = true; isInit = false;
    hasGammaInLepton = false; }

  // Destructor.
  virtual ~PDF() {}

  // Confirm that PDF has been set up (important for LHAPDF and H1 Pomeron).
  virtual bool isSetup() {return isSet;}

  // Dynamic choice of meson valence flavours for pi0, K0S, K0L, Pomeron.
  virtual void newValenceContent(int idVal1In, int idVal2In) {
    idVal1 = idVal1In; idVal2 = idVal2In;}

  // Allow extrapolation beyond boundaries. This is optional.
  virtual void setExtrapolate(bool) {}

  // Read out parton density.
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

  // Return number of members of this PDF family (LHAPDF6 only).
  virtual int nMembers() { return 1;}

  // Error envelope from PDF uncertainty.
  struct PDFEnvelope {
    double centralPDF, errplusPDF, errminusPDF, errsymmPDF, scalePDF;
    vector<double> pdfMemberVars;
    PDFEnvelope() : centralPDF(-1.0), errplusPDF(0.0), errminusPDF(0.0),
      errsymmPDF(0.0), scalePDF(-1.0), pdfMemberVars(0.0) {};
  };

  // Calculate PDF envelope.
  virtual void calcPDFEnvelope(int, double, double, int) {}
  virtual void calcPDFEnvelope(pair<int,int>, pair<double,double>, double,
    int) {}
  virtual PDFEnvelope getPDFEnvelope() { return PDFEnvelope(); }

  // Approximate photon PDFs by decoupling the scale and x-dependence.
  virtual double gammaPDFxDependence(int, double) { return 0.; }

  // Provide the reference scale for logarithmic Q^2 evolution for photons.
  virtual double gammaPDFRefScale(int) { return 0.; }

  // Sample the valence content for photons.
  virtual int sampleGammaValFlavor(double) { return 0.; }

  // The total x-integrated PDFs. Relevant for MPIs with photon beams.
  virtual double xfIntegratedTotal(double) { return 0.; }

  // Return the sampled value for x_gamma.
  virtual double xGamma() { return 1.; }

  // Keep track of pomeron momentum fraction.
  virtual void xPom(double = -1.0) {}

  // Return accurate and approximated photon fluxes and PDFs.
  virtual double xfFlux(int , double , double )   { return 0.; }
  virtual double xfApprox(int , double , double ) { return 0.; }
  virtual double xfGamma(int , double , double )  { return 0.; }
  virtual double intFluxApprox()                  { return 0.; }

  // Return the kinematical limits and sample Q2 and x.
  virtual double getXmin()              { return 0.; }
  virtual double getXhadr()             { return 0.; }
  virtual double sampleXgamma(double )  { return 0.; }
  virtual double sampleQ2gamma(double ) { return 0.; }

  // Normal PDFs unless gamma inside lepton -> an overestimate for sampling.
  virtual double xfMax(int id, double x, double Q2) { return xf( id, x, Q2); }

  // Normal PDFs unless gamma inside lepton -> Do not sample x_gamma.
  virtual double xfSame(int id, double x, double Q2) { return xf( id, x, Q2); }

  // Allow for new scaling factor for VMD PDFs.
  virtual void setVMDscale(double = 1.) {}

protected:

  // Allow the LHAPDF class to access these methods.
  friend class LHAPDF;

  // Store relevant quantities.
  int    idBeam, idBeamAbs, idSav, idVal1, idVal2;
  double xSav, Q2Sav;
  double xu, xd, xs, xubar, xdbar, xsbar, xc, xb, xg, xlepton, xgamma,
         xuVal, xuSea, xdVal, xdSea;
  bool   isSet, isInit;

  // More valence and sea flavors for photon PDFs.
  double xsVal, xcVal, xbVal, xsSea, xcSea, xbSea;

  // True if a photon beam inside a lepton beam, otherwise set false.
  bool hasGammaInLepton;

  // Resolve valence content for assumed meson. Possibly modified later.
  void setValenceContent();

  // Update parton densities.
  virtual void xfUpdate(int id, double x, double Q2) = 0;

  // Small routine for error printout, depending on infoPtr existing or not.
  void printErr(string errMsg, Info* infoPtr = 0) {
    if (infoPtr !=0) infoPtr->errorMsg(errMsg);
    else cout << errMsg << endl;
  }

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
    : PDF(idBeamIn), iFit(), alphaSorder(), alphaSnfmax(), mCharm(), mBottom(),
    alphaSQ0(), alphaSMZ(), distance(), tolerance(), xx(), qq(),
    c() {init( iFitIn,  xmlPath, infoPtr);}

  // Constructor with a stream.
  MSTWpdf(int idBeamIn, istream& is, Info* infoPtr = 0)
    : PDF(idBeamIn), iFit(), alphaSorder(), alphaSnfmax(), mCharm(), mBottom(),
    alphaSQ0(), alphaSMZ(), distance(), tolerance(), xx(), qq(),
    c() {init( is, infoPtr);}

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

  // Initialization through a stream.
  void init( istream& is, Info* infoPtr);

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
// Proton sets available:
// iFit = 1 : CTEQ6L
// iFit = 2 : CTEQ6L1
// iFit = 3 : CTEQ66.00 (NLO, central member)
// iFit = 4 : CT09MC1
// iFit = 5 : CT09MC2
// iFit = 6 : CT09MCS
// Pomeron sets available (uses same .pds file format as CTEQ6pdf) :
// iFit = 11: ACTWB14
// iFit = 12: ACTWD14
// iFit = 13: ACTWSG14
// iFit = 14: ACTWD19

class CTEQ6pdf : public PDF {

public:

  // Constructor.
  CTEQ6pdf(int idBeamIn = 2212, int iFitIn = 1, double rescaleIn = 1.,
    string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0)
    : PDF(idBeamIn), doExtraPol(false), iFit(), order(), nQuark(), nfMx(),
    mxVal(), nX(), nT(), nG(), iGridX(), iGridQ(), iGridLX(), iGridLQ(),
    rescale(rescaleIn), lambda(), mQ(), qIni(), qMax(), tv(), xMin(), xv(),
    upd(), xvpow(), xMinEps(), xMaxEps(), qMinEps(), qMaxEps(), fVec(),
    tConst(), xConst(), dlx(), xLast(),
    qLast() {init( iFitIn, xmlPath, infoPtr);}

  // Constructor with a stream.
  CTEQ6pdf(int idBeamIn, istream& is, bool isPdsGrid = false,
    Info* infoPtr = 0) : PDF(idBeamIn), doExtraPol(false), iFit(), order(),
    nQuark(), nfMx(), mxVal(), nX(), nT(), nG(), iGridX(), iGridQ(),
    iGridLX(), iGridLQ(), rescale(), lambda(), mQ(), qIni(), qMax(), tv(),
    xMin(), xv(), upd(), xvpow(), xMinEps(), xMaxEps(), qMinEps(), qMaxEps(),
    fVec(), tConst(), xConst(), dlx(), xLast(),
    qLast() { init( is, isPdsGrid, infoPtr); }

  // Allow extrapolation beyond boundaries. This is optional.
  void setExtrapolate(bool doExtraPolIn) {doExtraPol = doExtraPolIn;}

private:

  // Constants: could only be changed in the code itself.
  static const double EPSILON, XPOWER;

  // Data read in from grid file or set at initialization.
  bool   doExtraPol;
  int    iFit, order, nQuark, nfMx, mxVal, nX, nT, nG,
         iGridX, iGridQ, iGridLX, iGridLQ;
  double rescale, lambda, mQ[7], qIni, qMax, tv[26], xMin, xv[202], upd[57773],
         xvpow[202], xMinEps, xMaxEps, qMinEps, qMaxEps, fVec[5],
         tConst[9], xConst[9], dlx, xLast, qLast;

  // Initialization of data array.
  void init( int iFitIn, string xmlPath, Info* infoPtr);

  // Initialization through a stream.
  void init( istream& is, bool isPdsGrid, Info* infoPtr);

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
  ProtonPoint(int idBeamIn = 2212, Info* infoPtrIn = 0)
    : PDF(idBeamIn), infoPtr(infoPtrIn) {}

private:

  // Stored value for PDF choice.
  static const double ALPHAEM, Q2MAX, Q20, A, B, C;

  // Update PDF values.
  void xfUpdate(int , double x, double Q2);

  // phi function from Q2 integration.
  double phiFunc(double x, double Q);

  // Info and errors
  Info* infoPtr;

};

//==========================================================================

// Gives the GRV 1992 pi+ (leading order) parton distribution function set
// in parametrized form. Authors: Glueck, Reya and Vogt.

class GRVpiL : public PDF {

public:

  // Constructor.
  GRVpiL(int idBeamIn = 211, double rescaleIn = 1.) :
    PDF(idBeamIn) {rescale = rescaleIn;}

  // Allow for new rescaling factor of the PDF for VMD beams.
  void setVMDscale(double rescaleIn = 1.) {rescale = rescaleIn;}

private:

  // Rescaling of pion PDF for VMDs.
  double rescale;

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
    PomQuarkFrac(PomQuarkFracIn), PomStrangeSupp(PomStrangeSuppIn),
    normGluon(), normQuark() { init(); }

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
   : PDF(idBeamIn), doExtraPol(false), nx(), nQ2(), rescale(rescaleIn), xlow(),
    xupp(), dx(), Q2low(), Q2upp(), dQ2(), gluonGrid(), quarkGrid()
    { init( iFit, xmlPath, infoPtr); }

  // Constructor with a stream.
 PomH1FitAB(int idBeamIn, double rescaleIn, istream& is,
   Info* infoPtr = 0) : PDF(idBeamIn), doExtraPol(false), nx(), nQ2(),
    rescale(rescaleIn), xlow(),xupp(), dx(), Q2low(), Q2upp(), dQ2(),
    gluonGrid(), quarkGrid() { init( is, infoPtr); }

  // Allow extrapolation beyond boundaries. This is optional.
  void setExtrapolate(bool doExtraPolIn) {doExtraPol = doExtraPolIn;}

private:

  // Limits for grid in x, in Q2, and data in (x, Q2).
  bool   doExtraPol;
  int    nx, nQ2;
  double rescale, xlow, xupp, dx, Q2low, Q2upp, dQ2;
  double gluonGrid[100][30];
  double quarkGrid[100][30];

  // Initialization of data array.
  void init( int iFit, string xmlPath, Info* infoPtr);

  // Initialization through a stream.
  void init( istream& is, Info* infoPtr);

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
  PomH1Jets(int idBeamIn = 990, int iFit = 1, double rescaleIn = 1.,
    string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0)
    : PDF(idBeamIn), doExtraPol(false), rescale(rescaleIn), xGrid(), Q2Grid(),
    gluonGrid(), singletGrid(), charmGrid()
    {init( iFit, xmlPath, infoPtr);}

  // Constructor with a stream.
  PomH1Jets(int idBeamIn, double rescaleIn, istream& is,
    Info* infoPtr = 0) : PDF(idBeamIn), doExtraPol(false), rescale(rescaleIn),
    xGrid(), Q2Grid(), gluonGrid(), singletGrid(), charmGrid()
    { init( is, infoPtr); }

  // Allow extrapolation beyond boundaries. This is optional.
  void setExtrapolate(bool doExtraPolIn) {doExtraPol = doExtraPolIn;}

private:

  // Arrays for grid in x, in Q2, and data in (x, Q2).
  bool   doExtraPol;
  double rescale;
  double xGrid[100];
  double Q2Grid[88];
  double gluonGrid[100][88];
  double singletGrid[100][88];
  double charmGrid[100][88];

  // Initialization of data array.
  void init( int iFit, string xmlPath, Info* infoPtr);

  // Initialization through a stream.
  void init( istream& is, Info* infoPtr);

  // Update PDF values.
  void xfUpdate(int id, double x, double );

};

//==========================================================================

class PomHISASD : public PDF {

public:

  // Basic constructor
  PomHISASD(int idBeamIn, PDF* ppdf, Settings & settings,Info* infoPtrIn = 0)
    : PDF(idBeamIn), pPDFPtr(ppdf), xPomNow(-1.0), hixpow(4.0), newfac(1.0) {
    infoPtr = infoPtrIn;
    hixpow = settings.parm("PDF:PomHixSupp");
    if ( settings.mode("Angantyr:SASDmode") == 3 ) newfac =
      log(settings.parm("Beams:eCM")/settings.parm("Diffraction:mMinPert"));
    if ( settings.mode("Angantyr:SASDmode") == 4 ) newfac = 0.0;
  }

  // Delete also the proton PDF
  ~PomHISASD() { delete pPDFPtr; }

  // (re-)Set the x_pomeron value.
  void xPom(double xpom = -1.0) { xPomNow = xpom; }

private:

  // The proton PDF.
  PDF* pPDFPtr;

  // The momenum fraction if the corresponding pomeron .
  double xPomNow;

  // The high-x suppression power.
  double hixpow;

  // Special options.
  double newfac;

  // Report possible errors.
  Info* infoPtr;

  // Update PDF values.
  void xfUpdate(int , double x, double Q2);

};

//==========================================================================

// Gives electron (or muon, or tau) parton distribution.

class Lepton : public PDF {

public:

  // Constructor.
  Lepton(int idBeamIn = 11) : PDF(idBeamIn), m2Lep(), Q2maxGamma(), infoPtr(),
    rndmPtr() {}

  // Constructor with further info.
  Lepton(int idBeamIn, double Q2maxGammaIn, Info* infoPtrIn,
    Rndm* rndmPtrIn = 0) : PDF(idBeamIn), m2Lep() { Q2maxGamma = Q2maxGammaIn;
    infoPtr = infoPtrIn; rndmPtr = rndmPtrIn; }

  // Sample the Q2 value.
  double sampleQ2gamma(double Q2min)
    { return Q2min * pow(Q2maxGamma / Q2min, rndmPtr->flat()); }

private:

  // Constants: could only be changed in the code itself.
  static const double ALPHAEM, ME, MMU, MTAU;

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

  // The squared lepton mass, set at initialization.
  double m2Lep, Q2maxGamma;

  // Pointer to info, needed to get sqrt(s) to fix x_gamma limits.
  Info* infoPtr;

  // Pointer to random number generator, needed to sample Q2.
  Rndm* rndmPtr;

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
    : PDF(idBeamIn), iFit(), fNX(), fNQ2(), fPDFGrid(NULL), fXGrid(NULL),
    fLogXGrid(NULL), fQ2Grid(NULL), fLogQ2Grid(NULL), fRes(NULL) {
    init( iFitIn, xmlPath, infoPtr); };

  // Constructor with a stream.
  NNPDF(int idBeamIn, istream& is, Info* infoPtr = 0)
    : PDF(idBeamIn), iFit(), fNX(), fNQ2(), fPDFGrid(NULL), fXGrid(NULL),
    fLogXGrid(NULL), fQ2Grid(NULL), fLogQ2Grid(NULL), fRes(NULL) {
    init( is, infoPtr); };

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

  // Initialization through a stream.
  void init( istream& is, Info* infoPtr);

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

  // Return quark masses used in the PDF fit (LHAPDF6 only).
  int nMembers() {
    if(pdfPtr) return pdfPtr->nMembers(); else return 1;}


  // Calculate PDF envelope.
  void calcPDFEnvelope(int idNow, double xNow, double Q2Now, int valSea) {
    if (pdfPtr) pdfPtr->calcPDFEnvelope(idNow, xNow, Q2Now, valSea);}
  void calcPDFEnvelope(pair<int,int> idNows, pair<double,double> xNows,
    double Q2Now, int valSea) {
    if (pdfPtr) pdfPtr->calcPDFEnvelope(idNows,xNows,Q2Now,valSea);}
  PDFEnvelope getPDFEnvelope() { if (pdfPtr) return pdfPtr->getPDFEnvelope();
    else return PDFEnvelope(); }

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

  // The loaded LHAPDF object, info pointer, and plugin library and name.
  PDF   *pdfPtr;
  Info  *infoPtr;
  string libName;
  void  *lib;

};

//==========================================================================

// Gives the CJKL leading order parton distribution function set
// in parametrized form for the real photons. Authors: F.Cornet, P.Jankowski,
// M.Krawczyk and A.Lorca, Phys. Rev. D68: 014010, 2003.

class CJKL : public PDF {

public:

  // Constructor. Needs the randon number generator to sample valence content.
  CJKL(int idBeamIn = 22, Rndm* rndmPtrIn = 0 ) : PDF(idBeamIn) {
    rndmPtr = rndmPtrIn; }

  // Functions to approximate pdfs for ISR.
  double gammaPDFxDependence(int id, double);
  double gammaPDFRefScale(int);

  // Set the valence content for photons.
  int sampleGammaValFlavor(double Q2);

  // The total x-integrated PDFs. Relevant for MPIs with photon beams.
  double xfIntegratedTotal(double Q2);

private:

  // Parameters related to the fit.
  static const double ALPHAEM, Q02, Q2MIN, Q2REF, LAMBDA, MC, MB;

  // Pointer to random number generator used for valence sampling.
  Rndm *rndmPtr;

  // Update PDF values.
  void xfUpdate(int , double x, double Q2);

  // Functions for updating the point-like part.
  double pointlikeG(double x, double s);
  double pointlikeU(double x, double s);
  double pointlikeD(double x, double s);
  double pointlikeC(double x, double s, double Q2);
  double pointlikeB(double x, double s, double Q2);

  // Functions for updating the hadron-like part.
  double hadronlikeG(double x, double s);
  double hadronlikeSea(double x, double s);
  double hadronlikeVal(double x, double s);
  double hadronlikeC(double x, double s, double Q2);
  double hadronlikeB(double x, double s, double Q2);

};

//==========================================================================

// The LHAGrid1 can be used to read files in the LHAPDF6 lhagrid1 format,
// assuming that the same x grid is used for all Q subgrids.
// Results are not identical with LHAPDF6, owing to different interpolation.

class LHAGrid1 : public PDF {

public:

  // Constructor.
  LHAGrid1(int idBeamIn = 2212, string pdfWord = "void",
    string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0)
    : PDF(idBeamIn), doExtraPol(false), nx(), nq(), nqSub(), xMin(), xMax(),
    qMin(), qMax(), pdfVal(), pdfGrid(NULL),
    pdfSlope(NULL) { init( pdfWord, xmlPath, infoPtr); };

  // Constructor with a stream.
  LHAGrid1(int idBeamIn, istream& is, Info* infoPtr = 0)
    : PDF(idBeamIn), doExtraPol(false), nx(), nq(), nqSub(), xMin(), xMax(),
    qMin(), qMax(), pdfVal(), pdfGrid(NULL),
    pdfSlope(NULL) { init( is, infoPtr); };

  // Destructor.
  ~LHAGrid1() { if (pdfGrid) { for (int iid = 0; iid < 12; ++iid) {
    for (int ix = 0; ix < nx; ++ix) delete[] pdfGrid[iid][ix];
    delete[] pdfGrid[iid]; } delete[] pdfGrid; }
    if (pdfSlope) { for (int iid = 0; iid < 12; ++iid) delete[] pdfSlope[iid];
    delete[] pdfSlope;} };

  // Allow extrapolation beyond boundaries. This is optional.
  void setExtrapolate(bool doExtraPolIn) {doExtraPol = doExtraPolIn;}

private:

  // Variables to be set during code initialization.
  bool   doExtraPol;
  int    nx, nq, nqSub;
  vector<int> nqSum;
  double xMin, xMax, qMin, qMax, pdfVal[12];
  vector<double> xGrid, lnxGrid, qGrid, lnqGrid, qDiv;
  double*** pdfGrid;
  double** pdfSlope;

  // Initialization of data array.
  void init( string pdfSet, string xmlPath, Info* infoPtr);

  // Initialization through a stream.
  void init( istream& is, Info* infoPtr);

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

  // Interpolation in the grid for a given PDF flavour.
  void xfxevolve(double x, double Q2);

};

//==========================================================================

// Convolution with photon flux from leptons and photon PDFs.
// Photon flux from equivalent photon approximation (EPA).
// Contains a pointer to a photon PDF set and samples the
// convolution integral event-by-event basis.
// Includes also a overestimate for the PDF set in order to set up
// the phase-space sampling correctly.

class Lepton2gamma : public PDF {

public:

  // Constructor.
  Lepton2gamma(int idBeamIn, double m2leptonIn, double Q2maxGamma,
    PDF* gammaPDFPtrIn, Info* infoPtrIn, Rndm* rndmPtrIn) : PDF(idBeamIn),
    m2lepton(m2leptonIn), Q2max(Q2maxGamma), xGm(), sampleXgamma(true),
    gammaPDFPtr(gammaPDFPtrIn), rndmPtr(rndmPtrIn), infoPtr(infoPtrIn)
    { hasGammaInLepton = true; }

  // Overload the member function definitions where relevant.
  void xfUpdate(int id, double x, double Q2);
  double xGamma() { return xGm; }
  double xfMax(int id, double x, double Q2);
  double xfSame(int id, double x, double Q2);

  // Sample the Q2 value.
  double sampleQ2gamma(double Q2min)
    { return Q2min * pow(Q2max / Q2min, rndmPtr->flat()); }

private:

  // Parameters for convolution.
  static const double ALPHAEM, Q2MIN;
  double m2lepton, Q2max, xGm;

  // Sample new value for x_gamma.
  bool sampleXgamma;

  // Photon PDFs with the photon flux is convoluted with.
  PDF* gammaPDFPtr;

  // Pointer to random number generator used for sampling x_gamma.
  Rndm* rndmPtr;

  // Pointer to info, needed to get sqrt(s) to fix x_gamma limits.
  Info* infoPtr;

};

//==========================================================================

// Gives electron (or other lepton) parton distribution when unresolved.

class GammaPoint : public PDF {

public:

  // Constructor.
  GammaPoint(int idBeamIn = 22) : PDF(idBeamIn) {}

private:

  // Update PDF values in trivial way.
  void xfUpdate(int , double , double ) { xgamma = 1.;}

};

//==========================================================================

// Equivalent photon approximation for sampling with external photon flux.

class EPAexternal : public PDF {

public:

  // Constructor.
  EPAexternal(int idBeamIn, double m2In, PDF* gammaFluxPtrIn,
    PDF* gammaPDFPtrIn, Settings* settingsPtrIn, Info* infoPtrIn,
    Rndm* rndmPtrIn ) : PDF(idBeamIn), m2(m2In), Q2max(), Q2min(), xMax(),
    xMin(), xHadr(), norm(), xPow(), xCut(), norm1(), norm2(),
    integral1(), integral2(), bmhbarc(), gammaFluxPtr(gammaFluxPtrIn),
    gammaPDFPtr(gammaPDFPtrIn), rndmPtr(rndmPtrIn), infoPtr(infoPtrIn),
    settingsPtr(settingsPtrIn) { hasGammaInLepton = true; init(); }

  // Update PDFs.
  void xfUpdate(int , double x, double Q2);

  // External flux and photon PDFs, and approximated flux for sampling.
  double xfFlux(int id, double x, double Q2 = 1.);
  double xfGamma(int id, double x, double Q2);
  double xfApprox(int id, double x, double Q2);
  double intFluxApprox();

  // Kinematics.
  double getXmin()          { return xMin; }
  double getXhadr()         { return xHadr; }

  // Sampling of the x and Q2 according to differential flux.
  double sampleXgamma(double xMinIn);
  double sampleQ2gamma(double )
    { return Q2min * pow(Q2max / Q2min, rndmPtr->flat()); }

private:

  // Initialization.
  void init();

  // Kinematics.
  static const double ALPHAEM;
  double m2, Q2max, Q2min, xMax, xMin, xHadr, norm, xPow, xCut,
         norm1, norm2, integral1, integral2, bmhbarc;
  int    approxMode;

  // Photon Flux and PDF.
  PDF* gammaFluxPtr;
  PDF* gammaPDFPtr;

  // Pointer to random number generator used for sampling x_gamma.
  Rndm* rndmPtr;

  // Pointer to info, needed to get sqrt(s) to fix x_gamma limits.
  Info* infoPtr;

  // Pointer to settings to get Q2max.
  Settings* settingsPtr;

};

//==========================================================================

// A derived class for nuclear PDFs. Needs a pointer for (free) proton PDFs.

class nPDF : public PDF {

public:

  // Constructor.
  nPDF(int idBeamIn = 2212, PDF* protonPDFPtrIn = 0) : PDF(idBeamIn), ruv(),
    rdv(), ru(), rd(), rs(), rc(), rb(), rg(), a(), z(), za(), na(),
    protonPDFPtr() { initNPDF(protonPDFPtrIn); }

  // Update parton densities.
  void xfUpdate(int id, double x, double Q2);

  // Update nuclear modifications.
  virtual void rUpdate(int, double, double) = 0;

  // Initialize the nPDF-related members.
  void initNPDF(PDF* protonPDFPtrIn = 0);

  // Return the number of protons and nucleons.
  int getA() {return a;}
  int getZ() {return z;}

  // Set (and reset) the ratio of protons to nucleons to study nuclear
  // modifications of protons (= 1.0) and neutrons (= 0.0). By default Z/A.
  void setMode(double zaIn) { za = zaIn; na = 1. - za; }
  void resetMode() { za = double(z)/double(a); na = double(a-z)/double(a); }

protected:

  // The nuclear modifications for each flavour, modified by derived nPDF
  // classes.
  double ruv, rdv, ru, rd, rs, rc, rb, rg;

private:

  // The nuclear mass number and number of protons (charge) and normalized
  // number of protons and neutrons.
  int a, z;
  double za, na;

  // Pointer to (free) proton PDF.
  PDF* protonPDFPtr;

};

//==========================================================================

// Isospin modification with nuclear beam, i.e. no other modifications
// but correct number of protons and neutrons.

class Isospin : public nPDF {

public:

  // Constructor.
  Isospin(int idBeamIn = 2212, PDF* protonPDFPtrIn = 0)
    : nPDF(idBeamIn, protonPDFPtrIn) {}

  // Only the Isospin effect so no need to do anything here.
  void rUpdate(int , double , double ) {}
};

//==========================================================================

// Nuclear modifications from EPS09 fit.

class EPS09 : public nPDF {

public:

  // Constructor.
  EPS09(int idBeamIn = 2212, int iOrderIn = 1, int iSetIn = 1,
    string xmlPath = "../share/Pythia8/xmldoc/", PDF* protonPDFPtrIn = 0,
    Info* infoPtrIn = 0) : nPDF(idBeamIn, protonPDFPtrIn), iSet(), iOrder(),
    grid(), infoPtr(infoPtrIn) {init(iOrderIn, iSetIn, xmlPath);}

  // Update parton densities.
  void rUpdate(int id, double x, double Q2);

  // Use other than central set to study uncertainties.
  void setErrorSet(int iSetIn) {iSet = iSetIn;}

private:

  // Parameters related to the fit.
  static const double Q2MIN, Q2MAX, XMIN, XMAX, XCUT;
  static const int Q2STEPS, XSTEPS;

  // Set parameters and the grid.
  int iSet, iOrder;
  double grid[31][51][51][8];

  // Pointer to info for possible error messages.
  Info* infoPtr;

  // Initialize with given inputs.
  void init(int iOrderIn, int iSetIn, string xmlPath);

  // Interpolation algorithm.
  double polInt(double* fi, double* xi, int n, double x);
};

//==========================================================================

// Nuclear modifications from EPPS16 fit.

class EPPS16 : public nPDF {

public:

  // Constructor.
  EPPS16(int idBeamIn = 2212, int iSetIn = 1,
    string xmlPath = "../share/Pythia8/xmldoc/", PDF* protonPDFPtrIn = 0,
    Info* infoPtrIn = 0) : nPDF(idBeamIn, protonPDFPtrIn), iSet(), grid(),
    logQ2min(), loglogQ2maxmin(), logX2min(), infoPtr(infoPtrIn)
    { init(iSetIn, xmlPath); }

  // Update parton densities.
  void rUpdate(int id, double x, double Q2);

  // Use other than central set to study uncertainties.
  void setErrorSet(int iSetIn) {iSet = iSetIn;}

private:

  // Parameters related to the fit.
  static const double Q2MIN, Q2MAX, XMIN, XMAX, XCUT;
  static const int Q2STEPS, XSTEPS, NINTQ2, NINTX, NSETS;

  // Set parameters and the grid.
  int iSet;
  double grid[41][31][80][8];
  double logQ2min, loglogQ2maxmin, logX2min;

  // Pointer to info for possible error messages.
  Info* infoPtr;

  // Initialize with given inputs.
  void init(int iSetIn, string xmlPath);

  // Interpolation algorithm.
  double polInt(double* fi, double* xi, int n, double x);
};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PartonDistributions_H
