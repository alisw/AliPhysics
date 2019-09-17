// BeamParticle.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for information on incoming beams.
// ResolvedParton: an initiator or remnant in beam.
// BeamParticle: contains partons, parton densities, etc.

#ifndef Pythia8_BeamParticle_H
#define Pythia8_BeamParticle_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// This class holds info on a parton resolved inside the incoming beam,
// i.e. either an initiator (part of a hard or a multiparton interaction)
// or a remnant (part of the beam remnant treatment).

// The companion code is -1 from onset and for g, is -2 for an unmatched
// sea quark, is >= 0 for a matched sea quark, with the number giving the
// companion position, and is -3 for a valence quark.

// Rescattering partons properly do not belong here, but bookkeeping is
// simpler with them, so they are stored with companion code -10.

class ResolvedParton {

public:

  // Constructor.
  ResolvedParton( int iPosIn = 0, int idIn = 0, double xIn = 0.,
    int companionIn = -1) : iPosRes(iPosIn), idRes(idIn), xRes(xIn),
    companionRes(companionIn), xqCompRes(0.), mRes(0.), factorRes(1.),
    colRes(0), acolRes(0) { }

  // Set info on initiator or remnant parton.
  void iPos( int iPosIn) {iPosRes = iPosIn;}
  void id( int idIn) {idRes = idIn;}
  void x( double xIn) {xRes = xIn;}
  void update( int iPosIn, int idIn, double xIn) {iPosRes = iPosIn;
    idRes = idIn; xRes = xIn;}
  void companion( int companionIn) {companionRes = companionIn;}
  void xqCompanion( double xqCompIn) {xqCompRes = xqCompIn;}
  void p(Vec4 pIn) {pRes = pIn;}
  void px(double pxIn) {pRes.px(pxIn);}
  void py(double pyIn) {pRes.py(pyIn);}
  void pz(double pzIn) {pRes.pz(pzIn);}
  void e(double eIn) {pRes.e(eIn);}
  void m(double mIn) {mRes = mIn;}
  void col(int colIn) {colRes = colIn;}
  void acol(int acolIn) {acolRes = acolIn;}
  void cols(int colIn = 0,int acolIn = 0)
    {colRes = colIn; acolRes = acolIn;}
  void scalePT( double factorIn) {pRes.px(factorIn * pRes.px());
    pRes.py(factorIn * pRes.py()); factorRes *= factorIn;}
  void scaleX( double factorIn) {xRes *= factorIn;}

  // Get info on initiator or remnant parton.
  int    iPos()        const {return iPosRes;}
  int    id()          const {return idRes;}
  double x()           const {return xRes;}
  int    companion()   const {return companionRes;}
  bool   isValence()   const {return (companionRes == -3);}
  bool   isUnmatched() const {return (companionRes == -2);}
  bool   isCompanion() const {return (companionRes >= 0);}
  bool   isFromBeam()  const {return (companionRes > -10);}
  double xqCompanion() const {return xqCompRes;}
  Vec4   p()           const {return pRes;}
  double px()          const {return pRes.px();}
  double py()          const {return pRes.py();}
  double pz()          const {return pRes.pz();}
  double e()           const {return pRes.e();}
  double m()           const {return mRes;}
  double pT()          const {return pRes.pT();}
  double mT2()         const {return (mRes >= 0.)
    ? mRes*mRes + pRes.pT2() : - mRes*mRes + pRes.pT2();}
  double pPos()        const {return pRes.e() +  pRes.pz();}
  double pNeg()        const {return pRes.e() -  pRes.pz();}
  int    col()         const {return colRes;}
  int    acol()        const {return acolRes;}
  double pTfactor()    const {return factorRes;}
  bool hasCol()        const {return (idRes == 21 || (idRes > 0 && idRes < 9)
    || (-idRes > 1000 && -idRes < 10000 && (-idRes/10)%10 == 0));}
  bool hasAcol()       const {return (idRes == 21 || (-idRes > 0 && -idRes < 9)
    || (idRes > 1000 && idRes < 10000 && (idRes/10)%10 == 0));}

private:

  // Properties of a resolved parton.
  int    iPosRes, idRes;
  double xRes;
  // Companion code and distribution value, if any.
  int    companionRes;
  double xqCompRes;
  // Four-momentum and mass; for remnant kinematics construction.
  Vec4   pRes;
  double mRes, factorRes;
  // Colour codes.
  int   colRes, acolRes;

};

//==========================================================================

// This class holds info on a beam particle in the evolution of
// initial-state radiation and multiparton interactions.

class BeamParticle {

public:

  // Constructor.
  BeamParticle() : infoPtr(), particleDataPtr(), rndmPtr(), pdfBeamPtr(),
    pdfHardBeamPtr(), pdfUnresBeamPtr(), pdfBeamPtrSave(),
    pdfHardBeamPtrSave(), flavSelPtr(), allowJunction(), beamJunction(),
    maxValQuark(), companionPower(), valencePowerMeson(), valencePowerUinP(),
    valencePowerDinP(), valenceDiqEnhance(), pickQuarkNorm(), pickQuarkPower(),
    diffPrimKTwidth(), diffLargeMassSuppress(), beamSat(), gluonPower(),
    xGluonCutoff(), idBeam(), idBeamAbs(), idVMDBeam(), mBeam(), mVMDBeam(),
    scaleVMDBeam(), isUnresolvedBeam(), isLeptonBeam(), isHadronBeam(),
    isMesonBeam(), isBaryonBeam(), isGammaBeam(), nValKinds(), idVal(), nVal(),
    idSave(), iSkipSave(), nValLeft(), xqgTot(), xqVal(), xqgSea(),
    xqCompSum(), doISR(), doMPI(), doND(), isResolvedGamma(),
    hasResGammaInBeam(), isResUnres(), hasVMDstateInBeam(), pTminISR(),
    pTminMPI(), pT2gm2qqbar(), iGamVal(), iPosVal(), gammaMode(), xGm(),
    Q2gm(), kTgamma(), phiGamma(), resolved(), nInit(0), hasJunctionBeam(),
    junCol(), nJuncs(), nAjuncs(), nDiffJuncs(), allowBeamJunctions(),
    Q2ValFracSav(-1.), uValInt(), dValInt(), idVal1(), idVal2(), idVal3(),
    zRel(), pxRel(), pyRel() { }

  // Initialize data on a beam particle and save pointers.
  void init( int idIn, double pzIn, double eIn, double mIn,
    Info* infoPtrIn, Settings& settings, ParticleData* particleDataPtrIn,
    Rndm* rndmPtrIn, PDF* pdfInPtr, PDF* pdfHardInPtr, bool isUnresolvedIn,
    StringFlav* flavSelPtrIn);

  // Initialize only the two pdf pointers.
  void initPDFPtr(PDF* pdfInPtr, PDF* pdfHardInPtr) {
    pdfBeamPtr = pdfInPtr; pdfHardBeamPtr = pdfHardInPtr; }

  // Initialize additional PDF pointer for unresolved beam.
  void initUnres(PDF* pdfUnresInPtr);

  // For mesons like pi0 valence content varies from event to event.
  void newValenceContent();

  // Set new pZ and E, but keep the rest the same.
  void newPzE( double pzIn, double eIn) {pBeam = Vec4( 0., 0., pzIn, eIn);}

  // Set new mass. Used with photons when virtuality is sampled.
  void newM( double mIn) { mBeam = mIn; }

  // Member functions for output.
  int id()            const {return idBeam;}
  int  idVMD()        const {return idVMDBeam;}
  Vec4 p()            const {return pBeam;}
  double px()         const {return pBeam.px();}
  double py()         const {return pBeam.py();}
  double pz()         const {return pBeam.pz();}
  double e()          const {return pBeam.e();}
  double m()          const {return mBeam;}
  double mVMD()       const {return mVMDBeam;}
  double scaleVMD()   const {return scaleVMDBeam;}
  bool isLepton()     const {return isLeptonBeam;}
  bool isUnresolved() const {return isUnresolvedBeam;}
  // As hadrons here we only count those we know how to handle remnants for.
  bool isHadron()     const {return isHadronBeam;}
  bool isMeson()      const {return isMesonBeam;}
  bool isBaryon()     const {return isBaryonBeam;}
  bool isGamma()      const {return isGammaBeam;}
  bool hasResGamma()  const {return hasResGammaInBeam;}
  bool hasVMDstate()  const {return hasVMDstateInBeam;}

  // Maximum x remaining after previous MPI and ISR, plus safety margin.
  double xMax(int iSkip = -1);

  // Special hard-process parton distributions (can agree with standard ones).
  double xfHard(int idIn, double x, double Q2)
    {return pdfHardBeamPtr->xf(idIn, x, Q2);}

  // Overestimate for PDFs. Same as normal except photons inside leptons.
  double xfMax(int idIn, double x, double Q2)
    {return pdfHardBeamPtr->xfMax(idIn, x, Q2);}

  // Accurate and approximated photon flux and PDFs.
  double xfFlux(int idIn, double x, double Q2)
    {return pdfHardBeamPtr->xfFlux(idIn, x, Q2);}
  double xfApprox(int idIn, double x, double Q2)
    {return pdfHardBeamPtr->xfApprox(idIn, x, Q2);}
  double xfGamma(int idIn, double x, double Q2)
    {return pdfHardBeamPtr->xfGamma(idIn, x, Q2);}

  // Do not sample the x_gamma value to get correct cross section with
  // possible second call.
  double xfSame(int idIn, double x, double Q2)
    {return pdfHardBeamPtr->xfSame(idIn, x, Q2);}

  // Standard parton distributions.
  double xf(int idIn, double x, double Q2)
    {return pdfBeamPtr->xf(idIn, x, Q2);}

  // Ditto, split into valence and sea parts (where gluon counts as sea).
  double xfVal(int idIn, double x, double Q2)
    {return pdfBeamPtr->xfVal(idIn, x, Q2);}
  double xfSea(int idIn, double x, double Q2)
    {return pdfBeamPtr->xfSea(idIn, x, Q2);}

  // Rescaled parton distributions, as needed for MPI and ISR.
  // For ISR also allow split valence/sea, and only return relevant part.
  double xfMPI(int idIn, double x, double Q2)
    {return xfModified(-1, idIn, x, Q2);}
  double xfISR(int indexMPI, int idIn, double x, double Q2)
    {return xfModified( indexMPI, idIn, x, Q2);}

  // Check whether x and Q2 values fall inside the fit bounds (LHAPDF6 only).
  bool insideBounds(double x, double Q2)
    {return pdfBeamPtr->insideBounds(x,Q2);}

  // Access the running alpha_s of a PDF set (LHAPDF6 only).
  double alphaS(double Q2) {return pdfBeamPtr->alphaS(Q2);}

  // Return quark masses used in the PDF fit (LHAPDF6 only).
  double mQuarkPDF(int idIn) {return pdfBeamPtr->mQuarkPDF(idIn);}

  // Return number of members in PDF family (LHAPDF6 only).
  int nMembers() {return pdfBeamPtr->nMembers();}

  // Calculate envelope of PDF predictions
  void calcPDFEnvelope(int idNow, double xNow, double Q2Now, int valSea) {
    pdfBeamPtr->calcPDFEnvelope(idNow,xNow,Q2Now,valSea);}
  void calcPDFEnvelope(pair<int,int> idNows, pair<double,double> xNows,
    double Q2Now, int valSea) {
    pdfBeamPtr->calcPDFEnvelope(idNows,xNows,Q2Now,valSea);}
  PDF::PDFEnvelope getPDFEnvelope() { return pdfBeamPtr->getPDFEnvelope(); }

  // Decide whether chosen quark is valence, sea or companion.
  int pickValSeaComp();

  // Initialize kind of incoming beam particle.
  void initBeamKind();

  // Overload index operator to access a resolved parton from the list.
  ResolvedParton& operator[](int i) {return resolved[i];}
  const ResolvedParton& operator[](int i) const {return resolved[i];}

  // Total number of partons extracted from beam, and initiators only.
  int size() const {return resolved.size();}
  int sizeInit() const {return nInit;}

  // Clear list of resolved partons.
  void clear() {resolved.resize(0); nInit = 0;}

  // Reset variables related to photon beam.
  void resetGamma() {iGamVal = -1; iPosVal = -1; pT2gm2qqbar = 0.;
    isResolvedGamma = (gammaMode == 1) ? true : false;}

  // Reset variables related to photon beam inside a lepton.
  void resetGammaInLepton() {xGm = 1.; kTgamma = 0.; phiGamma = 0.;}

  // Add a resolved parton to list.
  int append( int iPos, int idIn, double x, int companion = -1)
    {resolved.push_back( ResolvedParton( iPos, idIn, x, companion) );
    return resolved.size() - 1;}

  // Remove the last particle from the beam. Reset companion code if needed.
  void popBack() { int iComp = resolved.back().companion();
    resolved.pop_back(); if ( iComp >= 0 ) { iSkipSave = iComp;
      idSave = resolved[iComp].id(); pickValSeaComp(); } }

  // Print extracted parton list; for debug mainly.
  void list() const;

  // How many different flavours, and how many quarks of given flavour.
  int nValenceKinds() const {return nValKinds;}
  int nValence(int idIn) const {for (int i = 0; i < nValKinds; ++i)
      if (idIn == idVal[i]) return nVal[i];
    return 0;}

  // Test whether a lepton is to be considered as unresolved.
  bool isUnresolvedLepton();

  // Add extra remnant flavours to make valence and sea come out right.
  bool remnantFlavours(Event& event, bool isDIS = false);

  // Correlate all initiators and remnants to make a colour singlet.
  bool remnantColours(Event& event, vector<int>& colFrom,
    vector<int>& colTo);

  // Pick unrescaled x of remnant parton (valence or sea).
  double xRemnant(int i);

  // Tell whether a junction has been resolved, and its junction colours.
  bool hasJunction() const {return hasJunctionBeam;}
  int junctionCol(int i) const {return junCol[i];}
  void junctionCol(int i, int col) {junCol[i] = col;}

  // For a diffractive system, decide whether to kick out gluon or quark.
  bool pickGluon(double mDiff);

  // Pick a valence quark at random, and provide the remaining flavour.
  int pickValence();
  int pickRemnant() const {return idVal2;}

  // Share lightcone momentum between two remnants in a diffractive system.
  // At the same time generate a relative pT for the two.
  double zShare( double mDiff, double m1, double m2);
  double pxShare() const {return pxRel;}
  double pyShare() const {return pyRel;}

  // Add extra remnant flavours to make valence and sea come out right.
  bool remnantFlavoursNew(Event& event);

  // Find the colour setup of the removed partons from the scatterings.
  void findColSetup(Event& event);

  // Set initial colours.
  void setInitialCol(Event & event);

  // Update colours.
  void updateCol(vector<pair<int,int> > colourChanges);

  vector<pair <int,int> > getColUpdates() {return colUpdates;}

  // Set valence content for photon beams and position of first valence quark.
  bool gammaInitiatorIsVal(int iResolved, int id, double x, double Q2);
  bool gammaInitiatorIsVal(int iResolved, double Q2);
  int  getGammaValFlavour() { return abs(idVal[0]); }
  int  gammaValSeaComp(int iResolved);
  void posVal(int iPosValIn)          { iPosVal = iPosValIn; }
  void gamVal(int iGamValIn)          { iGamVal = iGamValIn; }
  int  gamVal()                       { return iGamVal; }

  // Set and get the state (resolved and/or unresolved) of photon beam.
  void resolvedGamma(bool isResolved) { isResolvedGamma = isResolved; }
  bool resolvedGamma()                { return isResolvedGamma; }
  void setGammaMode(int gammaModeIn);
  int  getGammaMode()                 { return gammaMode; }
  bool isResolvedUnresolved()         { return isResUnres; }

  // Set state of VMD inside gamma.
  void setVMDstate(bool isVMDIn, int idIn, double mIn, double scaleIn,
    bool reassignState = false) {
    hasVMDstateInBeam = isVMDIn;
    idVMDBeam         = idIn;
    mVMDBeam          = mIn;
    scaleVMDBeam      = scaleIn;
    if (reassignState) {
      idBeam = idVMDBeam;
      mBeam  = mVMDBeam;
      pdfBeamPtr->setVMDscale(scaleVMDBeam);
    }
  }

  // Store the pT2 value of gamma->qqbar splitting.
  void   pT2gamma2qqbar(double pT2in) { pT2gm2qqbar = pT2in; }
  double pT2gamma2qqbar()             { return pT2gm2qqbar; }

  // Store the pT value for the latest MPI.
  void   pTMPI(double pTminMPIin)     { pTminMPI = pTminMPIin; }

  // Check whether room for beam remnants.
  bool roomFor1Remnant(double eCM);
  bool roomFor1Remnant(int id1, double x1, double eCM);
  bool roomFor2Remnants(int id1, double x1, double eCM);
  bool roomForRemnants(BeamParticle beamOther);

  // Evaluate the remnant mass with initiator idIn.
  double remnantMass(int idIn);

  // Functions to approximate pdfs for ISR.
  double gammaPDFxDependence(int flavour, double x)
    { return pdfBeamPtr->gammaPDFxDependence(flavour, x); }
  double gammaPDFRefScale(int flavour)
    { return pdfBeamPtr->gammaPDFRefScale(flavour); }
  double xIntegratedPDFs(double Q2)
    { return pdfBeamPtr->xfIntegratedTotal(Q2); }

  // Save the x_gamma value after latest PDF call or set it later if ND.
  void xGammaPDF()            { xGm = pdfHardBeamPtr->xGamma(); }
  void xGamma(double xGmIn)   { xGm = xGmIn; }
  void Q2Gamma(double Q2GmIn) { Q2gm = Q2GmIn; }
  void newGammaKTPhi(double kTIn, double phiIn)
    { kTgamma = kTIn; phiGamma = phiIn; }

  // Get the kinematic limits for photons emitted by the beam.
  double xGammaMin()          { return pdfHardBeamPtr->getXmin(); }
  double xGammaHadr()         { return pdfHardBeamPtr->getXhadr(); }
  double gammaFluxIntApprox() { return pdfHardBeamPtr->intFluxApprox(); }

  // Get the kinematics related photons form lepton beams.
  double xGamma()   const { return xGm; }
  double Q2Gamma()  const { return Q2gm; }
  double gammaKTx() const { return kTgamma*cos(phiGamma); }
  double gammaKTy() const { return kTgamma*sin(phiGamma); }
  double gammaKT()  const { return kTgamma; }
  double gammaPhi() const { return phiGamma; }

  // Keep track of pomeron momentum fraction.
  void xPom(double xpom = -1.0)
    { if ( pdfBeamPtr ) pdfBeamPtr->xPom(xpom); }

  // Sample x and Q2 for emitted photons according to flux.
  double sampleXgamma(double xMinIn)
    { xGm = pdfHardBeamPtr->sampleXgamma(xMinIn); return xGm; }
  double sampleQ2gamma(double Q2min)
    { Q2gm = pdfHardBeamPtr->sampleQ2gamma(Q2min); return Q2gm;}

private:

  // Constants: could only be changed in the code itself.
  static const double XMINUNRESOLVED, POMERONMASS, XMAXCOMPANION, TINYZREL;
  static const int NMAX, NRANDOMTRIES;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Pointers to PDF sets.
  PDF*          pdfBeamPtr;
  PDF*          pdfHardBeamPtr;

  // Pointer to unresolved PDF and two others to save the resolved ptrs.
  PDF*          pdfUnresBeamPtr;
  PDF*          pdfBeamPtrSave;
  PDF*          pdfHardBeamPtrSave;

  // Pointer to class for flavour generation.
  StringFlav*   flavSelPtr;

  // Initialization data, normally only set once.
  bool   allowJunction, beamJunction;
  int    maxValQuark, companionPower;
  double valencePowerMeson, valencePowerUinP, valencePowerDinP,
         valenceDiqEnhance, pickQuarkNorm, pickQuarkPower,
         diffPrimKTwidth, diffLargeMassSuppress, beamSat, gluonPower,
         xGluonCutoff;

  // Basic properties of a beam particle.
  int    idBeam, idBeamAbs, idVMDBeam;
  Vec4   pBeam;
  double mBeam, mVMDBeam, scaleVMDBeam;

  // Beam kind. Valence flavour content for hadrons.
  bool   isUnresolvedBeam, isLeptonBeam, isHadronBeam, isMesonBeam,
         isBaryonBeam, isGammaBeam;
  int    nValKinds, idVal[3], nVal[3];

  // Current parton density, by valence, sea and companion.
  int    idSave, iSkipSave, nValLeft[3];
  double xqgTot, xqVal, xqgSea, xqCompSum;

  // Variables related to photon beams (also inside lepton).
  bool   doISR, doMPI, doND, isResolvedGamma, hasResGammaInBeam,
         isResUnres, hasVMDstateInBeam;
  double pTminISR, pTminMPI, pT2gm2qqbar;
  int    iGamVal, iPosVal, gammaMode;

  // Variables for photon from lepton.
  double xGm, Q2gm, kTgamma, phiGamma;

  // The list of resolved partons.
  vector<ResolvedParton> resolved;

  // Status after all initiators have been accounted for. Junction content.
  int    nInit;
  bool   hasJunctionBeam;
  int    junCol[3];

  // Variables for new colour reconnection;
  pair <int,int> colSetup;
  vector<int> acols, cols;
  vector<bool> usedCol,usedAcol;
  vector< pair<int,int> > colUpdates;
  int nJuncs, nAjuncs, nDiffJuncs;
  bool allowBeamJunctions;

  // Routine to calculate pdf's given previous interactions.
  double xfModified( int iSkip, int idIn, double x, double Q2);

  // Fraction of hadron momentum sitting in a valence quark distribution.
  double xValFrac(int j, double Q2);
  double Q2ValFracSav, uValInt, dValInt;

  // Fraction of hadron momentum sitting in a companion quark distribution.
  double xCompFrac(double xs);

  // Value of companion quark PDF, also given the sea quark x.
  double xCompDist(double xc, double xs);

  // Valence quark subdivision for diffractive systems.
  int    idVal1, idVal2, idVal3;
  double zRel, pxRel, pyRel;

  // Update a single (anti-) colour of the event.
  void updateSingleCol(int oldCol, int newCol);

  // Find a single (anti-) colour in the beam,
  // if a beam remnant is set the new colour.
  int findSingleCol(Event& event, bool isAcol, bool useHardScatters);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_BeamParticle_H
