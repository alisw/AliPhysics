// MultipartonInteractions.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main classes for multiparton interactions physics.
// SigmaMultiparton stores allowed processes by in-flavour combination.
// MultipartonInteractions: generates multiparton parton-parton interactions.

#ifndef Pythia8_MultipartonInteractions_H
#define Pythia8_MultipartonInteractions_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PartonVertex.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/SigmaTotal.h"
#include "Pythia8/SigmaProcess.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/UserHooks.h"

namespace Pythia8 {

//==========================================================================

// SigmaMultiparton is a helper class to MultipartonInteractions.
// It packs pointers to the allowed processes for different
// flavour combinations and levels of ambition.

class SigmaMultiparton {

public:

  // Constructor.
  SigmaMultiparton() : nChan(), needMasses(), useNarrowBW3(), useNarrowBW4(),
    m3Fix(), m4Fix(), sHatMin(), sigmaT(), sigmaU(), sigmaTval(), sigmaUval(),
    sigmaTsum(), sigmaUsum(), pickOther(), pickedU(), particleDataPtr(),
    rndmPtr() {}

  // Destructor.
  ~SigmaMultiparton() {
    for (int i = 0; i < int(sigmaT.size()); ++i) delete sigmaT[i];
    for (int i = 0; i < int(sigmaU.size()); ++i) delete sigmaU[i];}

  // Initialize list of processes.
  bool init(int inState, int processLevel, Info* infoPtr,
    Settings* settingsPtr, ParticleData* particleDataPtr, Rndm* rndmPtrIn,
    BeamParticle* beamAPtr, BeamParticle* beamBPtr, Couplings* couplingsPtr);

  // Calculate cross section summed over possibilities.
  double sigma( int id1, int id2, double x1, double x2, double sHat,
    double tHat, double uHat, double alpS, double alpEM,
    bool restore = false, bool pickOtherIn = false);

  // Return whether the other, rare processes were selected.
  bool pickedOther() {return pickOther;}

  // Return one subprocess, picked according to relative cross sections.
  SigmaProcess* sigmaSel();
  bool swapTU() {return pickedU;}

  // Return code or name of a specified process, for statistics table.
  int    nProc() const {return nChan;}
  int    codeProc(int iProc) const {return sigmaT[iProc]->code();}
  string nameProc(int iProc) const {return sigmaT[iProc]->name();}

private:

  // Constants: could only be changed in the code itself.
  static const double MASSMARGIN, OTHERFRAC;

  // Number of processes. Some use massive matrix elements.
  int            nChan;
  vector<bool>   needMasses, useNarrowBW3, useNarrowBW4;
  vector<double> m3Fix, m4Fix, sHatMin;

  // Vector of process list, one for t-channel and one for u-channel.
  vector<SigmaProcess*> sigmaT, sigmaU;

  // Values of cross sections in process list above.
  vector<double> sigmaTval, sigmaUval;
  double         sigmaTsum, sigmaUsum;
  bool           pickOther, pickedU;

  // Pointers to the particle data and random number generator.
  ParticleData*  particleDataPtr;
  Rndm*          rndmPtr;

};

//==========================================================================

// The MultipartonInteractions class contains the main methods for the
// generation of multiparton parton-parton interactions in hadronic collisions.

class MultipartonInteractions {

public:

  // Constructor.
  MultipartonInteractions() : allowRescatter(), allowDoubleRes(), canVetoMPI(),
    doPartonVertex(), doVarEcm(), pTmaxMatch(), alphaSorder(), alphaEMorder(),
    alphaSnfmax(), bProfile(), processLevel(), bSelScale(), rescatterMode(),
    nQuarkIn(), nSample(), enhanceScreening(), pT0paramMode(), alphaSvalue(),
    Kfactor(), pT0Ref(), ecmRef(), ecmPow(), pTmin(), coreRadius(),
    coreFraction(), expPow(), ySepResc(), deltaYResc(), sigmaPomP(), mPomP(),
    pPomP(), mMaxPertDiff(), mMinPertDiff(), a1(), a0now(), a02now(),
    bstepNow(), a2max(), b2now(), enhanceBmax(), enhanceBnow(), id1Save(),
    id2Save(), pT2Save(), x1Save(), x2Save(), sHatSave(), tHatSave(),
    uHatSave(), alpSsave(), alpEMsave(), pT2FacSave(), pT2RenSave(),
    xPDF1nowSave(), xPDF2nowSave(), scaleLimitPTsave(), dSigmaDtSelSave(),
    vsc1(), vsc2(), hasBaryonBeams(), hasPomeronBeams(), hasLowPow(),
    globalRecoilFSR(), iDiffSys(), nMaxGlobalRecoilFSR(), bSelHard(), eCM(),
    sCM(), pT0(), pT20(), pT2min(), pTmax(), pT2max(), pT20R(), pT20minR(),
    pT20maxR(), pT20min0maxR(), pT2maxmin(), sigmaND(), pT4dSigmaMax(),
    pT4dProbMax(), dSigmaApprox(), sigmaInt(), sudExpPT(), zeroIntCorr(),
    normOverlap(), nAvg(), kNow(), normPi(), bAvg(), bDiv(), probLowB(),
    radius2B(), radius2C(), fracA(), fracB(), fracC(), fracAhigh(),
    fracBhigh(), fracChigh(), fracABChigh(), expRev(), cDiv(), cMax(),
    enhanceBavg(), bIsSet(false), bSetInFirst(), isAtLowB(), pickOtherSel(),
    id1(), id2(), i1Sel(), i2Sel(), id1Sel(), id2Sel(), bNow(), enhanceB(),
    pT2(), pT2shift(), pT2Ren(), pT2Fac(), x1(), x2(), xT(), xT2(), tau(),
    y(), sHat(), tHat(), uHat(), alpS(), alpEM(), xPDF1now(), xPDF2now(),
    dSigmaSum(), x1Sel(), x2Sel(), sHatSel(), tHatSel(), uHatSel(), nStep(),
    iStepFrom(), iStepTo(), eCMsave(), eStepMin(), eStepMax(), eStepSize(),
    eStepSave(), eStepFrom(), eStepTo(), pT0Save(), pT4dSigmaMaxSave(),
    pT4dProbMaxSave(), sigmaIntSave(), sudExpPTSave(), zeroIntCorrSave(),
    normOverlapSave(), kNowSave(), bAvgSave(), bDivSave(), probLowBSave(),
    fracAhighSave(), fracBhighSave(), fracChighSave(), fracABChighSave(),
    cDivSave(), cMaxSave(), beamOffset(), mGmGmMin(), mGmGmMax(), hasGamma(),
    isGammaGamma(), isGammaHadron(), isHadronGamma(), infoPtr(), rndmPtr(),
    beamAPtr(), beamBPtr(), couplingsPtr(), partonSystemsPtr(), sigmaTotPtr(),
    userHooksPtr(), partonVertexPtr(), sigma2Sel(), dSigmaDtSel() {}

  // Initialize the generation process for given beams.
  bool init( bool doMPIinit, int iDiffSysIn, Info* infoPtrIn,
    Settings& settings, ParticleData* particleDataPtr, Rndm* rndmPtrIn,
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    Couplings* couplingsPtrIn, PartonSystems* partonSystemsPtrIn,
    SigmaTotal* sigmaTotPtrIn, UserHooks* userHooksPtrIn,
    PartonVertex* partonVertexPtrIn, bool hasGammaIn = false);

  // Reset impact parameter choice and update the CM energy.
  void reset();

  // Select first = hardest pT in nondiffractive process.
  void pTfirst();

  // Set up kinematics for first = hardest pT in nondiffractive process.
  void setupFirstSys( Event& process);

  // Find whether to limit maximum scale of emissions.
  // Provide sum pT / 2 as potential limit where relevant.
  bool limitPTmax( Event& event);
  double scaleLimitPT() const {return scaleLimitPTsave;}

  // Prepare system for evolution.
  void prepare(Event& event, double pTscale = 1000., bool rehashB = false) {
    if (!bSetInFirst) overlapNext(event, pTscale, rehashB); }

  // Select next pT in downwards evolution.
  double pTnext( double pTbegAll, double pTendAll, Event& event);

  // Set up kinematics of acceptable interaction.
  bool scatter( Event& event);

  // Set "empty" values to avoid query of undefined quantities.
  void setEmpty() {pT2Ren = alpS = alpEM = x1 = x2 = pT2Fac
    = xPDF1now = xPDF2now = 0.; bIsSet = false;}

  // Get some information on current interaction.
  double Q2Ren()      const {return pT2Ren;}
  double alphaSH()    const {return alpS;}
  double alphaEMH()   const {return alpEM;}
  double x1H()        const {return x1;}
  double x2H()        const {return x2;}
  double Q2Fac()      const {return pT2Fac;}
  double pdf1()       const {return (id1 == 21 ? 4./9. : 1.) * xPDF1now;}
  double pdf2()       const {return (id2 == 21 ? 4./9. : 1.) * xPDF2now;}
  double bMPI()       const {return (bIsSet) ? bNow : 0.;}
  double enhanceMPI() const {return (bIsSet) ? enhanceB / zeroIntCorr : 1.;}
  double enhanceMPIavg() const {return enhanceBavg;}

  // For x-dependent matter profile, return incoming valence/sea
  // decision from trial interactions.
  int    getVSC1()   const {return vsc1;}
  int    getVSC2()   const {return vsc2;}

  // Set the offset wrt. to normal beam particle positions for hard diffraction
  // and for photon beam from lepton.
  int  getBeamOffset()       const {return beamOffset;}
  void setBeamOffset(int offsetIn) {beamOffset = offsetIn;}

  // Update and print statistics on number of processes.
  // Note: currently only valid for nondiffractive systems, not diffraction??
  void accumulate() { int iBeg = (infoPtr->isNonDiffractive()) ? 0 : 1;
    for (int i = iBeg; i < infoPtr->nMPI(); ++i)
    ++nGen[ infoPtr->codeMPI(i) ];}
  void statistics(bool resetStat = false);
  void resetStatistics() { for ( map<int, int>::iterator iter = nGen.begin();
    iter != nGen.end(); ++iter) iter->second = 0; }

private:

  // Constants: could only be changed in the code itself.
  static const bool   SHIFTFACSCALE, PREPICKRESCATTER;
  static const double SIGMAFUDGE, RPT20, PT0STEP, SIGMASTEP, PT0MIN,
                      EXPPOWMIN, PROBATLOWB, BSTEP, BMAX, EXPMAX,
                      KCONVERGE, CONVERT2MB, ROOTMIN, ECMDEV, WTACCWARN,
                      SIGMAMBLIMIT;

  // Initialization data, read from Settings.
  bool   allowRescatter, allowDoubleRes, canVetoMPI, doPartonVertex, doVarEcm;
  int    pTmaxMatch, alphaSorder, alphaEMorder, alphaSnfmax, bProfile,
         processLevel, bSelScale, rescatterMode, nQuarkIn, nSample,
         enhanceScreening, pT0paramMode;
  double alphaSvalue, Kfactor, pT0Ref, ecmRef, ecmPow, pTmin, coreRadius,
         coreFraction, expPow, ySepResc, deltaYResc, sigmaPomP, mPomP, pPomP,
         mMaxPertDiff, mMinPertDiff;

  // x-dependent matter profile:
  // Constants.
  static const int    XDEP_BBIN;
  static const double XDEP_A0, XDEP_A1, XDEP_BSTEP, XDEP_BSTEPINC,
                      XDEP_CUTOFF, XDEP_WARN, XDEP_SMB2FM;

  // Table of Int( dSigma/dX * O(b, X), dX ) in bins of b for fast
  // integration. Only needed during initialisation.
  vector <double> sigmaIntWgt, sigmaSumWgt;

  // a1:             value of a1 constant, taken from settings database.
  // a0now (a02now): tuned value of a0 (squared value).
  // bstepNow:       current size of bins in b.
  // a2max:          given an xMin, a maximal (squared) value of the
  //                 width, to be used in overestimate Omax(b).
  // enhanceBmax,    retain enhanceB as enhancement factor for the hardest
  // enhanceBnow:    interaction. Use enhanceBmax as overestimate for fastPT2,
  //                 and enhanceBnow to store the value for the current
  //                 interaction.
  double a1, a0now, a02now, bstepNow, a2max, b2now, enhanceBmax, enhanceBnow;

  // Storage for trial interactions.
  int    id1Save, id2Save;
  double pT2Save, x1Save, x2Save, sHatSave, tHatSave, uHatSave,
         alpSsave, alpEMsave, pT2FacSave, pT2RenSave, xPDF1nowSave,
         xPDF2nowSave, scaleLimitPTsave;
  SigmaProcess *dSigmaDtSelSave;

  // vsc1, vsc2:     for minimum-bias events with trial interaction, store
  //                 decision on whether hard interaction was valence or sea.
  int    vsc1, vsc2;

  // Other initialization data.
  bool   hasBaryonBeams, hasPomeronBeams, hasLowPow, globalRecoilFSR;
  int    iDiffSys, nMaxGlobalRecoilFSR, bSelHard;
  double eCM, sCM, pT0, pT20, pT2min, pTmax, pT2max, pT20R, pT20minR,
         pT20maxR, pT20min0maxR, pT2maxmin, sigmaND, pT4dSigmaMax,
         pT4dProbMax, dSigmaApprox, sigmaInt, sudExpPT[101],
         zeroIntCorr, normOverlap, nAvg, kNow, normPi, bAvg, bDiv,
         probLowB, radius2B, radius2C, fracA, fracB, fracC, fracAhigh,
         fracBhigh, fracChigh, fracABChigh, expRev, cDiv, cMax,
         enhanceBavg;

  // Properties specific to current system.
  bool   bIsSet, bSetInFirst, isAtLowB, pickOtherSel;
  int    id1, id2, i1Sel, i2Sel, id1Sel, id2Sel;
  double bNow, enhanceB, pT2, pT2shift, pT2Ren, pT2Fac, x1, x2, xT, xT2,
         tau, y, sHat, tHat, uHat, alpS, alpEM, xPDF1now, xPDF2now,
         dSigmaSum, x1Sel, x2Sel, sHatSel, tHatSel, uHatSel;

  // Stored values for mass interpolation for diffractive systems.
  int    nStep, iStepFrom, iStepTo;
  double eCMsave, eStepMin, eStepMax, eStepSize, eStepSave, eStepFrom, eStepTo,
         pT0Save[20], pT4dSigmaMaxSave[20], pT4dProbMaxSave[20],
         sigmaIntSave[20], sudExpPTSave[20][101], zeroIntCorrSave[20],
         normOverlapSave[20], kNowSave[20], bAvgSave[20], bDivSave[20],
         probLowBSave[20], fracAhighSave[20], fracBhighSave[20],
         fracChighSave[20], fracABChighSave[20], cDivSave[20], cMaxSave[20];

  // Beam offset wrt. normal situation and other photon-related parameters.
  int    beamOffset;
  double mGmGmMin, mGmGmMax;
  bool   hasGamma, isGammaGamma, isGammaHadron, isHadronGamma;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to the two incoming beams.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;

  // Pointers to Standard Model couplings.
  Couplings*     couplingsPtr;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Pointer to total cross section parametrization.
  SigmaTotal*    sigmaTotPtr;

  // Pointer to user hooks.
  UserHooks*     userHooksPtr;

  // Pointer to assign space-time vertices during parton evolution.
  PartonVertex*  partonVertexPtr;

  // Collections of parton-level 2 -> 2 cross sections. Selected one.
  SigmaMultiparton  sigma2gg, sigma2qg, sigma2qqbarSame, sigma2qq;
  SigmaMultiparton* sigma2Sel;
  SigmaProcess*  dSigmaDtSel;

  // Statistics on generated 2 -> 2 processes.
  map<int, int>  nGen;

  // alphaStrong and alphaEM calculations.
  AlphaStrong    alphaS;
  AlphaEM        alphaEM;

  // Scattered partons.
  vector<int>    scatteredA, scatteredB;

  // Determine constant in d(Prob)/d(pT2) < const / (pT2 + r * pT20)^2.
  void upperEnvelope();

  // Integrate the parton-parton interaction cross section.
  void jetCrossSection();

  // Evaluate "Sudakov form factor" for not having a harder interaction.
  double sudakov(double pT2sud, double enhance = 1.);

  // Do a quick evolution towards the next smaller pT.
  double fastPT2( double pT2beg);

  // Calculate the actual cross section, either for the first interaction
  // (including at initialization) or for any subsequent in the sequence.
  double sigmaPT2scatter(bool isFirst = false);

  // Find the partons that may rescatter.
  void findScatteredPartons( Event& event);

  // Calculate the actual cross section for a rescattering.
  double sigmaPT2rescatter( Event& event);

  // Calculate factor relating matter overlap and interaction rate.
  void overlapInit();

  // Pick impact parameter and interaction rate enhancement,
  // either before the first interaction (for nondiffractive) or after it.
  void overlapFirst();
  void overlapNext(Event& event, double pTscale, bool rehashB);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_MultipartonInteractions_H
