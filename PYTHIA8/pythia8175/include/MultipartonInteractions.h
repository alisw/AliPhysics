// MultipartonInteractions.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main classes for multiparton interactions physics.
// SigmaMultiparton stores allowed processes by in-flavour combination.
// MultipartonInteractions: generates multiparton parton-parton interactions.

#ifndef Pythia8_MultipartonInteractions_H
#define Pythia8_MultipartonInteractions_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "Info.h"
#include "PartonSystems.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "SigmaTotal.h"
#include "SigmaProcess.h"
#include "StandardModel.h"
#include "UserHooks.h"

namespace Pythia8 {
 
//==========================================================================

// SigmaMultiparton is a helper class to MultipartonInteractions.
// It packs pointers to the allowed processes for different 
// flavour combinations and levels of ambition.

class SigmaMultiparton {

public:

  // Constructor.
  SigmaMultiparton() : nChan(0),
    sigmaTsum(0), sigmaUsum(0),
    pickOther(0), pickedU(0),
    rndmPtr(0x0) {}
  
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
  vector<bool>   needMasses;
  vector<double> m3Fix, m4Fix, sHatMin;

  // Vector of process list, one for t-channel and one for u-channel.
  vector<SigmaProcess*> sigmaT, sigmaU;

  // Values of cross sections in process list above.
  vector<double> sigmaTval, sigmaUval;
  double         sigmaTsum, sigmaUsum;
  bool           pickOther, pickedU;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;
  
};
 
//==========================================================================

// The MultipartonInteractions class contains the main methods for the 
// generation of multiparton parton-parton interactions in hadronic collisions.

class MultipartonInteractions {

public:

  // Constructor.
 MultipartonInteractions() :
  allowRescatter(false), allowDoubleRes(false), canVetoMPI(false),
  pTmaxMatch(0), alphaSorder(0), alphaEMorder(0), bProfile(0), processLevel(0), 
  bSelScale(0), rescatterMode(0), nQuarkIn(0), nSample(0), enhanceScreening(0),
  alphaSvalue(0), Kfactor(0), pT0Ref(0), ecmRef(0), ecmPow(0), pTmin(0), coreRadius(0), 
  coreFraction(0), expPow(0), ySepResc(0), deltaYResc(0), sigmaPomP(0), mPomP(0), pPomP(0), 
  mMaxPertDiff(0), mMinPertDiff(0),
  a1(0), a0now(0), a02now(0), bstepNow(0), a2max(0), b2now(0), enhanceBmax(0), enhanceBnow(0),
  id1Save(0), id2Save(0),
  pT2Save(0), x1Save(0), x2Save(0), sHatSave(0), tHatSave(0), uHatSave(0),
  alpSsave(0), alpEMsave(0), pT2FacSave(0), pT2RenSave(0), xPDF1nowSave(0),
  xPDF2nowSave(0),
  dSigmaDtSelSave(0x0),
  vsc1(0), vsc2(0),
  hasBaryonBeams(false), hasLowPow(false), globalRecoilFSR(false),
  iDiffSys(0), nMaxGlobalRecoilFSR(0),
  eCM(0), sCM(0), pT0(0), pT20(0), pT2min(0), pTmax(0), pT2max(0), pT20R(0), pT20minR(0), 
  pT20maxR(0), pT20min0maxR(0), pT2maxmin(0), sigmaND(0), pT4dSigmaMax(0), 
  pT4dProbMax(0), dSigmaApprox(0), sigmaInt(0),
  zeroIntCorr(0), normOverlap(0), nAvg(0), kNow(0), normPi(0), bAvg(0), bDiv(0), 
  probLowB(0), radius2B(0), radius2C(0), fracA(0), fracB(0), fracC(0), fracAhigh(0), 
  fracBhigh(0), fracChigh(0), fracABChigh(0), expRev(0), cDiv(0), cMax(0),
  bIsSet(false), bSetInFirst(false), isAtLowB(false), pickOtherSel(false),
  id1(0), id2(0), i1Sel(0), i2Sel(0), id1Sel(0), id2Sel(0),
  bNow(0), enhanceB(0), pT2(0), pT2shift(0), pT2Ren(0), pT2Fac(0), x1(0), x2(0), xT(0), xT2(0), 
  tau(0), y(0), sHat(0), tHat(0), uHat(0), alpS(0), alpEM(0), xPDF1now(0), xPDF2now(0),
  dSigmaSum(0), x1Sel(0), x2Sel(0), sHatSel(0), tHatSel(0), uHatSel(0),
  nStep(0), iStepFrom(0), iStepTo(0), 
  eCMsave(0), eStepSize(0), eStepSave(0), eStepFrom(0), eStepTo(0),
  infoPtr(0x0),
  rndmPtr(0x0),
  beamAPtr(0x0),
  beamBPtr(0x0),
  couplingsPtr(0x0),
  partonSystemsPtr(0x0),
  sigmaTotPtr(0x0),
  userHooksPtr(0x0),
  sigma2gg(), sigma2qg(), sigma2qqbarSame(), sigma2qq(),
  sigma2Sel(0x0),
  dSigmaDtSel(0x0),
  alphaS(),
  alphaEM() {
    for (int i=0; i<101; i++)  sudExpPT[i] = 0;
    for (int i=0; i<5; i++) {
      pT0Save[i] = 0; 
      pT4dSigmaMaxSave[i] = 0; pT4dProbMaxSave[i] = 0; sigmaIntSave[i] = 0; 
      
      for (int j=0; j<101; j++) sudExpPTSave[i][j] = 0;
      zeroIntCorrSave[i] = 0; normOverlapSave[i] = 0; 
      kNowSave[i] = 0; bAvgSave[i] = 0; bDivSave[i] = 0; probLowBSave[i] = 0; 
      fracAhighSave[i] = 0; fracBhighSave[i] = 0; fracChighSave[i] = 0; 
      fracABChighSave[i] = 0; cDivSave[i] = 0; cMaxSave[i] = 0;
    }
  }

  // Initialize the generation process for given beams.
  bool init( bool doMPIinit, int iDiffSysIn, Info* infoPtrIn, 
    Settings& settings, ParticleData* particleDataPtr, Rndm* rndmPtrIn, 
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, 
    Couplings* couplingsPtrIn, PartonSystems* partonSystemsPtrIn, 
    SigmaTotal* sigmaTotPtrIn, UserHooks* userHooksPtrIn,
    ostream& os = cout);

  // Reset impact parameter choice and update the CM energy.
  void reset();

  // Select first = hardest pT in minbias process.
  void pTfirst(); 

  // Set up kinematics for first = hardest pT in minbias process.
  void setupFirstSys( Event& process);

  // Find whether to limit maximum scale of emissions.
  bool limitPTmax( Event& event);

  // Prepare system for evolution.
  void prepare(Event& event, double pTscale = 1000.) {
    if (!bSetInFirst) overlapNext(event, pTscale); }

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
  double pdf1()       const {return xPDF1now;}
  double pdf2()       const {return xPDF2now;}
  double bMPI()       const {return (bIsSet) ? bNow / bAvg : 0.;}
  double enhanceMPI() const {return (bIsSet) ? enhanceB / zeroIntCorr : 1.;}

  // For x-dependent matter profile, return incoming valence/sea
  // decision from trial interactions.
  int    getVSC1()   const {return vsc1;}
  int    getVSC2()   const {return vsc2;}

  // Update and print statistics on number of processes.
  // Note: currently only valid for MB systems, not diffraction??
  void accumulate() { int iBeg = (infoPtr->isMinBias()) ? 0 : 1; 
    for (int i = iBeg; i < infoPtr->nMPI(); ++i) 
    ++nGen[ infoPtr->codeMPI(i) ];}
  void statistics(bool resetStat = false, ostream& os = cout);
  void resetStatistics() { for ( map<int, int>::iterator iter = nGen.begin(); 
    iter != nGen.end(); ++iter) iter->second = 0; } 
  
private: 

  // Constants: could only be changed in the code itself.
  static const bool   SHIFTFACSCALE, PREPICKRESCATTER;
  static const double SIGMAFUDGE, RPT20, PT0STEP, SIGMASTEP, PT0MIN,
                      EXPPOWMIN, PROBATLOWB, BSTEP, BMAX, EXPMAX, 
                      KCONVERGE, CONVERT2MB, ROOTMIN, ECMDEV, WTACCWARN;

  // Initialization data, read from Settings.
  bool   allowRescatter, allowDoubleRes, canVetoMPI;
  int    pTmaxMatch, alphaSorder, alphaEMorder, bProfile, processLevel, 
         bSelScale, rescatterMode, nQuarkIn, nSample, enhanceScreening;
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
         xPDF2nowSave;
  SigmaProcess *dSigmaDtSelSave;

  // vsc1, vsc2:     for minimum-bias events with trial interaction, store
  //                 decision on whether hard interaction was valence or sea.
  int    vsc1, vsc2;

  // Other initialization data.
  bool   hasBaryonBeams, hasLowPow, globalRecoilFSR;
  int    iDiffSys, nMaxGlobalRecoilFSR;
  double eCM, sCM, pT0, pT20, pT2min, pTmax, pT2max, pT20R, pT20minR, 
         pT20maxR, pT20min0maxR, pT2maxmin, sigmaND, pT4dSigmaMax, 
         pT4dProbMax, dSigmaApprox, sigmaInt, sudExpPT[101], 
         zeroIntCorr, normOverlap, nAvg, kNow, normPi, bAvg, bDiv, 
         probLowB, radius2B, radius2C, fracA, fracB, fracC, fracAhigh, 
         fracBhigh, fracChigh, fracABChigh, expRev, cDiv, cMax;

  // Properties specific to current system.
  bool   bIsSet, bSetInFirst, isAtLowB, pickOtherSel;
  int    id1, id2, i1Sel, i2Sel, id1Sel, id2Sel;
  double bNow, enhanceB, pT2, pT2shift, pT2Ren, pT2Fac, x1, x2, xT, xT2, 
         tau, y, sHat, tHat, uHat, alpS, alpEM, xPDF1now, xPDF2now,
         dSigmaSum, x1Sel, x2Sel, sHatSel, tHatSel, uHatSel;

  // Stored values for mass interpolation for diffractive systems.
  int    nStep, iStepFrom, iStepTo; 
  double eCMsave, eStepSize, eStepSave, eStepFrom, eStepTo, pT0Save[5], 
         pT4dSigmaMaxSave[5], pT4dProbMaxSave[5], sigmaIntSave[5], 
         sudExpPTSave[5][101], zeroIntCorrSave[5], normOverlapSave[5], 
         kNowSave[5], bAvgSave[5], bDivSave[5], probLowBSave[5], 
         fracAhighSave[5], fracBhighSave[5], fracChighSave[5], 
         fracABChighSave[5], cDivSave[5], cMaxSave[5];

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
  // either before the first interaction (for minbias) or after it.
  void overlapFirst();
  void overlapNext(Event& event, double pTscale);

};
 
//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_MultipartonInteractions_H
