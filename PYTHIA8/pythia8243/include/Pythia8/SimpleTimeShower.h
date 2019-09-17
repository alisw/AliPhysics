// SimpleTimeShower.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the original simple timelike final-state showers.
// TimeDipoleEnd: data on a radiating dipole end in FSR.
// SimpleTimeShower: handles the showering description.

#ifndef Pythia8_SimpleTimeShower_H
#define Pythia8_SimpleTimeShower_H

#include "Pythia8/TimeShower.h"
#include "Pythia8/SimpleWeakShowerMEs.h"

namespace Pythia8 {

//==========================================================================

// Data on radiating dipole ends; only used inside SimpleTimeShower class.

class TimeDipoleEnd {

public:

  // Constructors.
  TimeDipoleEnd() : iRadiator(-1), iRecoiler(-1), pTmax(0.), colType(0),
    chgType(0), gamType(0), weakType(0), isrType(0), system(0), systemRec(0),
    MEtype(0), iMEpartner(-1), weakPol(0), isOctetOnium(false),
    isHiddenValley(false), colvType(0), MEmix(0.), MEorder(true),
    MEsplit(true), MEgluinoRec(false), isFlexible(false), flavour(), iAunt(),
    mRad(), m2Rad(), mRec(), m2Rec(), mDip(), m2Dip(), m2DipCorr(), pT2(),
    m2(), z(), mFlavour(), asymPol(), flexFactor(), pAccept() { }
  TimeDipoleEnd(int iRadiatorIn, int iRecoilerIn, double pTmaxIn = 0.,
    int colIn = 0, int chgIn = 0, int gamIn = 0, int weakTypeIn = 0,
    int isrIn = 0, int systemIn = 0, int MEtypeIn = 0, int iMEpartnerIn = -1,
    int weakPolIn = 0, bool isOctetOniumIn = false,
    bool isHiddenValleyIn = false, int colvTypeIn = 0, double MEmixIn = 0.,
    bool MEorderIn = true, bool MEsplitIn = true, bool MEgluinoRecIn = false,
    bool isFlexibleIn = false) :
    iRadiator(iRadiatorIn), iRecoiler(iRecoilerIn), pTmax(pTmaxIn),
    colType(colIn), chgType(chgIn), gamType(gamIn), weakType(weakTypeIn),
    isrType(isrIn), system(systemIn), systemRec(systemIn), MEtype(MEtypeIn),
    iMEpartner(iMEpartnerIn), weakPol(weakPolIn), isOctetOnium(isOctetOniumIn),
    isHiddenValley(isHiddenValleyIn), colvType(colvTypeIn), MEmix(MEmixIn),
    MEorder (MEorderIn), MEsplit(MEsplitIn), MEgluinoRec(MEgluinoRecIn),
    isFlexible(isFlexibleIn), flavour(), iAunt(), mRad(), m2Rad(), mRec(),
    m2Rec(), mDip(), m2Dip(), m2DipCorr(), pT2(), m2(), z(), mFlavour(),
    asymPol(), flexFactor(), pAccept()  { }

  // Basic properties related to dipole and matrix element corrections.
  int    iRadiator, iRecoiler;
  double pTmax;
  int    colType, chgType, gamType, weakType, isrType, system, systemRec,
         MEtype, iMEpartner, weakPol;
  bool   isOctetOnium, isHiddenValley;
  int    colvType;
  double MEmix;
  bool   MEorder, MEsplit, MEgluinoRec, isFlexible;

  // Properties specific to current trial emission.
  int    flavour, iAunt;
  double mRad, m2Rad, mRec, m2Rec, mDip, m2Dip, m2DipCorr,
         pT2, m2, z, mFlavour, asymPol, flexFactor, pAccept;

};

//==========================================================================

// The SimpleTimeShower class does timelike showers.

class SimpleTimeShower : public TimeShower {

public:

  // Constructor.
  SimpleTimeShower() : hasWeaklyRadiated(), iSysSel(), pTmaxFudge(),
    pTLastBranch(), doQCDshower(), doQEDshowerByQ(), doQEDshowerByL(),
    doQEDshowerByOther(), doQEDshowerByGamma(), doWeakShower(),
    doMEcorrections(), doMEextended(), doMEafterFirst(), doPhiPolAsym(),
    doPhiPolAsymHard(), doInterleave(), allowBeamRecoil(), dampenBeamRecoil(),
    recoilToColoured(), useFixedFacScale(), allowRescatter(),
    canVetoEmission(), doHVshower(), brokenHVsym(), globalRecoil(),
    useLocalRecoilNow(), doSecondHard(), hasUserHooks(), singleWeakEmission(),
    alphaSuseCMW(), vetoWeakJets(), allowMPIdipole(), weakExternal(),
    recoilDeadCone(), doDipoleRecoil(), doPartonVertex(), pTmaxMatch(),
    pTdampMatch(), alphaSorder(), alphaSnfmax(), nGluonToQuark(),
    weightGluonToQuark(), alphaEMorder(), nGammaToQuark(), nGammaToLepton(),
    nCHV(), idHV(), alphaHVorder(), nMaxGlobalRecoil(), weakMode(),
    pTdampFudge(), mc(), mb(), m2c(), m2b(), renormMultFac(), factorMultFac(),
    fixedFacScale2(), alphaSvalue(), alphaS2pi(), Lambda3flav(), Lambda4flav(),
    Lambda5flav(), Lambda3flav2(), Lambda4flav2(), Lambda5flav2(),
    scaleGluonToQuark(), extraGluonToQuark(), pTcolCutMin(), pTcolCut(),
    pT2colCut(), pTchgQCut(), pT2chgQCut(), pTchgLCut(), pT2chgLCut(),
    pTweakCut(), pT2weakCut(), mMaxGamma(), m2MaxGamma(), octetOniumFraction(),
    octetOniumColFac(), mZ(), gammaZ(), thetaWRat(), mW(), gammaW(), CFHV(),
    nFlavHV(), alphaHVfix(), LambdaHV(), pThvCut(), pT2hvCut(), mHV(),
    pTmaxFudgeMPI(), weakEnhancement(), vetoWeakDeltaR2(), twoHard(),
    dopTlimit1(), dopTlimit2(), dopTdamp(), pT2damp(), kRad(), kEmt(),
    pdfScale2(), doTrialNow(), canEnhanceEmission(), canEnhanceTrial(),
    canEnhanceET(), doUncertaintiesNow(), dipSel(), iDipSel(), nHard(),
    nFinalBorn(), nMaxGlobalBranch(), nGlobal(), globalRecoilMode(),
    limitMUQ(), weakHardSize() { beamOffset = 0;}

  // Destructor.
  virtual ~SimpleTimeShower() {}

  // Initialize alphaStrong and related pTmin parameters.
  virtual void init( BeamParticle* beamAPtrIn = 0,
    BeamParticle* beamBPtrIn = 0);

  // Find whether to limit maximum scale of emissions, and whether to dampen.
  virtual bool limitPTmax( Event& event, double Q2Fac = 0.,
    double Q2Ren = 0.);

  // Top-level routine to do a full time-like shower in resonance decay.
  virtual int shower( int iBeg, int iEnd, Event& event, double pTmax,
    int nBranchMax = 0);

  // Top-level routine for QED radiation in hadronic decay to two leptons.
  virtual int showerQED( int i1, int i2, Event& event, double pTmax);

  // Global recoil: reset counters and store locations of outgoing partons.
  virtual void prepareGlobal( Event& event);

  // Prepare system for evolution after each new interaction; identify ME.
  virtual void prepare( int iSys, Event& event, bool limitPTmaxIn = true);

  // Update dipole list after a multiparton interactions rescattering.
  virtual void rescatterUpdate( int iSys, Event& event);

  // Update dipole list after each ISR emission.
  virtual void update( int iSys, Event& event, bool hasWeakRad = false);

  // Select next pT in downwards evolution.
  virtual double pTnext( Event& event, double pTbegAll, double pTendAll,
    bool isFirstTrial = false, bool doTrialIn = false);

  // ME corrections and kinematics that may give failure.
  virtual bool branch( Event& event, bool isInterleaved = false);

  // Print dipole list; for debug mainly.
  virtual void list() const;

  // Initialize data members for calculation of uncertainty bands.
  virtual bool initUncertainties();

  // Tell whether FSR has done a weak emission.
  virtual bool getHasWeaklyRadiated() {return hasWeaklyRadiated;}

  // Tell which system was the last processed one.
  virtual int system() const {return iSysSel;}

  // Potential enhancement factor of pTmax scale for hardest emission.
  virtual double enhancePTmax() {return pTmaxFudge;}

  // Provide the pT scale of the last branching in the above shower.
  virtual double pTLastInShower() {return pTLastBranch;}

private:

  // Constants: could only be changed in the code itself.
  static const double MCMIN, MBMIN, SIMPLIFYROOT, XMARGIN, XMARGINCOMB,
         TINYPDF, LARGEM2, THRESHM2, LAMBDA3MARGIN, WEAKPSWEIGHT, WG2QEXTRA,
         REJECTFACTOR, PROBLIMIT;
  // Rescatter: try to fix up recoil between systems
  static const bool   FIXRESCATTER, VETONEGENERGY;
  static const double MAXVIRTUALITYFRACTION, MAXNEGENERGYFRACTION;

  // Store properties to be returned by methods.
  bool   hasWeaklyRadiated;
  int    iSysSel;
  double pTmaxFudge, pTLastBranch;

  // Initialization data, normally only set once.
  bool   doQCDshower, doQEDshowerByQ, doQEDshowerByL, doQEDshowerByOther,
         doQEDshowerByGamma, doWeakShower, doMEcorrections, doMEextended,
         doMEafterFirst, doPhiPolAsym, doPhiPolAsymHard, doInterleave,
         allowBeamRecoil, dampenBeamRecoil, recoilToColoured, useFixedFacScale,
         allowRescatter, canVetoEmission, doHVshower, brokenHVsym,
         globalRecoil, useLocalRecoilNow, doSecondHard, hasUserHooks,
         singleWeakEmission, alphaSuseCMW, vetoWeakJets, allowMPIdipole,
         weakExternal, recoilDeadCone,  doDipoleRecoil, doPartonVertex;
  int    pTmaxMatch, pTdampMatch, alphaSorder, alphaSnfmax, nGluonToQuark,
         weightGluonToQuark, alphaEMorder, nGammaToQuark, nGammaToLepton,
         nCHV, idHV, alphaHVorder, nMaxGlobalRecoil, weakMode;
  double pTdampFudge, mc, mb, m2c, m2b, renormMultFac, factorMultFac,
         fixedFacScale2, alphaSvalue, alphaS2pi, Lambda3flav, Lambda4flav,
         Lambda5flav, Lambda3flav2, Lambda4flav2, Lambda5flav2,
         scaleGluonToQuark, extraGluonToQuark, pTcolCutMin, pTcolCut,
         pT2colCut, pTchgQCut, pT2chgQCut, pTchgLCut, pT2chgLCut,
         pTweakCut, pT2weakCut, mMaxGamma, m2MaxGamma, octetOniumFraction,
         octetOniumColFac, mZ, gammaZ, thetaWRat, mW, gammaW, CFHV, nFlavHV,
         alphaHVfix, LambdaHV, pThvCut, pT2hvCut, mHV, pTmaxFudgeMPI,
         weakEnhancement, vetoWeakDeltaR2;

  // alphaStrong and alphaEM calculations.
  AlphaStrong alphaS;
  AlphaEM     alphaEM;

  // Weak matrix elements used for corrections both of ISR and FSR.
  SimpleWeakShowerMEs  simpleWeakShowerMEs;

  // Some current values.
  bool   twoHard, dopTlimit1, dopTlimit2, dopTdamp;
  double pT2damp, kRad, kEmt, pdfScale2;

  // Bookkeeping of enhanced  actual or trial emissions (see EPJC (2013) 73).
  bool doTrialNow, canEnhanceEmission, canEnhanceTrial, canEnhanceET,
       doUncertaintiesNow;
  string splittingNameNow, splittingNameSel;
  map< double, pair<string,double> > enhanceFactors;
  void storeEnhanceFactor(double pT2, string name, double enhanceFactorIn)
    { enhanceFactors.insert(make_pair(pT2,make_pair(name,enhanceFactorIn)));}

  // All dipole ends and a pointer to the selected hardest dipole end.
  vector<TimeDipoleEnd> dipEnd;
  TimeDipoleEnd* dipSel;
  int iDipSel;

  // Setup a dipole end, either QCD, QED/photon, weak or Hidden Valley one.
  void setupQCDdip( int iSys, int i, int colTag,  int colSign, Event& event,
    bool isOctetOnium = false, bool limitPTmaxIn = true);
  void setupQEDdip( int iSys, int i, int chgType, int gamType, Event& event,
    bool limitPTmaxIn = true);
  void setupWeakdip( int iSys, int i,int weakType, Event& event,
    bool limitPTmaxIn = true);
  // Special setup for weak dipoles if already specified in info ptr.
  void setupWeakdipExternal(Event& event, bool limitPTmaxIn = true);
  void setupHVdip( int iSys, int i, Event& event, bool limitPTmaxIn = true);

  // Evolve a QCD dipole end.
  void pT2nextQCD( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

  // Evolve a QED dipole end, either charged or photon.
  void pT2nextQED( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

  // Evolve a weak dipole end.
  void pT2nextWeak( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

  // Evolve a Hidden Valley dipole end.
  void pT2nextHV( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& );

  // Find kind of QCD ME correction.
  void findMEtype( Event& event, TimeDipoleEnd& dip);

  // Find type of particle; used by findMEtype.
  int findMEparticle( int id, bool isHiddenColour = false);

  // Find mixture of V and A in gamma/Z: energy- and flavour-dependent.
  double gammaZmix( Event& event, int iRes, int iDau1, int iDau2);

  // Set up to calculate QCD ME correction with calcMEcorr.
  double findMEcorr(TimeDipoleEnd* dip, Particle& rad, Particle& partner,
   Particle& emt, bool cutEdge = true);

  // Set up to calculate weak ME corrections for t-channel processes.
  double findMEcorrWeak(TimeDipoleEnd* dip, Vec4 rad, Vec4 rec,
   Vec4 emt, Vec4 p3, Vec4 p4, Vec4 radBef, Vec4 recBef);

  // Calculate value of QCD ME correction.
  double calcMEcorr( int kind, int combiIn, double mixIn, double x1,
    double x2, double r1, double r2, double r3 = 0., bool cutEdge = true);

  // Find coefficient of azimuthal asymmetry from gluon polarization.
  void findAsymPol( Event& event, TimeDipoleEnd* dip);

  // Rescatter: propagate dipole recoil to internal lines connecting systems.
  bool rescatterPropagateRecoil( Event& event, Vec4& pNew);

  // Properties stored for (some) global recoil schemes.
  // Vectors of event indices defining the hard process.
  vector<int> hardPartons;
  // Number of partons in current hard event, number of partons in Born-type
  // hard event (to distinguish between S and H), maximally allowed number of
  // global recoil branchings.
  int nHard, nFinalBorn, nMaxGlobalBranch;
  // Number of proposed splittings in hard scattering systems.
  map<int,int> nProposed;
  // Number of splittings with global recoil (currently only 1).
  int nGlobal, globalRecoilMode;
  // Switch to constrain recoiling system.
  bool limitMUQ;

  // Calculate uncertainty-band weights for accepted/rejected trial branching.
  // Usage: calcUncertainties( accept, pAccept, enhance, vp, dip,
  //   radPtr, emtPtr, recPtr).
  void calcUncertainties(bool , double , double , double ,
    TimeDipoleEnd* , Particle* , Particle* , Particle* );

  // 2 -> 2 information needed for the external weak setup.
  vector<Vec4> weakMomenta;
  vector<int> weak2to2lines;
  int weakHardSize;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SimpleTimeShower_H
