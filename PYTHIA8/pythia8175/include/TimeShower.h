// TimeShower.h is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the timelike final-state showers.
// TimeDipoleEnd: data on a radiating dipole end.
// TimeShower: handles the showering description.

#ifndef Pythia8_TimeShower_H
#define Pythia8_TimeShower_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "Info.h"
#include "ParticleData.h"
#include "PartonSystems.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"
#include "UserHooks.h"
#include "MergingHooks.h"

namespace Pythia8 {

//==========================================================================

// Data on radiating dipole ends; only used inside TimeShower class.

class TimeDipoleEnd {

public:

  // Constructors.
  TimeDipoleEnd() : iRadiator(-1), iRecoiler(-1), pTmax(0.), colType(0), 
    chgType(0), gamType(0), isrType(0), system(0), systemRec(0), MEtype(0), 
    iMEpartner(-1), isOctetOnium(false), isHiddenValley(false), colvType(0),
    MEmix(0.), MEorder(true), MEsplit(true), MEgluinoRec(false),
    isFlexible(false), flavour(0), iAunt(0), 
    mRad(0.), m2Rad(0.), mRec(0.), m2Rec(0.), mDip(0.), m2Dip(0.), 
    m2DipCorr(0.), 
    pT2(0.), m2(0.), z(0.), mFlavour(0.), asymPol(0.), flexFactor(0.)
    { }  
  TimeDipoleEnd(int iRadiatorIn, int iRecoilerIn, double pTmaxIn = 0., 
    int colIn = 0, int chgIn = 0, int gamIn = 0, int isrIn = 0, 
    int systemIn = 0, int MEtypeIn = 0, int iMEpartnerIn = -1, 
    bool isOctetOniumIn = false, bool isHiddenValleyIn = false,
    int colvTypeIn = 0, double MEmixIn = 0., bool MEorderIn = true, 
    bool MEsplitIn = true, bool MEgluinoRecIn = false, 
    bool isFlexibleIn = false) : 
    iRadiator(iRadiatorIn), iRecoiler(iRecoilerIn), pTmax(pTmaxIn), 
    colType(colIn), chgType(chgIn), gamType(gamIn), isrType(isrIn), 
    system(systemIn), systemRec(systemIn) , MEtype(MEtypeIn), 
    iMEpartner(iMEpartnerIn), isOctetOnium(isOctetOniumIn), 
    isHiddenValley(isHiddenValleyIn), colvType(colvTypeIn), MEmix(MEmixIn), 
    MEorder (MEorderIn), MEsplit(MEsplitIn), MEgluinoRec(MEgluinoRecIn),
      isFlexible(isFlexibleIn),
      flavour(0), iAunt(0), 
      mRad(0.), m2Rad(0.), mRec(0.), m2Rec(0.), mDip(0.), m2Dip(0.), 
      m2DipCorr(0.), 
      pT2(0.), m2(0.), z(0.), mFlavour(0.), asymPol(0.), flexFactor(0.)
      { }

  // Basic properties related to dipole and matrix element corrections.
  int    iRadiator, iRecoiler;
  double pTmax;
  int    colType, chgType, gamType, isrType, system, systemRec,
         MEtype, iMEpartner;
  bool   isOctetOnium, isHiddenValley;
  int    colvType;
  double MEmix;
  bool   MEorder, MEsplit, MEgluinoRec;
  bool   isFlexible;

  // Properties specific to current trial emission.
  int    flavour, iAunt;
  double mRad, m2Rad, mRec, m2Rec, mDip, m2Dip, m2DipCorr, 
         pT2, m2, z, mFlavour, asymPol, flexFactor;   
  
} ;

//==========================================================================

// The TimeShower class does timelike showers.

class TimeShower {

public:

  // Constructor.
  TimeShower() {beamOffset = 0;}

  // Destructor.
  virtual ~TimeShower() {}

  // Initialize various pointers. 
  // (Separated from rest of init since not virtual.)
  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn, 
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    CoupSM* coupSMPtrIn, PartonSystems* partonSystemsPtrIn, 
    UserHooks* userHooksPtrIn, MergingHooks* mergingHooksPtrIn = 0) {
    infoPtr = infoPtrIn; settingsPtr = settingsPtrIn; 
    particleDataPtr = particleDataPtrIn; rndmPtr = rndmPtrIn; 
    coupSMPtr = coupSMPtrIn; partonSystemsPtr = partonSystemsPtrIn; 
    userHooksPtr = userHooksPtrIn; mergingHooksPtr = mergingHooksPtrIn;}

  // Initialize alphaStrong and related pTmin parameters.
  virtual void init( BeamParticle* beamAPtrIn = 0, 
    BeamParticle* beamBPtrIn = 0);

  // New beams possible for handling of hard diffraction. (Not virtual.)
  void reassignBeamPtrs( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    int beamOffsetIn = 0) {beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn;
    beamOffset = beamOffsetIn;}

  // Find whether to limit maximum scale of emissions, and whether to dampen.
  virtual bool limitPTmax( Event& event, double Q2Fac = 0., 
    double Q2Ren = 0.);

  // Potential enhancement factor of pTmax scale for hardest emission.
  double enhancePTmax() {return pTmaxFudge;}

  // Top-level routine to do a full time-like shower in resonance decay.
  virtual int shower( int iBeg, int iEnd, Event& event, double pTmax,
    int nBranchMax = 0);

  // Top-level routine for QED radiation in hadronic decay to two leptons.
  virtual int showerQED( int i1, int i2, Event& event, double pTmax);

  // Provide the pT scale of the last branching in the above shower.
  double pTLastInShower() {return pTLastBranch;}

  // Prepare system for evolution after each new interaction; identify ME.
  virtual void prepare( int iSys, Event& event, bool limitPTmaxIn = true);

  // Update dipole list after a multiparton interactions rescattering.
  virtual void rescatterUpdate( int iSys, Event& event);

  // Update dipole list after each ISR emission.  
  virtual void update( int iSys, Event& event);

  // Select next pT in downwards evolution.
  virtual double pTnext( Event& event, double pTbegAll, double pTendAll);

  // ME corrections and kinematics that may give failure.
  virtual bool branch( Event& event, bool isInterleaved = false); 

  // Tell which system was the last processed one.
  int system() const {return iSysSel;}; 

  // Print dipole list; for debug mainly.
  virtual void list( ostream& os = cout) const;

protected:

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the settings database.
  Settings*      settingsPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointer to Standard Model couplings.
  CoupSM*        coupSMPtr;

  // Pointers to the two incoming beams. Offset their location in event.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  int            beamOffset;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Pointer to userHooks object for user interaction with program.
  UserHooks*     userHooksPtr;

  // Store properties to be returned by methods.
  int    iSysSel;
  double pTmaxFudge, pTLastBranch;

private:

  // Constants: could only be changed in the code itself.
  static const double SIMPLIFYROOT, XMARGIN, XMARGINCOMB, TINYPDF, LARGEM2, 
                      THRESHM2, LAMBDA3MARGIN;
  // Rescatter: try to fix up recoil between systems
  static const bool   FIXRESCATTER, VETONEGENERGY;
  static const double MAXVIRTUALITYFRACTION, MAXNEGENERGYFRACTION;

  // Initialization data, normally only set once.
  bool   doQCDshower, doQEDshowerByQ, doQEDshowerByL, doQEDshowerByGamma, 
         doMEcorrections, doMEafterFirst, doPhiPolAsym, doInterleave, 
         allowBeamRecoil, dampenBeamRecoil, recoilToColoured, 
         allowRescatter, canVetoEmission, doHVshower, brokenHVsym,
         globalRecoil, useLocalRecoilNow, doSecondHard;
  int    pTmaxMatch, pTdampMatch, alphaSorder, nGluonToQuark, 
         alphaEMorder, nGammaToQuark, nGammaToLepton, nCHV, idHV,
         nMaxGlobalRecoil;
  double pTdampFudge, mc, mb, m2c, m2b, renormMultFac, factorMultFac,
         alphaSvalue, alphaS2pi, Lambda3flav, Lambda4flav, Lambda5flav, 
         Lambda3flav2, Lambda4flav2, Lambda5flav2, pTcolCutMin, pTcolCut, 
         pT2colCut, pTchgQCut, pT2chgQCut, pTchgLCut, pT2chgLCut, 
         mMaxGamma, m2MaxGamma, octetOniumFraction, octetOniumColFac, 
         mZ, gammaZ, thetaWRat, CFHV, alphaHVfix, pThvCut, pT2hvCut, 
         mHV, pTmaxFudgeMPI;

  // alphaStrong and alphaEM calculations.
  AlphaStrong alphaS;
  AlphaEM     alphaEM;

  // Some current values.
  bool   dopTlimit1, dopTlimit2, dopTdamp;
  double pT2damp, kRad, kEmt;

  // All dipole ends and a pointer to the selected hardest dipole end.
  vector<TimeDipoleEnd> dipEnd;
  TimeDipoleEnd* dipSel;
  int iDipSel;

  // Setup a dipole end, either QCD, QED/photon or Hidden Valley one.
  void setupQCDdip( int iSys, int i, int colTag,  int colSign, Event& event,
    bool isOctetOnium = false, bool limitPTmaxIn = true);
  void setupQEDdip( int iSys, int i, int chgType, int gamType, Event& event, 
    bool limitPTmaxIn = true); 
  void setupHVdip( int iSys, int i, Event& event, bool limitPTmaxIn = true); 

  // Evolve a QCD dipole end. 
  void pT2nextQCD( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

  // Evolve a QED dipole end, either charged or photon. 
  void pT2nextQED( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
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

  // Calculate value of QCD ME correction.
  double calcMEcorr( int kind, int combiIn, double mixIn, double x1, 
    double x2, double r1, double r2, double r3 = 0., bool cutEdge = true);

  // Find coefficient of azimuthal asymmetry from gluon polarization.
  void findAsymPol( Event& event, TimeDipoleEnd* dip);

  // Rescatter: propagate dipole recoil to internal lines connecting systems.
  bool rescatterPropagateRecoil( Event& event, Vec4& pNew);

  // Pointer to MergingHooks object for NLO merging.
  MergingHooks* mergingHooksPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_TimeShower_H

