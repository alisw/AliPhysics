// SimpleSpaceShower.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the original simple spacelike initial-state showers.
// SpaceDipoleEnd: data on a radiating dipole end in ISR.
// SimpleSpaceShower: handles the showering description.

#ifndef Pythia8_SimpleSpaceShower_H
#define Pythia8_SimpleSpaceShower_H

#include "Pythia8/SpaceShower.h"
#include "Pythia8/SimpleWeakShowerMEs.h"

namespace Pythia8 {

//==========================================================================

// Data on radiating dipole ends, only used inside SimpleSpaceShower.

class SpaceDipoleEnd {

public:

  // Constructor.
  SpaceDipoleEnd( int systemIn = 0, int sideIn = 0, int iRadiatorIn = 0,
    int iRecoilerIn = 0, double pTmaxIn = 0., int colTypeIn = 0,
    int chgTypeIn = 0, int weakTypeIn = 0,  int MEtypeIn = 0,
    bool normalRecoilIn = true, int weakPolIn = 0,
    int iColPartnerIn = 0, int idColPartnerIn = 0) :
    system(systemIn), side(sideIn), iRadiator(iRadiatorIn),
    iRecoiler(iRecoilerIn), pTmax(pTmaxIn), colType(colTypeIn),
    chgType(chgTypeIn), weakType(weakTypeIn), MEtype(MEtypeIn),
    normalRecoil(normalRecoilIn), weakPol(weakPolIn),
    iColPartner(iColPartnerIn), idColPartner(idColPartnerIn),
    nBranch(0), idDaughter(), idMother(), idSister(), iFinPol(), x1(), x2(),
    m2Dip(), pT2(), z(), xMo(), Q2(), mSister(), m2Sister(), pT2corr(),
    pT2Old(0.), zOld(0.5), asymPol(), m2IF(), mColPartner(),
    pAccept() { }

  // Store values for trial emission.
  void store( int idDaughterIn, int idMotherIn, int idSisterIn,
    double x1In, double x2In, double m2DipIn, double pT2In, double zIn,
    double xMoIn, double Q2In, double mSisterIn, double m2SisterIn,
    double pT2corrIn, int iColPartnerIn, double m2IFIn, double mColPartnerIn)
    {idDaughter = idDaughterIn; idMother = idMotherIn;
    idSister = idSisterIn; x1 = x1In; x2 = x2In; m2Dip = m2DipIn;
    pT2 = pT2In; z = zIn; xMo = xMoIn; Q2 = Q2In; mSister = mSisterIn;
    m2Sister = m2SisterIn; pT2corr = pT2corrIn; iColPartner = iColPartnerIn;
    m2IF = m2IFIn; mColPartner = mColPartnerIn;}

  // Basic properties related to evolution and matrix element corrections.
  int    system, side, iRadiator, iRecoiler;
  double pTmax;
  int    colType, chgType, weakType, MEtype;
  bool   normalRecoil;
  int    weakPol, iColPartner, idColPartner;

  // Properties specific to current trial emission.
  int    nBranch, idDaughter, idMother, idSister, iFinPol;
  double x1, x2, m2Dip, pT2, z, xMo, Q2, mSister, m2Sister, pT2corr,
         pT2Old, zOld, asymPol, m2IF, mColPartner;

  // Properties needed for the evaluation of parameter variations
  double pAccept;

} ;

//==========================================================================

// The SimpleSpaceShower class does spacelike showers.

class SimpleSpaceShower : public SpaceShower {

public:

  // Constructor.
  SimpleSpaceShower() : rescatterFail(), gamma2qqbar(), hasWeaklyRadiated(),
    iSysSel(), pTmaxFudge(), doQCDshower(), doQEDshowerByQ(), doQEDshowerByL(),
    useSamePTasMPI(), doWeakShower(), doMEcorrections(), doMEafterFirst(),
    doPhiPolAsym(), doPhiPolAsymHard(), doPhiIntAsym(), doRapidityOrder(),
    useFixedFacScale(), doSecondHard(), canVetoEmission(), hasUserHooks(),
    alphaSuseCMW(), singleWeakEmission(), vetoWeakJets(), weakExternal(),
    doRapidityOrderMPI(), doMPI(), doDipoleRecoil(), doPartonVertex(),
    pTmaxMatch(), pTdampMatch(), alphaSorder(), alphaSnfmax(), alphaEMorder(),
    nQuarkIn(), enhanceScreening(), weakMode(), pT0paramMode(), pTdampFudge(),
    mc(), mb(), m2c(), m2b(), renormMultFac(), factorMultFac(),
    fixedFacScale2(), alphaSvalue(), alphaS2pi(), Lambda3flav(), Lambda4flav(),
    Lambda5flav(), Lambda3flav2(), Lambda4flav2(), Lambda5flav2(), pT0Ref(),
    ecmRef(), ecmPow(), pTmin(), sCM(), eCM(), pT0(), pTminChgQ(), pTminChgL(),
    pT20(), pT2min(), pT2minChgQ(), pT2minChgL(), pTweakCut(), pT2weakCut(),
    pTmaxFudgeMPI(), strengthIntAsym(), weakEnhancement(), mZ(), gammaZ(),
    thetaWRat(), mW(), gammaW(), weakMaxWt(), vetoWeakDeltaR2(), sideA(),
    twoHard(), dopTlimit1(), dopTlimit2(), dopTdamp(),  tChannel(),
    doUncertaintiesNow(), iNow(), iRec(), idDaughter(), nRad(), idResFirst(),
    idResSecond(), xDaughter(), x1Now(), x2Now(), m2ColPair(), mColPartner(),
    m2ColPartner(), m2Dip(), m2Rec(), pT2damp(), pTbegRef(), pdfScale2(),
    doTrialNow(), canEnhanceEmission(), canEnhanceTrial(), canEnhanceET(),
    iDipNow(), iSysNow(), dipEndNow(), iDipSel(), dipEndSel() {
    beamOffset = 0;}

  // Destructor.
  virtual ~SimpleSpaceShower() {}

  // Initialize generation. Possibility to force re-initialization by hand.
  virtual void init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn);

  // Find whether to limit maximum scale of emissions, and whether to dampen.
  virtual bool limitPTmax( Event& event, double Q2Fac = 0.,
    double Q2Ren = 0.);

  // Prepare system for evolution; identify ME.
  virtual void prepare( int iSys, Event& event, bool limitPTmaxIn = true);

  // Update dipole list after each FSR emission.
  virtual void update( int iSys, Event& event, bool hasWeakRad = false);

  // Select next pT in downwards evolution.
  virtual double pTnext( Event& event, double pTbegAll, double pTendAll,
    int nRadIn = -1, bool doTrialIn = false);

  // ME corrections and kinematics that may give failure.
  virtual bool branch( Event& event);

  // Print dipole list; for debug mainly.
  virtual void list() const;

  // Initialize data members for calculation of uncertainty bands.
  virtual bool initUncertainties();

  // Flag for failure in branch(...) that will force a retry of parton level.
  virtual bool doRestart() const {return rescatterFail;}

  // Tell if latest scattering was a gamma->qqbar.
  virtual bool wasGamma2qqbar() { return gamma2qqbar; }

  // Tell whether ISR has done a weak emission.
  virtual bool getHasWeaklyRadiated() {return hasWeaklyRadiated;}

  // Tell which system was the last processed one.
  virtual int system() const {return iSysSel;}

  // Potential enhancement factor of pTmax scale for hardest emission.
  virtual double enhancePTmax() const {return pTmaxFudge;}

private:

  // Constants: could only be changed in the code itself.
  static const int    MAXLOOPTINYPDF;
  static const double MCMIN, MBMIN, CTHRESHOLD, BTHRESHOLD, EVALPDFSTEP,
         TINYPDF, TINYKERNELPDF, TINYPT2, HEAVYPT2EVOL, HEAVYXEVOL,
         EXTRASPACEQ, LAMBDA3MARGIN, PT2MINWARN, LEPTONXMIN, LEPTONXMAX,
         LEPTONPT2MIN, LEPTONFUDGE, WEAKPSWEIGHT, HEADROOMQ2Q, HEADROOMQ2G,
         HEADROOMG2G, HEADROOMG2Q, HEADROOMHQG, REJECTFACTOR, PROBLIMIT;

  // Store properties to be returned by methods.
  bool   rescatterFail, gamma2qqbar, hasWeaklyRadiated;
  int    iSysSel;
  double pTmaxFudge;

  // Initialization data, normally only set once.
  bool   doQCDshower, doQEDshowerByQ, doQEDshowerByL, useSamePTasMPI,
         doWeakShower, doMEcorrections, doMEafterFirst, doPhiPolAsym,
         doPhiPolAsymHard, doPhiIntAsym, doRapidityOrder, useFixedFacScale,
         doSecondHard, canVetoEmission, hasUserHooks, alphaSuseCMW,
         singleWeakEmission, vetoWeakJets, weakExternal, doRapidityOrderMPI,
         doMPI, doDipoleRecoil, doPartonVertex;
  int    pTmaxMatch, pTdampMatch, alphaSorder, alphaSnfmax, alphaEMorder,
         nQuarkIn, enhanceScreening, weakMode, pT0paramMode;
  double pTdampFudge, mc, mb, m2c, m2b, renormMultFac, factorMultFac,
         fixedFacScale2, alphaSvalue, alphaS2pi, Lambda3flav, Lambda4flav,
         Lambda5flav, Lambda3flav2, Lambda4flav2, Lambda5flav2, pT0Ref,
         ecmRef, ecmPow, pTmin, sCM, eCM, pT0, pTminChgQ, pTminChgL, pT20,
         pT2min, pT2minChgQ, pT2minChgL, pTweakCut, pT2weakCut, pTmaxFudgeMPI,
         strengthIntAsym, weakEnhancement, mZ, gammaZ, thetaWRat, mW, gammaW,
         weakMaxWt, vetoWeakDeltaR2;

  // alphaStrong and alphaEM calculations.
  AlphaStrong alphaS;
  AlphaEM alphaEM;

  // Weak matrix elements used for corrections both of ISR and FSR.
  SimpleWeakShowerMEs  simpleWeakShowerMEs;

  // Some current values.
  bool   sideA, twoHard, dopTlimit1, dopTlimit2, dopTdamp, tChannel,
         doUncertaintiesNow;
  int    iNow, iRec, idDaughter, nRad, idResFirst, idResSecond;
  double xDaughter, x1Now, x2Now, m2ColPair, mColPartner, m2ColPartner,
         m2Dip, m2Rec, pT2damp, pTbegRef, pdfScale2;

  // Bookkeeping of enhanced  actual or trial emissions (see EPJC (2013) 73).
  bool doTrialNow, canEnhanceEmission, canEnhanceTrial, canEnhanceET;
  string splittingNameNow, splittingNameSel;
  map< double, pair<string,double> > enhanceFactors;
  void storeEnhanceFactor(double pT2, string name, double enhanceFactorIn)
    { enhanceFactors.insert(make_pair(pT2,make_pair(name,enhanceFactorIn)));}

  // List of emissions in different sides in different systems:
  vector<int> nRadA,nRadB;

  // All dipole ends
  vector<SpaceDipoleEnd> dipEnd;

  // List of 2 -> 2 momenta for external weak setup.
  vector<Vec4> weakMomenta;

  // Pointers to the current and hardest (so far) dipole ends.
  int iDipNow, iSysNow;
  SpaceDipoleEnd* dipEndNow;
  int iDipSel;
  SpaceDipoleEnd* dipEndSel;

  // Evolve a QCD dipole end.
  void pT2nextQCD( double pT2begDip, double pT2endDip);

  // Evolve a QCD and QED dipole end near heavy quark threshold region.
  void pT2nearThreshold( BeamParticle& beam, double m2Massive,
    double m2Threshold, double xMaxAbs, double zMinAbs,
    double zMaxMassive, int iColPartner);

  // Evolve a QED dipole end.
  void pT2nextQED( double pT2begDip, double pT2endDip);

  // Evolve a Weak dipole end.
  void pT2nextWeak( double pT2begDip, double pT2endDip);

  // Find class of ME correction.
  int findMEtype( int iSys, Event& event, bool weakRadiation = false);

  // Provide maximum of expected ME weight; for preweighting of evolution.
  double calcMEmax( int MEtype, int idMother, int idDaughterIn);

  // Provide actual ME weight for current branching.
  double calcMEcorr(int MEtype, int idMother, int idDaughterIn, double M2,
    double z, double Q2,double m2Sister);

  // Provide actual ME weight for t-channel weak emissions.
  double calcMEcorrWeak(int MEtype, double m2, double z,
    double pT2, Vec4 pMother, Vec4 pB, Vec4 pDaughter,
    Vec4 pB0, Vec4 p1, Vec4 p2, Vec4 pSister);

  // Find coefficient of azimuthal asymmetry from gluon polarization.
  void findAsymPol( Event& event, SpaceDipoleEnd* dip);

  // Find a possible colour partner in the case of dipole recoil.
  int findColPartner(Event& event, int iSideA, int iSideB, int iSystem);

  // Calculate uncertainty-band weights for accepted/rejected trial branching.
  void calcUncertainties(bool accept, double pAcceptIn, double pT20in,
    double enhance, double vp, SpaceDipoleEnd* dip, Particle* motherPtr,
    Particle* sisterPtr);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SimpleSpaceShower_H
