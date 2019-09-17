// PhaseSpace.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for phase space generators in kinematics selection.
// PhaseSpace: base class for phase space generators.
// Base class for derived classes> 2 ->1 , 2 -> 2, 2 -> 2 elastic/diffractive,
// 2 -> 2 nondiffractive, 2 -> 3, Les Houches.

#ifndef Pythia8_PhaseSpace_H
#define Pythia8_PhaseSpace_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Info.h"
#include "Pythia8/LesHouches.h"
#include "Pythia8/MultipartonInteractions.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/SigmaProcess.h"
#include "Pythia8/SigmaTotal.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/UserHooks.h"
#include "Pythia8/GammaKinematics.h"

namespace Pythia8 {

//==========================================================================

// Forward reference to the UserHooks class.
class UserHooks;

//==========================================================================

// PhaseSpace is a base class for  phase space generators
// used in the selection of hard-process kinematics.

class PhaseSpace {

public:

  // Destructor.
  virtual ~PhaseSpace() {}

  // Perform simple initialization and store pointers.
  void init(bool isFirst, SigmaProcess* sigmaProcessPtrIn,
    Info* infoPtrIn, Settings* settingsPtrIn, ParticleData* particleDataPtrIn,
    Rndm* rndmPtrIn, BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    Couplings* couplingsPtrIn, SigmaTotal* sigmaTotPtrIn,
    UserHooks* userHooksPtrIn);

  // Update the CM energy of the event.
  void newECM(double eCMin) {eCM = eCMin; s = eCM * eCM;}

  // Store or replace Les Houches pointer.
  void setLHAPtr(LHAup* lhaUpPtrIn) {lhaUpPtr = lhaUpPtrIn;}

  // A pure virtual method, wherein an optimization procedure
  // is used to determine how phase space should be sampled.
  virtual bool setupSampling() = 0;

  // A pure virtual method, wherein a trial event kinematics
  // is to be selected in the derived class.
  virtual bool trialKin(bool inEvent = true, bool repeatSame = false) = 0;

  // A pure virtual method, wherein the accepted event kinematics
  // is to be constructed in the derived class.
  virtual bool finalKin() = 0;

  // Allow for nonisotropic decays when ME's available.
  void   decayKinematics( Event& process);

  // Give back current or maximum cross section, or set latter.
  double sigmaNow() const {return sigmaNw;}
  double sigmaMax() const {return sigmaMx;}
  double biasSelectionWeight()  const {return biasWt;}
  bool   newSigmaMax() const {return newSigmaMx;}
  void   setSigmaMax(double sigmaMaxIn) {sigmaMx = sigmaMaxIn;}

  // For Les Houches with negative event weight needs
  virtual double sigmaSumSigned() const {return sigmaMx;}

  // Give back constructed four-vectors and known masses.
  Vec4   p(int i)   const {return pH[i];}
  double m(int i)   const {return mH[i];}

  // Reset the four-momentum.
  void   setP(int i, Vec4 pNew) { pH[i] = pNew; }

  // Give back other event properties.
  double ecm()      const {return eCM;}
  double x1()       const {return x1H;}
  double x2()       const {return x2H;}
  double sHat()     const {return sH;}
  double tHat()     const {return tH;}
  double uHat()     const {return uH;}
  double pTHat()    const {return pTH;}
  double thetaHat() const {return theta;}
  double phiHat()   const {return phi;}
  double runBW3()   const {return runBW3H;}
  double runBW4()   const {return runBW4H;}
  double runBW5()   const {return runBW5H;}

  // Inform whether beam particles are resolved in partons or scatter directly.
  virtual bool isResolved() const {return true;}

  // Functions to rescale momenta and cross section for new sHat
  // Currently implemented only for PhaseSpace2to2tauyz class.
  virtual void rescaleSigma( double){}
  virtual void rescaleMomenta( double){}

  // Calculate the weight for over-estimated cross section.
  virtual double weightGammaPDFApprox(){ return 1.;}

  // Set the GammaKinematics pointer needed for soft photoproduction.
  virtual void setGammaKinPtr( GammaKinematics* gammaKinPtrIn) {
    gammaKinPtr = gammaKinPtrIn; }

protected:

  // Constructor.
  PhaseSpace() : sigmaProcessPtr(), infoPtr(), settingsPtr(),
    particleDataPtr(), rndmPtr(), beamAPtr(), beamBPtr(), couplingsPtr(),
    sigmaTotPtr(), userHooksPtr(), lhaUpPtr(), gammaKinPtr(),
    useBreitWigners(), doEnergySpread(), showSearch(), showViolation(),
    increaseMaximum(), hasQ2Min(), gmZmodeGlobal(), mHatGlobalMin(),
    mHatGlobalMax(), pTHatGlobalMin(), pTHatGlobalMax(), Q2GlobalMin(),
    pTHatMinDiverge(), minWidthBreitWigners(), minWidthNarrowBW(), idA(),
    idB(), idAgm(), idBgm(), mA(), mB(), eCM(), s(), sigmaMxGm(),
    hasLeptonBeamA(), hasLeptonBeamB(), hasOneLeptonBeam(),
    hasTwoLeptonBeams(), hasPointGammaA(), hasPointGammaB(),
    hasOnePointParticle(), hasTwoPointParticles(), hasGamma(),
    hasVMD(), newSigmaMx(),
    canModifySigma(), canBiasSelection(), canBias2Sel(), gmZmode(),
    bias2SelPow(), bias2SelRef(), wtBW(), sigmaNw(),
    sigmaMx(), sigmaPos(), sigmaNeg(), biasWt(), mHatMin(), mHatMax(),
    sHatMin(), sHatMax(), pTHatMin(), pTHatMax(), pT2HatMin(), pT2HatMax(),
    x1H(), x2H(), m3(), m4(), m5(), s3(), s4(), s5(), mHat(), sH(), tH(), uH(),
    pAbs(), p2Abs(), pTH(), theta(), phi(), betaZ(), mH(), idResA(), idResB(),
    mResA(), mResB(), GammaResA(), GammaResB(), tauResA(), tauResB(),
    widResA(), widResB(), sameResMass(), useMirrorWeight(), hasNegZ(),
    hasPosZ(), tau(), y(), z(), tauMin(), tauMax(), yMax(), zMin(), zMax(),
    ratio34(), unity34(), zNeg(), zPos(), wtTau(), wtY(), wtZ(), wt3Body(),
    runBW3H(), runBW4H(), runBW5H(), intTau0(), intTau1(), intTau2(),
    intTau3(), intTau4(), intTau5(), intTau6(), intY0(), intY12(), intY34(),
    intY56(), mTchan1(), sTchan1(), mTchan2(), sTchan2(), frac3Flat(),
    frac3Pow1(), frac3Pow2(), zNegMin(), zNegMax(), zPosMin(), zPosMax(),
    nTau(), nY(), nZ(), tauCoef(), yCoef(), zCoef(), tauCoefSum(), yCoefSum(),
    zCoefSum(), useBW(), useNarrowBW(), idMass(), mPeak(), sPeak(), mWidth(),
    mMin(), mMax(), mw(), wmRat(), mLower(), mUpper(), sLower(), sUpper(),
    fracFlatS(), fracFlatM(), fracInv(), fracInv2(), atanLower(), atanUpper(),
    intBW(), intFlatS(), intFlatM(), intInv(), intInv2() {}

  // Constants: could only be changed in the code itself.
  static const int    NMAXTRY, NTRY3BODY;
  static const double SAFETYMARGIN, TINY, EVENFRAC, SAMESIGMA, MRESMINABS,
                      WIDTHMARGIN, SAMEMASS, MASSMARGIN, EXTRABWWTMAX,
                      THRESHOLDSIZE, THRESHOLDSTEP, YRANGEMARGIN, LEPTONXMIN,
                      LEPTONXMAX, LEPTONXLOGMIN, LEPTONXLOGMAX, LEPTONTAUMIN,
                      SHATMINZ, PT2RATMINZ, WTCORRECTION[11];

  // Pointer to cross section.
  SigmaProcess* sigmaProcessPtr;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the settings database.
  Settings*     settingsPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Pointers to incoming beams.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // Pointer to Standard Model couplings.
  Couplings*    couplingsPtr;

  // Pointer to the total/elastic/diffractive cross section object.
  SigmaTotal*   sigmaTotPtr;

  // Pointer to userHooks object for user interaction with program.
  UserHooks*    userHooksPtr;

  // Pointer to LHAup for generating external events.
  LHAup*        lhaUpPtr;

  // Pointer to object that samples photon kinematics from leptons.
  GammaKinematics* gammaKinPtr;

  // Initialization data, normally only set once.
  bool   useBreitWigners, doEnergySpread, showSearch, showViolation,
         increaseMaximum, hasQ2Min;
  int    gmZmodeGlobal;
  double mHatGlobalMin, mHatGlobalMax, pTHatGlobalMin, pTHatGlobalMax,
         Q2GlobalMin, pTHatMinDiverge, minWidthBreitWigners, minWidthNarrowBW;

  // Information on incoming beams.
  int    idA, idB, idAgm, idBgm;
  double mA, mB, eCM, s, sigmaMxGm;
  bool   hasLeptonBeamA, hasLeptonBeamB, hasOneLeptonBeam, hasTwoLeptonBeams,
         hasPointGammaA, hasPointGammaB, hasOnePointParticle,
         hasTwoPointParticles, hasGamma, hasVMD;

  // Cross section information.
  bool   newSigmaMx, canModifySigma, canBiasSelection, canBias2Sel;
  int    gmZmode;
  double bias2SelPow, bias2SelRef, wtBW, sigmaNw, sigmaMx, sigmaPos,
         sigmaNeg, biasWt;

  // Process-specific kinematics properties, almost always available.
  double mHatMin, mHatMax, sHatMin, sHatMax, pTHatMin, pTHatMax,
         pT2HatMin, pT2HatMax;

  // Event-specific kinematics properties, almost always available.
  double x1H, x2H, m3, m4, m5, s3, s4, s5, mHat, sH, tH, uH, pAbs, p2Abs,
         pTH, theta, phi, betaZ;
  Vec4   pH[12];
  double mH[12];

  // Reselect decay products momenta isotropically in phase space.
  void decayKinematicsStep( Event& process, int iRes);

  // Much common code for normal 2 -> 1, 2 -> 2 and 2 -> 3 cases:

  // Determine how phase space should be sampled.
  void setup3Body();
  bool setupSampling123(bool is2, bool is3);

  // Select a trial kinematics phase space point.
  bool trialKin123(bool is2, bool is3, bool inEvent = true);

  // Presence and properties of any s-channel resonances.
  int    idResA, idResB;
  double mResA, mResB, GammaResA, GammaResB, tauResA, tauResB, widResA,
         widResB;
  bool   sameResMass;

  // Kinematics properties specific to 2 -> 1/2/3.
  bool   useMirrorWeight, hasNegZ, hasPosZ;
  double tau, y, z, tauMin, tauMax, yMax, zMin, zMax, ratio34, unity34,
         zNeg, zPos, wtTau, wtY, wtZ, wt3Body, runBW3H, runBW4H, runBW5H,
         intTau0, intTau1, intTau2, intTau3, intTau4, intTau5, intTau6,
         intY0, intY12, intY34, intY56, mTchan1, sTchan1, mTchan2, sTchan2,
         frac3Flat, frac3Pow1, frac3Pow2, zNegMin, zNegMax, zPosMin, zPosMax;
  Vec4   p3cm, p4cm, p5cm;

  // Coefficients for optimized selection in 2 -> 1/2/3.
  int    nTau, nY, nZ;
  double tauCoef[8], yCoef[8], zCoef[8], tauCoefSum[8], yCoefSum[8],
         zCoefSum[8];

  // Calculate kinematical limits for 2 -> 1/2/3.
  bool limitTau(bool is2, bool is3);
  bool limitY();
  bool limitZ();

  // Select kinematical variable between defined limits for 2 -> 1/2/3.
  void selectTau(int iTau, double tauVal, bool is2);
  void selectY(int iY, double yVal);
  void selectZ(int iZ, double zVal);
  bool select3Body();

  // Solve equation system for better phase space coefficients in 2 -> 1/2/3.
  void solveSys( int n, int bin[8], double vec[8], double mat[8][8],
    double coef[8]);

  // Properties specific to resonance mass selection in 2 -> 2 and 2 -> 3.
  bool   useBW[6], useNarrowBW[6];
  int    idMass[6];
  double mPeak[6], sPeak[6], mWidth[6], mMin[6], mMax[6], mw[6], wmRat[6],
         mLower[6], mUpper[6], sLower[6], sUpper[6], fracFlatS[6],
         fracFlatM[6], fracInv[6], fracInv2[6], atanLower[6], atanUpper[6],
         intBW[6], intFlatS[6], intFlatM[6], intInv[6], intInv2[6];

  // Setup mass selection for one resonance at a time. Split in two parts.
  void   setupMass1(int iM);
  void   setupMass2(int iM, double distToThresh);

  // Do mass selection and find the associated weight.
  void   trialMass(int iM);
  double weightMass(int iM);

  // Standard methods to find t range of a 2 -> 2 process
  // and to check whether a given t value is in that range.
  pair<double,double> tRange( double sIn, double s1In, double s2In,
    double s3In,  double s4In) {
    double lambda12 = pow2( sIn - s1In - s2In) - 4. * s1In * s2In;
    double lambda34 = pow2( sIn - s3In - s4In) - 4. * s3In * s4In;
    if (lambda12 < 0. || lambda34 < 0.) return make_pair( 0., 0.);
    double tLow = -0.5 * (sIn - (s1In + s2In + s3In + s4In) + (s1In - s2In)
      * (s3In - s4In) / sIn + sqrtpos(lambda12 *  lambda34) / sIn);
    double tUpp = ( (s3In - s1In) * (s4In - s2In) + (s1In + s4In - s2In - s3In)
      * (s1In * s4In - s2In * s3In) / sIn ) / tLow;
    return make_pair( tLow, tUpp); }
  bool tInRange( double tIn, double sIn, double s1In, double s2In,
    double s3In,  double s4In) {
    pair<double, double> tRng = tRange( sIn, s1In, s2In, s3In, s4In);
    return (tIn > tRng.first && tIn < tRng.second); }

  // The error function erf(x) should normally be in your math library,
  // but if not uncomment this simple parametrization by Sergei Winitzki.
  //double erf(double x) { double x2 = x * x; double kx2 = 0.147 * x2;
  //  double tmp = sqrt(1. - exp(-x2 * (4./M_PI + kx2) / (1. + kx2)));
  //  return ((x >= 0.) ? tmp : -tmp); }

};

//==========================================================================

// A derived class with 2 -> 1 kinematics set up in tau, y.

class PhaseSpace2to1tauy : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to1tauy() {}

  // Optimize subsequent kinematics selection.
  virtual bool setupSampling() {if (!setupMass()) return false;
    return setupSampling123(false, false);}

  // Construct the trial kinematics.
  virtual bool trialKin(bool inEvent = true, bool = false) {wtBW = 1.;
    return trialKin123(false, false, inEvent);}

  // Construct the final event kinematics.
  virtual bool finalKin();

private:

  // Set up allowed mass range.
  bool setupMass();

};

//==========================================================================

// A derived class with 2 -> 2 kinematics set up in tau, y, z = cos(theta).

class PhaseSpace2to2tauyz : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to2tauyz() {}

  // Optimize subsequent kinematics selection.
  virtual bool setupSampling() {if (!setupMasses()) return false;
    return setupSampling123(true, false);}

  // Construct the trial kinematics.
  virtual bool trialKin(bool inEvent = true, bool = false) {
    if (!trialMasses()) return false;
    return trialKin123(true, false, inEvent);}

  // Construct the final event kinematics.
  virtual bool finalKin();

  // Rescales the momenta of incoming and outgoing partons according to
  // new sHat.
  virtual void rescaleMomenta ( double sHatNew);

  // Recalculates cross section with rescaled sHat.
  virtual void rescaleSigma ( double sHatNew);

  // Calculate the weight for over-estimated cross section.
  virtual double weightGammaPDFApprox();

private:

  // Set up for fixed or Breit-Wigner mass selection.
  bool setupMasses();

  // Select fixed or Breit-Wigner-distributed masses.
  bool trialMasses();

  // Pick off-shell initialization masses when on-shell not allowed.
  bool constrainedM3M4();
  bool constrainedM3();
  bool constrainedM4();

};

//==========================================================================

// A derived class with 2 -> 2 kinematics set up for elastic scattering.

class PhaseSpace2to2elastic : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to2elastic() : isOneExp(), useCoulomb(), s1(), s2(), alphaEM0(),
    lambda12S(), lambda12(), lambda34(), tempA(), tempB(), tempC(),
    tLow(), tUpp(), bSlope1(), bSlope2(), sigRef1(), sigRef2(),
    sigRef(), sigNorm1(), sigNorm2(), sigNorm3(), sigNormSum(), rel2() {}

  // Construct the trial or final event kinematics.
  virtual bool setupSampling();
  virtual bool trialKin(bool inEvent = true, bool = false);
  virtual bool finalKin();

  // Are beam particles resolved in partons or scatter directly?
  virtual bool isResolved() const {return false;}

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRY;
  static const double HBARC2, BNARROW, BWIDE, WIDEFRAC, TOFFSET;

  // Kinematics properties specific to 2 -> 2 elastic.
  bool   isOneExp, useCoulomb;
  double s1, s2, alphaEM0, lambda12S, lambda12, lambda34, tempA, tempB, tempC,
         tLow, tUpp, bSlope1, bSlope2, sigRef1, sigRef2, sigRef,
         sigNorm1, sigNorm2, sigNorm3, sigNormSum, rel2;

};

//==========================================================================

// A derived class with 2 -> 2 kinematics set up for diffractive scattering.

class PhaseSpace2to2diffractive : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to2diffractive(bool isDiffAin = false, bool isDiffBin = false)
    : isDiffA(isDiffAin), isDiffB(isDiffBin), splitxit(), m3ElDiff(),
    m4ElDiff(), s1(), s2(), xiMin(), xiMax(), xiNow(), sigNow(), sigMax(),
    sigMaxNow(), lambda12(), lambda34(), bNow(), tempA(), tempB(), tempC(),
    tLow(), tUpp(), tWeight(), fWid1(), fWid2(), fWid3(), fWid4(), fbWid1(),
    fbWid2(), fbWid3(), fbWid4(), fbWid1234() {isSD = !isDiffA || !isDiffB;}

  // Construct the trial or final event kinematics.
  virtual bool setupSampling();
  virtual bool trialKin(bool inEvent = true, bool = false);
  virtual bool finalKin();

  // Are beam particles resolved in partons or scatter directly?
  virtual bool isResolved() const {return false;}

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRY;
  static const double BWID1, BWID2, BWID3, BWID4, FWID1SD, FWID2SD, FWID3SD,
                      FWID4SD, FWID1DD, FWID2DD, FWID3DD, FWID4DD, MAXFUDGESD,
                      MAXFUDGEDD, MAXFUDGET, DIFFMASSMARGIN, SPROTON;

  // Initialization data.
  bool   isDiffA, isDiffB, isSD, splitxit;

  // Kinematics properties specific to 2 -> 2 diffraction.
  double m3ElDiff, m4ElDiff, s1, s2, xiMin, xiMax, xiNow, sigNow, sigMax,
         sigMaxNow, lambda12, lambda34, bNow, tempA, tempB, tempC,
         tLow, tUpp, tWeight, fWid1, fWid2, fWid3, fWid4, fbWid1, fbWid2,
         fbWid3, fbWid4, fbWid1234;

};

//==========================================================================

// A derived class with 2 -> 3 kinematics set up for central diffractive
// scattering.

class PhaseSpace2to3diffractive : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to3diffractive() : splitxit(), s1(), s2(), m5min(), s5min(), m5(),
    sigNow(), sigMax(), sigMaxNow(), xiMin(), xi1(), xi2(), fWid1(), fWid2(),
    fWid3(), fbWid1(), fbWid2(), fbWid3(), fbWid123() {}

  // Construct the trial or final event kinematics.
  virtual bool setupSampling();
  virtual bool trialKin(bool inEvent = true, bool = false);
  virtual bool finalKin();

  // Are beam particles resolved in partons or scatter directly?
  virtual bool isResolved() const {return false;}

 private:

  // Constants: could only be changed in the code itself.
  static const int    NTRY;
  static const double BWID1, BWID2, BWID3, BWID4, FWID1, FWID2, FWID3,
                      MAXFUDGECD, MAXFUDGET, DIFFMASSMARGIN;

  // Initialization data.
  bool   splitxit;

  // Local variables to calculate DPE kinematics.
  double s1, s2, m5min, s5min, m5, sigNow, sigMax, sigMaxNow, xiMin, xi1, xi2,
         fWid1, fWid2, fWid3, fbWid1, fbWid2, fbWid3, fbWid123;
  Vec4   p1, p2, p3, p4, p5;

};

//==========================================================================

class PhaseSpace2to2nondiffractive : public PhaseSpace {

// A derived class for nondiffractive events. Hardly does anything, since
// the real action is taken care of by the MultipartonInteractions class.

public:

  // Constructor.
  PhaseSpace2to2nondiffractive() {}

  // Construct the trial or final event kinematics.
  virtual bool setupSampling();
  virtual bool trialKin( bool , bool = false);
  virtual bool finalKin() { if (hasGamma) gammaKinPtr->finalize();
    return true;}

private:

};

//==========================================================================

// A derived class with 2 -> 3 kinematics 1 + 2 -> 3 + 4 + 5 set up in
// tau, y, pT2_4, pT2_5, phi_4, phi_5 and y_3 (partial cylindrical symmetry).

class PhaseSpace2to3tauycyl : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to3tauycyl() {}

  // Optimize subsequent kinematics selection.
  virtual bool setupSampling() {if (!setupMasses()) return false;
    setup3Body(); return setupSampling123(false, true);}

  // Construct the trial kinematics.
  virtual bool trialKin(bool inEvent = true, bool = false) {
    if (!trialMasses()) return false;
    return trialKin123(false, true, inEvent);}

  // Construct the final event kinematics.
  virtual bool finalKin();

private:

  // Constants: could only be changed in the code itself.
  static const int    NITERNR;

  // Set up for fixed or Breit-Wigner mass selection.
  bool setupMasses();

  // Select fixed or Breit-Wigner-distributed masses.
  bool trialMasses();

};

//==========================================================================

// A derived class with 2 -> 3 kinematics 1 + 2 -> 3 + 4 + 5 set up in
// y3, y4, y5, pT2_3, pT2_5, phi_3 and phi_5, and with R separation cut.
// Intended specifically for (essentially massless) 2 -> 3 QCD processes.

class PhaseSpace2to3yyycyl : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to3yyycyl() : pTHat3Min(), pTHat3Max(), pTHat5Min(), pTHat5Max(),
    RsepMin(), R2sepMin(), hasBaryonBeams(), pT3Min(), pT3Max(), pT5Min(),
    pT5Max(), y3Max(), y4Max(), y5Max(), pT3(), pT4(), pT5(), phi3(), phi4(),
    phi5(), y3(), y4(), y5(), dphi() {}

  // Optimize subsequent kinematics selection.
  virtual bool setupSampling();

  // Construct the trial kinematics.
  virtual bool trialKin(bool inEvent = true, bool = false);

  // Construct the final event kinematics.
  virtual bool finalKin();

private:

  // Phase space cuts specifically for 2 -> 3 QCD processes.
  double pTHat3Min, pTHat3Max, pTHat5Min, pTHat5Max, RsepMin, R2sepMin;
  bool   hasBaryonBeams;

  // Event kinematics choices.
  double pT3Min, pT3Max, pT5Min, pT5Max, y3Max, y4Max, y5Max,
         pT3, pT4, pT5, phi3, phi4, phi5, y3, y4, y5, dphi;
  Vec4   pInSum;

};

//==========================================================================

// A derived class for Les Houches events.

class PhaseSpaceLHA : public PhaseSpace {

public:

  // Constructor.
  PhaseSpaceLHA() : strategy(), stratAbs(), nProc(), idProcSave(0),
    xMaxAbsSum(), xSecSgnSum(), sigmaSgn() {}

  // Find maximal cross section for comparison with internal processes.
  virtual bool setupSampling();

  // Construct the next process, by interface to Les Houches class.
  virtual bool trialKin( bool , bool repeatSame = false);

  // Set scale, alpha_s and alpha_em if not done.
  virtual bool finalKin() {sigmaProcessPtr->setScale(); return true;}

  // For Les Houches with negative event weight needs
  virtual double sigmaSumSigned() const {return sigmaSgn;}

private:

  // Constants.
  static const double CONVERTPB2MB;

  // Local properties.
  int    strategy, stratAbs, nProc, idProcSave;
  double xMaxAbsSum, xSecSgnSum, sigmaSgn;
  vector<int>    idProc;
  vector<double> xMaxAbsProc;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PhaseSpace_H
