// ResonanceWidths.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for resonance properties: dynamical widths etc.
// ResonanceWidths: base class for all resonances.
// ResonanceGmZ, ...: derived classes for individual resonances.

#ifndef Pythia8_ResonanceWidths_H
#define Pythia8_ResonanceWidths_H

#include "Pythia8/Basics.h"
#include "Pythia8/Info.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"

namespace Pythia8 {

//==========================================================================

// Forward references to ParticleData and StandardModel classes.
class DecayChannel;
class ParticleData;
class ParticleDataEntry;
class Couplings;

//==========================================================================

// The ResonanceWidths is the base class. Also used for generic resonaces.

class ResonanceWidths {

public:

  // Destructor.
  virtual ~ResonanceWidths() {}

  // Set up standard properties.
  void initBasic(int idResIn, bool isGenericIn = false) {
    idRes = idResIn; isGeneric = isGenericIn;}

  // Calculate and store partial and total widths at the nominal mass.
  virtual bool init(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, Couplings* couplingsPtrIn);

  // Return identity of particle species.
  int id() const {return idRes;}

  // Calculate the total/open width for given mass, charge and instate.
  double width(int idSgn, double mHatIn, int idInFlavIn = 0,
    bool openOnly = false, bool setBR = false, int idOutFlav1 = 0,
    int idOutFlav2 = 0);

  // Special case to calculate open final-state width.
  double widthOpen(int idSgn, double mHatIn, int idIn = 0) {
    return width( idSgn, mHatIn, idIn, true, false);}

  // Special case to store open final-state widths for channel selection.
  double widthStore(int idSgn, double mHatIn, int idIn = 0) {
    return width( idSgn, mHatIn, idIn, true, true);}

  // Return fraction of width open for particle and antiparticle.
  double openFrac(int idSgn) {return (idSgn > 0) ? openPos : openNeg;}

  // Return forced rescaling factor of resonance width.
  double widthRescaleFactor() {return forceFactor;}

  // Special case to calculate one final-state width.
  // Currently only used for Higgs -> qqbar, g g or gamma gamma.
  double widthChan(double mHatIn, int idOutFlav1, int idOutFlav2) {
    return width( 1, mHatIn, 0, false, false, idOutFlav1, idOutFlav2);}

protected:

  // Constructor.
  ResonanceWidths() : idRes(), hasAntiRes(), doForceWidth(), isGeneric(),
    allowCalcWidth(), minWidth(), minThreshold(), mRes(), GammaRes(), m2Res(),
    GamMRat(), openPos(), openNeg(), forceFactor(), iChannel(), onMode(),
    meMode(), mult(), id1(), id2(), id3(), id1Abs(), id2Abs(), id3Abs(),
    idInFlav(), widNow(), mHat(), mf1(), mf2(), mf3(), mr1(), mr2(), mr3(),
    ps(), kinFac(), alpEM(), alpS(), colQ(), preFac(), particlePtr(),
    infoPtr(), settingsPtr(), particleDataPtr(), couplingsPtr() {}

  // Constants: could only be changed in the code itself.
  static const int    NPOINT;
  static const double MASSMIN, MASSMARGIN;

  // Particle properties always present.
  int    idRes, hasAntiRes;
  bool   doForceWidth, isGeneric, allowCalcWidth;
  double minWidth, minThreshold, mRes, GammaRes, m2Res, GamMRat,
         openPos, openNeg, forceFactor;

  // Properties for currently studied decay channel(s).
  int    iChannel, onMode, meMode, mult, id1, id2, id3, id1Abs,
         id2Abs, id3Abs, idInFlav;
  double widNow, mHat, mf1, mf2, mf3, mr1, mr2, mr3, ps, kinFac,
         alpEM, alpS, colQ, preFac;

  // Pointer to properties of the particle species.
  ParticleDataEntry* particlePtr;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the settings database.
  Settings*     settingsPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointers to Standard Model and SUSY couplings.
  Couplings*    couplingsPtr;

  // Initialize constants.
  virtual void initConstants() {}

  // Virtual methods to handle model-specific (non-SM) part of initialization
  // for use by derived classes that implement additional models (eg SUSY).
  virtual bool initBSM() {return true;}
  virtual bool allowCalc() {return true;}

  // Calculate various common prefactors for the current mass.
  // Optional argument calledFromInit only used for Z0.
  virtual void calcPreFac(bool = false) {}

  // Calculate width for currently considered channel.
  // Optional argument calledFromInit only used for Z0.
  virtual void calcWidth(bool = false) {}

  // Simple routines for matrix-element integration over Breit-Wigners.
  double numInt1BW(double mHatIn, double m1, double Gamma1, double mMin1,
    double m2, int psMode = 1);
  double numInt2BW(double mHatIn, double m1, double Gamma1, double mMin1,
    double m2, double Gamma2, double mMin2, int psMode = 1);

};

//==========================================================================

// The ResonanceGeneric class handles a generic resonance.
// Only needs a constructor and allowCalc = false; for the rest uses
// defaults in base class.

class ResonanceGeneric : public ResonanceWidths {

public:

  // Constructor.
  ResonanceGeneric(int idResIn) {initBasic(idResIn, true);}

  // By default, assume no dedicated code exists to compute width.
  virtual bool allowCalc() {return false;}

};

//==========================================================================

// The ResonanceGmZ class handles the gamma*/Z0 resonance.

class ResonanceGmZ : public ResonanceWidths {

public:

  // Constructor.
  ResonanceGmZ(int idResIn) : gmZmode(), thetaWRat(), ei2(), eivi(), vi2ai2(),
    gamNorm(), intNorm(), resNorm() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  int    gmZmode;
  double thetaWRat, ei2, eivi, vi2ai2, gamNorm, intNorm, resNorm;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

};

//==========================================================================

// The ResonanceW class handles the W+- resonance.

class ResonanceW : public ResonanceWidths {

public:

  // Constructor.
  ResonanceW(int idResIn) : thetaWRat(), alpEM() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double thetaWRat, alpEM;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceTop class handles the top/antitop resonance.

class ResonanceTop : public ResonanceWidths {

public:

  // Constructor.
  ResonanceTop(int idResIn) : thetaWRat(), m2W(), tanBeta(), tan2Beta(),
    mbRun() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double thetaWRat, m2W, tanBeta, tan2Beta, mbRun;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceFour class handles fourth-generation resonances.

class ResonanceFour : public ResonanceWidths {

public:

  // Constructor.
  ResonanceFour(int idResIn) : thetaWRat(), m2W() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double thetaWRat, m2W;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceH class handles the SM and BSM Higgs resonance.
// higgsType = 0 : SM H; = 1: h^0/H_1; = 2 : H^0/H_2; = 3 : A^0/A_3.

class ResonanceH : public ResonanceWidths {

public:

  // Constructor.
  ResonanceH(int higgsTypeIn, int idResIn) : higgsType(higgsTypeIn),
    useCubicWidth(), useRunLoopMass(), useNLOWidths(), sin2tW(), cos2tW(),
    mT(), mZ(), mW(), mHchg(), GammaT(), GammaZ(), GammaW(), rescAlpS(),
    rescColQ(), coup2d(), coup2u(), coup2l(), coup2Z(), coup2W(), coup2Hchg(),
    coup2H1H1(), coup2A3A3(), coup2H1Z(), coup2A3Z(), coup2A3H1(),
    coup2HchgW(), mLowT(), mStepT(), mLowZ(), mStepZ(), mLowW(), mStepW(),
    kinFacT(), kinFacZ(), kinFacW() {initBasic(idResIn);}

private:

  // Constants: could only be changed in the code itself.
  static const double MASSMINWZ, MASSMINT, GAMMAMARGIN;

  // Higgs type in current instance.
  int    higgsType;

  // Locally stored properties and couplings.
  bool   useCubicWidth, useRunLoopMass, useNLOWidths;
  double sin2tW, cos2tW, mT, mZ, mW, mHchg, GammaT, GammaZ, GammaW,
         rescAlpS, rescColQ, coup2d, coup2u, coup2l, coup2Z, coup2W,
         coup2Hchg, coup2H1H1, coup2A3A3, coup2H1Z, coup2A3Z, coup2A3H1,
         coup2HchgW, mLowT, mStepT, mLowZ, mStepZ, mLowW, mStepW,
         kinFacT[101], kinFacZ[101], kinFacW[101];

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

  // Sum up loop contributions in Higgs -> g + g.
  double eta2gg();

  // Sum up loop contributions in Higgs -> gamma + gamma.
  double eta2gaga();

  // Sum up loop contributions in Higgs -> gamma + Z0.
  double eta2gaZ();

};

//==========================================================================

// The ResonanceHchg class handles the H+- resonance.

class ResonanceHchg : public ResonanceWidths {

public:

  // Constructor.
  ResonanceHchg(int idResIn) : useCubicWidth(), thetaWRat(), mW(), tanBeta(),
    tan2Beta(), coup2H1W() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  bool   useCubicWidth;
  double thetaWRat, mW, tanBeta, tan2Beta, coup2H1W;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceZprime class handles the gamma*/Z0 /Z'^0 resonance.

class ResonanceZprime : public ResonanceWidths {

public:

  // Constructor.
  ResonanceZprime(int idResIn) : gmZmode(), maxZpGen(), sin2tW(), cos2tW(),
    thetaWRat(), mZ(), GammaZ(), m2Z(), GamMRatZ(), afZp(), vfZp(), coupZpWW(),
    ei2(), eivi(), vai2(), eivpi(), vaivapi(), vapi2(), gamNorm(), gamZNorm(),
    ZNorm(), gamZpNorm(), ZZpNorm(), ZpNorm() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  int    gmZmode, maxZpGen;
  double sin2tW, cos2tW, thetaWRat, mZ, GammaZ, m2Z, GamMRatZ, afZp[20],
         vfZp[20], coupZpWW, ei2, eivi, vai2, eivpi, vaivapi, vapi2,
         gamNorm, gamZNorm, ZNorm, gamZpNorm, ZZpNorm, ZpNorm;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

};

//==========================================================================

// The ResonanceWprime class handles the W'+- resonance.

class ResonanceWprime : public ResonanceWidths {

public:

  // Constructor.
  ResonanceWprime(int idResIn) : thetaWRat(), cos2tW(), alpEM(), aqWp(),
    vqWp(), alWp(), vlWp(), coupWpWZ() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double thetaWRat, cos2tW, alpEM, aqWp, vqWp, alWp, vlWp, coupWpWZ;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceRhorizontal class handles the R^0 resonance.

class ResonanceRhorizontal : public ResonanceWidths {

public:

  // Constructor.
  ResonanceRhorizontal(int idResIn) : thetaWRat() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double thetaWRat;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceExcited class handles excited-fermion resonances.

class ResonanceExcited : public ResonanceWidths {

public:

  // Constructor.
  ResonanceExcited(int idResIn) : Lambda(), coupF(), coupFprime(), coupFcol(),
    contactDec(), sin2tW(), cos2tW() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double Lambda, coupF, coupFprime, coupFcol, contactDec, sin2tW, cos2tW;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceGraviton class handles the excited Graviton resonance.

class ResonanceGraviton : public ResonanceWidths {

public:

  // Constructor.
  ResonanceGraviton(int idResIn) : eDsmbulk(), eDvlvl(), kappaMG(),
    eDcoupling() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  bool   eDsmbulk, eDvlvl;
  double kappaMG;

  // Couplings between graviton and SM (map from particle id to coupling).
  double eDcoupling[27];

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceKKgluon class handles the g^*/KK-gluon^* resonance.

class ResonanceKKgluon : public ResonanceWidths {

public:

  // Constructor.
  ResonanceKKgluon(int idResIn) : normSM(), normInt(), normKK(), eDgv(),
    eDga(), interfMode() {initBasic(idResIn);}

private:

  // Locally stored properties.
  double normSM, normInt, normKK;

  // Couplings between kk gluon and SM (indexed by particle id).
  // Helicity dependent couplings. Use vector/axial-vector
  // couplings internally, gv/ga = 0.5 * (gL +/- gR).
  double eDgv[10], eDga[10];

  // Interference parameter.
  int interfMode;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool calledFromInit = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

};

//==========================================================================

// The ResonanceLeptoquark class handles the LQ/LQbar resonance.

class ResonanceLeptoquark : public ResonanceWidths {

public:

  // Constructor.
  ResonanceLeptoquark(int idResIn) : kCoup() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double kCoup;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceNuRight class handles righthanded Majorana neutrinos.

class ResonanceNuRight : public ResonanceWidths {

public:

  // Constructor.
  ResonanceNuRight(int idResIn) : thetaWRat(), mWR() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double thetaWRat, mWR;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceZRight class handles the Z_R^0 resonance.

class ResonanceZRight : public ResonanceWidths {

public:

  // Constructor.
  ResonanceZRight(int idResIn) : sin2tW(), thetaWRat() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double sin2tW, thetaWRat;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceWRight class handles the W_R+- resonance.

class ResonanceWRight : public ResonanceWidths {

public:

  // Constructor.
  ResonanceWRight(int idResIn) : thetaWRat() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double thetaWRat;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceHchgchgLeft class handles the H++/H-- (left) resonance.

class ResonanceHchgchgLeft : public ResonanceWidths {

public:

  // Constructor.
  ResonanceHchgchgLeft(int idResIn) : yukawa(), gL(), vL(),
    mW() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  double yukawa[4][4], gL, vL, mW;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

// The ResonanceHchgchgRight class handles the H++/H-- (right) resonance.

class ResonanceHchgchgRight : public ResonanceWidths {

public:

  // Constructor.
  ResonanceHchgchgRight(int idResIn) : idWR(), yukawa(),
    gR() {initBasic(idResIn);}

private:

  // Locally stored properties and couplings.
  int    idWR;
  double yukawa[4][4], gR;

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool = false);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ResonanceWidths_H
