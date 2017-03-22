// ParticleData.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// DecayChannel, ParticleDataEntry and ParticleData classes.

#include "Pythia8/ParticleData.h"
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyResonanceWidths.h"

// Allow string and character manipulation.
#include <cctype>

namespace Pythia8 {

//==========================================================================

// DecayChannel class.
// This class holds info on a single decay channel.

//--------------------------------------------------------------------------

// Check whether id1 occurs anywhere in product list.

bool DecayChannel::contains(int id1) const {

  bool found1 = false;
  for (int i = 0; i < nProd; ++i) if (prod[i] == id1) found1 = true;
  return found1;

}

//--------------------------------------------------------------------------

// Check whether id1 and id2 occur anywhere in product list.
// iF ID1 == ID2 then two copies of this particle must be present.

bool DecayChannel::contains(int id1, int id2) const {

  bool found1 = false;
  bool found2 = false;
  for (int i = 0; i < nProd; ++i) {
    if (!found1 && prod[i] == id1) {found1 = true; continue;}
    if (!found2 && prod[i] == id2) {found2 = true; continue;}
  }
  return found1 && found2;

}

//--------------------------------------------------------------------------

// Check whether id1, id2 and id3 occur anywhere in product list.
// iF ID1 == ID2 then two copies of this particle must be present, etc.

bool DecayChannel::contains(int id1, int id2, int id3) const {

  bool found1 = false;
  bool found2 = false;
  bool found3 = false;
  for (int i = 0; i < nProd; ++i) {
    if (!found1 && prod[i] == id1) {found1 = true; continue;}
    if (!found2 && prod[i] == id2) {found2 = true; continue;}
    if (!found3 && prod[i] == id3) {found3 = true; continue;}
  }
  return found1 && found2 && found3;

}

//==========================================================================

// ParticleDataEntry class.
// This class holds info on a single particle species.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// A particle is invisible if it has neither strong nor electric charge,
// and is not made up of constituents that have it. Only relevant for
// long-lived particles. This list may need to be extended.
const int ParticleDataEntry::INVISIBLENUMBER = 52;
const int ParticleDataEntry::INVISIBLETABLE[52] = {
       12,      14,      16,      18,      23,      25,      32,      33,
       35,      36,      39,      41,      45,      46, 1000012, 1000014,
  1000016, 1000018, 1000022, 1000023, 1000025, 1000035, 1000045, 1000039,
  2000012, 2000014, 2000016, 2000018, 4900012, 4900014, 4900016, 4900021,
  4900022, 4900101, 4900102, 4900103, 4900104, 4900105, 4900106, 4900107,
  4900108, 4900111, 4900113, 4900211, 4900213, 4900991, 5000039, 5100039,
  9900012, 9900014, 9900016, 9900023 };

// For some particles we know it is necessary to switch off width,
// although they do have one, so do not warn.
const int ParticleDataEntry::KNOWNNOWIDTH[3] = {10313, 10323, 10333};

// Particles with a read-in tau0 (in mm/c) below this mayDecay by default.
const double ParticleDataEntry::MAXTAU0FORDECAY = 1000.;

// Particles with a read-in m0 above this isResonance by default.
const double ParticleDataEntry::MINMASSRESONANCE = 20.;

// Narrow states are assigned nominal mass.
const double ParticleDataEntry::NARROWMASS       = 1e-6;

// Constituent quark masses (d, u, s, c, b, -, -, -, g).
const double ParticleDataEntry::CONSTITUENTMASSTABLE[10]
  = {0., 0.325, 0.325, 0.50, 1.60, 5.00, 0., 0., 0., 0.7};

//--------------------------------------------------------------------------

// Destructor: delete any ResonanceWidths object.

ParticleDataEntry::~ParticleDataEntry() {
  if (resonancePtr != 0) delete resonancePtr;
}

//--------------------------------------------------------------------------

// Set initial default values for some quantities.

void ParticleDataEntry::setDefaults() {

  // A particle is a resonance if it is heavy enough.
  isResonanceSave     = (m0Save > MINMASSRESONANCE);

  // A particle may decay if it is shortlived enough.
  mayDecaySave        = (tau0Save < MAXTAU0FORDECAY);

  // A particle by default has no external decays.
  doExternalDecaySave = false;

  // A particle is invisible if in current table of such.
  isVisibleSave = true;
  for (int i = 0; i < INVISIBLENUMBER; ++i)
    if (idSave == INVISIBLETABLE[i]) isVisibleSave = false;

  // Normally a resonance should not have width forced to fixed value.
  doForceWidthSave  = false;

  // Set up constituent masses.
  setConstituentMass();

  // No Breit-Wigner mass selection before initialized.
  modeBWnow = 0;

}

//--------------------------------------------------------------------------

// Find out if a particle is a hadron.
// Only covers normal hadrons, not e.g. R-hadrons.

bool ParticleDataEntry::isHadron() const {

  if (idSave <= 100 || (idSave >= 1000000 && idSave <= 9000000)
    || idSave >= 9900000) return false;
  if (idSave == 130 || idSave == 310) return true;
  if (idSave%10 == 0 || (idSave/10)%10 == 0 || (idSave/100)%10 == 0)
    return false;
  return true;

}

//--------------------------------------------------------------------------

// Find out if a particle is a meson.
// Only covers normal hadrons, not e.g. R-hadrons.

bool ParticleDataEntry::isMeson() const {

  if (idSave <= 100 || (idSave >= 1000000 && idSave <= 9000000)
    || idSave >= 9900000) return false;
  if (idSave == 130 || idSave == 310) return true;
  if (idSave%10 == 0 || (idSave/10)%10 == 0 || (idSave/100)%10 == 0
    || (idSave/1000)%10 != 0) return false;
  return true;

}

//--------------------------------------------------------------------------

// Find out if a particle is a baryon.
// Only covers normal hadrons, not e.g. R-hadrons.

bool ParticleDataEntry::isBaryon() const {

  if (idSave <= 1000 || (idSave >= 1000000 && idSave <= 9000000)
    || idSave >= 9900000) return false;
  if (idSave%10 == 0 || (idSave/10)%10 == 0 || (idSave/100)%10 == 0
    || (idSave/1000)%10 == 0) return false;
  return true;


}

//--------------------------------------------------------------------------

// Extract the heaviest (= largest id)  quark in a hadron.

int ParticleDataEntry::heaviestQuark(int idIn) const {

  if (!isHadron()) return 0;
  int hQ = 0;

  // Meson.
  if ( (idSave/1000)%10 == 0 ) {
    hQ = (idSave/100)%10;
    if (idSave == 130) hQ = 3;
    if (hQ%2 == 1) hQ = -hQ;

  // Baryon.
  } else hQ = (idSave/1000)%10;

  // Done.
  return (idIn > 0) ? hQ : -hQ;

}

//--------------------------------------------------------------------------

// Calculate three times baryon number, i.e. net quark - antiquark number.

int ParticleDataEntry::baryonNumberType(int idIn) const {

  // Quarks.
  if (isQuark()) return (idIn > 0) ? 1 : -1;

  // Diquarks
  if (isDiquark()) return (idIn > 0) ? 2 : -2;

  // Baryons.
  if (isBaryon()) return (idIn > 0) ? 3 : -3;

  // Done.
  return 0;

}

//--------------------------------------------------------------------------

// Prepare the Breit-Wigner mass selection by precalculating
// frequently-used expressions.

void ParticleDataEntry::initBWmass() {

  // Find Breit-Wigner mode for current particle.
  modeBWnow = particleDataPtr->modeBreitWigner;
  if ( m0Save < NARROWMASS ) mWidthSave = 0.;
  if ( mWidthSave < NARROWMASS || (mMaxSave > mMinSave
    && mMaxSave - mMinSave < NARROWMASS) ) modeBWnow = 0;
  if (modeBWnow == 0) return;

  // Find atan expressions to be used in random mass selection.
  if (modeBWnow < 3) {
    atanLow = atan( 2. * (mMinSave - m0Save) / mWidthSave );
    double atanHigh = (mMaxSave > mMinSave)
      ? atan( 2. * (mMaxSave - m0Save) / mWidthSave ) : 0.5 * M_PI;
    atanDif = atanHigh - atanLow;
  } else {
    atanLow = atan( (pow2(mMinSave) - pow2(m0Save))
      / (m0Save * mWidthSave) );
    double atanHigh = (mMaxSave > mMinSave)
      ? atan( (pow2(mMaxSave) - pow2(m0Save)) / (m0Save * mWidthSave) )
      : 0.5 * M_PI;
    atanDif = atanHigh - atanLow;
  }

  // Done if no threshold factor.
  if (modeBWnow%2 == 1) return;

  // Find average mass threshold for threshold-factor correction.
  double bRatSum = 0.;
  double mThrSum = 0;
  for (int i = 0; i < int(channels.size()); ++i)
  if (channels[i].onMode() > 0) {
    bRatSum += channels[i].bRatio();
    double mChannelSum = 0.;
    for (int j = 0; j < channels[i].multiplicity(); ++j)
      mChannelSum += particleDataPtr->m0( channels[i].product(j) );
    mThrSum += channels[i].bRatio() * mChannelSum;
  }
  mThr = (bRatSum == 0.) ? 0. : mThrSum / bRatSum;

  // Switch off Breit-Wigner if very close to threshold.
  if (mThr + NARROWMASS > m0Save) {
    modeBWnow = 0;
    bool knownProblem = false;
    for (int i = 0; i < 3; ++i) if (idSave == KNOWNNOWIDTH[i])
      knownProblem = true;
    if (!knownProblem) {
      ostringstream osWarn;
      osWarn << "for id = " << idSave;
      particleDataPtr->infoPtr->errorMsg("Warning in ParticleDataEntry::"
        "initBWmass: switching off width", osWarn.str(), true);
    }
  }

}

//--------------------------------------------------------------------------

// Function to give mass of a particle, either at the nominal value
// or picked according to a (linear or quadratic) Breit-Wigner.

double ParticleDataEntry::mSel() {

  // Nominal value. (Width check should not be needed, but just in case.)
  if (modeBWnow == 0 || mWidthSave < NARROWMASS) return m0Save;
  double mNow, m2Now;

  // Mass according to a Breit-Wigner linear in m.
  if (modeBWnow == 1) {
     mNow = m0Save + 0.5 * mWidthSave
       * tan( atanLow + atanDif * particleDataPtr->rndmPtr->flat() );

  // Ditto, but make Gamma proportional to sqrt(m^2 - m_threshold^2).
  } else if (modeBWnow == 2) {
    double mWidthNow, fixBW, runBW;
    double m0ThrS = m0Save*m0Save - mThr*mThr;
    do {
      mNow = m0Save + 0.5 * mWidthSave
        * tan( atanLow + atanDif * particleDataPtr->rndmPtr->flat() );
      mWidthNow = mWidthSave * sqrtpos( (mNow*mNow - mThr*mThr) / m0ThrS );
      fixBW = mWidthSave / (pow2(mNow - m0Save) + pow2(0.5 * mWidthSave));
      runBW = mWidthNow / (pow2(mNow - m0Save) + pow2(0.5 * mWidthNow));
    } while (runBW < particleDataPtr->rndmPtr->flat()
      * particleDataPtr->maxEnhanceBW * fixBW);

  // Mass according to a Breit-Wigner quadratic in m.
  } else if (modeBWnow == 3) {
    m2Now = m0Save*m0Save + m0Save * mWidthSave
      * tan( atanLow + atanDif * particleDataPtr->rndmPtr->flat() );
    mNow = sqrtpos( m2Now);

  // Ditto, but m_0 Gamma_0 -> m Gamma(m) with threshold factor as above.
  } else {
    double mwNow, fixBW, runBW;
    double m2Ref = m0Save * m0Save;
    double mwRef = m0Save * mWidthSave;
    double m2Thr = mThr * mThr;
    do {
      m2Now = m2Ref + mwRef * tan( atanLow + atanDif
        * particleDataPtr->rndmPtr->flat() );
      mNow = sqrtpos( m2Now);
      mwNow = mNow * mWidthSave
        * sqrtpos( (m2Now - m2Thr) / (m2Ref - m2Thr) );
      fixBW = mwRef / (pow2(m2Now - m2Ref) + pow2(mwRef));
      runBW = mwNow / (pow2(m2Now - m2Ref) + pow2(mwNow));
    } while (runBW < particleDataPtr->rndmPtr->flat()
      * particleDataPtr->maxEnhanceBW * fixBW);
  }

  // Done.
  return mNow;
}

//--------------------------------------------------------------------------

// Function to calculate running mass at given mass scale.

double ParticleDataEntry::mRun(double mHat) {

  // Except for six quarks return nominal mass.
  if (idSave > 6) return m0Save;
  double mQRun = particleDataPtr->mQRun[idSave];
  double Lam5  = particleDataPtr->Lambda5Run;

  // For d, u, s quarks start running at 2 GeV (RPP 2006 p. 505).
  if (idSave < 4) return mQRun * pow ( log(2. / Lam5)
    / log(max(2., mHat) / Lam5), 12./23.);

  // For c, b and t quarks start running at respective mass.
  return mQRun * pow ( log(mQRun / Lam5)
    / log(max(mQRun, mHat) / Lam5), 12./23.);

}

//--------------------------------------------------------------------------

// Rescale all branching ratios to assure normalization to unity.

void ParticleDataEntry::rescaleBR(double newSumBR) {

  // Sum up branching ratios. Find rescaling factor. Rescale.
  double oldSumBR = 0.;
  for ( int i = 0; i < int(channels.size()); ++ i)
    oldSumBR += channels[i].bRatio();
  double rescaleFactor = newSumBR / oldSumBR;
  for ( int i = 0; i < int(channels.size()); ++ i)
    channels[i].rescaleBR(rescaleFactor);

}

//--------------------------------------------------------------------------

// Prepare to pick a decay channel.

bool ParticleDataEntry::preparePick(int idSgn, double mHat, int idInFlav) {

  // Reset sum of allowed widths/branching ratios.
  currentBRSum = 0.;

  // For resonances the widths are calculated dynamically.
  if (isResonanceSave && resonancePtr != 0) {
    resonancePtr->widthStore(idSgn, mHat, idInFlav);
    for (int i = 0; i < int(channels.size()); ++i)
      currentBRSum += channels[i].currentBR();

  // Else use normal fixed branching ratios.
  } else {
    int onMode;
    double currentBRNow;
    for (int i = 0; i < int(channels.size()); ++i) {
      onMode = channels[i].onMode();
      currentBRNow = 0.;
      if ( idSgn > 0 && (onMode == 1 || onMode == 2) )
        currentBRNow = channels[i].bRatio();
      else if ( idSgn < 0 && (onMode == 1 || onMode == 3) )
        currentBRNow = channels[i].bRatio();
      channels[i].currentBR(currentBRNow);
      currentBRSum += currentBRNow;
    }
  }

  // Failure if no channels found with positive branching ratios.
  return (currentBRSum > 0.);

}

//--------------------------------------------------------------------------

// Pick a decay channel according to branching ratios from preparePick.

DecayChannel& ParticleDataEntry::pickChannel() {

  // Find channel in table.
  int size = channels.size();
  double rndmBR = currentBRSum * particleDataPtr->rndmPtr->flat();
  int i = -1;
  do rndmBR -= channels[++i].currentBR();
  while (rndmBR > 0. && i < size);

  // Emergency if no channel found. Done.
  if (i == size) i = 0;
  return channels[i];

}

//--------------------------------------------------------------------------

// Access methods stored in ResonanceWidths. Could have been
// inline in .h, except for problems with forward declarations.

void ParticleDataEntry::setResonancePtr(
  ResonanceWidths* resonancePtrIn) {
  if (resonancePtr == resonancePtrIn) return;
  if (resonancePtr != 0) delete resonancePtr;
  resonancePtr = resonancePtrIn;
}

void ParticleDataEntry::resInit(Info* infoPtrIn, Settings* settingsPtrIn,
  ParticleData* particleDataPtrIn, Couplings* couplingsPtrIn) {
  if (resonancePtr != 0) resonancePtr->init(infoPtrIn, settingsPtrIn,
  particleDataPtrIn, couplingsPtrIn);
}

double ParticleDataEntry::resWidth(int idSgn, double mHat, int idIn,
  bool openOnly, bool setBR) {
  return (resonancePtr != 0) ? resonancePtr->width( idSgn, mHat,
    idIn, openOnly, setBR) : 0.;
}

double ParticleDataEntry::resWidthOpen(int idSgn, double mHat, int idIn) {
  return (resonancePtr != 0) ? resonancePtr->widthOpen( idSgn, mHat, idIn)
  : 0.;
}

double ParticleDataEntry::resWidthStore(int idSgn, double mHat, int idIn) {
  return (resonancePtr != 0) ? resonancePtr->widthStore( idSgn, mHat, idIn)
  : 0.;
}

double ParticleDataEntry::resOpenFrac(int idSgn) {
  return (resonancePtr != 0) ? resonancePtr->openFrac(idSgn) : 1.;
}

double ParticleDataEntry::resWidthRescaleFactor() {
  return (resonancePtr != 0) ? resonancePtr->widthRescaleFactor() : 1.;
}

double ParticleDataEntry::resWidthChan(double mHat, int idAbs1,
  int idAbs2) {
  return (resonancePtr != 0) ? resonancePtr->widthChan( mHat, idAbs1,
    idAbs2) : 0.;
}

//--------------------------------------------------------------------------

// Constituent masses for (d, u, s, c, b) quarks and diquarks.
// Hardcoded in CONSTITUENTMASSTABLE so that they are not overwritten
// by mistake, and separated from the "normal" masses.
// Called both by setDefaults and setM0 so kept as separate method.

void ParticleDataEntry::setConstituentMass() {

  // Equate with the normal masses as default guess.
  constituentMassSave = m0Save;

  // Quark masses trivial. Also gluon mass.
  if (idSave < 6) constituentMassSave = CONSTITUENTMASSTABLE[idSave];
  if (idSave == 21) constituentMassSave = CONSTITUENTMASSTABLE[9];

  // Diquarks as simple sum of constituent quarks.
  if (idSave > 1000 && idSave < 10000 && (idSave/10)%10 == 0) {
    int id1 = idSave/1000;
    int id2 = (idSave/100)%10;
    if (id1 <6 && id2 < 6) constituentMassSave
      = CONSTITUENTMASSTABLE[id1] + CONSTITUENTMASSTABLE[id2];
  }

}

//==========================================================================

// ParticleData class.
// This class holds a map of all ParticleDataEntries,
// each entry containing info on a particle species.

//--------------------------------------------------------------------------

// Get data to be distributed among particles during setup.
// Note: this routine is called twice. Firstly from init(...), but
// the data should not be used at that point, so is likely overkill.
// Secondly, from initWidths, after user has had time to change.

void ParticleData::initCommon() {

  // Mass generation: fixed mass or linear/quadratic Breit-Wigner.
  modeBreitWigner = settingsPtr->mode("ParticleData:modeBreitWigner");

  // Maximum tail enhancement when adding threshold factor to Breit-Wigner.
  maxEnhanceBW    = settingsPtr->parm("ParticleData:maxEnhanceBW");

  // Find initial MSbar masses for five light flavours.
  mQRun[1]        = settingsPtr->parm("ParticleData:mdRun");
  mQRun[2]        = settingsPtr->parm("ParticleData:muRun");
  mQRun[3]        = settingsPtr->parm("ParticleData:msRun");
  mQRun[4]        = settingsPtr->parm("ParticleData:mcRun");
  mQRun[5]        = settingsPtr->parm("ParticleData:mbRun");
  mQRun[6]        = settingsPtr->parm("ParticleData:mtRun");

  // Find Lambda5 value to use in running of MSbar masses.
  double alphaSvalue = settingsPtr->parm("ParticleData:alphaSvalueMRun");
  AlphaStrong alphaS;
  alphaS.init( alphaSvalue, 1, 5, false);
  Lambda5Run = alphaS.Lambda5();

}

//--------------------------------------------------------------------------

// Initialize pointer for particles to the full database, the Breit-Wigners
// of normal hadrons and the ResonanceWidths of resonances. For the latter
// the order of initialization is essential to get secondary widths right.

void ParticleData::initWidths( vector<ResonanceWidths*> resonancePtrs) {

  // Initialize some common data.
  initCommon();

  // Pointer to database and Breit-Wigner mass initialization for each
  // particle.
  ResonanceWidths* resonancePtr = 0;
  for (map<int, ParticleDataEntry>::iterator pdtEntry = pdt.begin();
    pdtEntry != pdt.end(); ++pdtEntry) {
    ParticleDataEntry& pdtNow = pdtEntry->second;
    pdtNow.initBWmass();

    // Remove any existing resonances.
    resonancePtr = pdtNow.getResonancePtr();
    if (resonancePtr != 0) pdtNow.setResonancePtr(0);
  }

  // Begin set up new resonance objects.
  // Z0, W+- and top are almost always needed.
  resonancePtr = new ResonanceGmZ(23);
  setResonancePtr( 23, resonancePtr);
  resonancePtr = new ResonanceW(24);
  setResonancePtr( 24, resonancePtr);
  resonancePtr = new ResonanceTop(6);
  setResonancePtr(  6, resonancePtr);

  // Higgs in SM.
  if (!settingsPtr->flag("Higgs:useBSM")) {
    resonancePtr = new ResonanceH(0, 25);
    setResonancePtr( 25, resonancePtr);

  // Higgses in BSM.
  } else {
    resonancePtr = new ResonanceH(1, 25);
    setResonancePtr( 25, resonancePtr);
    resonancePtr = new ResonanceH(2, 35);
    setResonancePtr( 35, resonancePtr);
    resonancePtr = new ResonanceH(3, 36);
    setResonancePtr( 36, resonancePtr);
    resonancePtr = new ResonanceHchg(37);
    setResonancePtr( 37, resonancePtr);
    resonancePtr = new ResonanceH(4, 45);
    setResonancePtr( 45, resonancePtr);
    resonancePtr = new ResonanceH(5, 46);
    setResonancePtr( 46, resonancePtr);
  }

  // A fourth generation: b', t', tau', nu'_tau.
  resonancePtr = new ResonanceFour(7);
  setResonancePtr( 7, resonancePtr);
  resonancePtr = new ResonanceFour(8);
  setResonancePtr( 8, resonancePtr);
  resonancePtr = new ResonanceFour(17);
  setResonancePtr( 17, resonancePtr);
  resonancePtr = new ResonanceFour(18);
  setResonancePtr( 18, resonancePtr);

  // New gauge bosons: Z', W', R.
  resonancePtr = new ResonanceZprime(32);
  setResonancePtr( 32, resonancePtr);
  resonancePtr = new ResonanceWprime(34);
  setResonancePtr( 34, resonancePtr);
  resonancePtr = new ResonanceRhorizontal(41);
  setResonancePtr( 41, resonancePtr);

  // A leptoquark.
  resonancePtr = new ResonanceLeptoquark(42);
  setResonancePtr( 42, resonancePtr);

  // 93 = Z0copy and 94 = W+-copy used to pick decay channels
  // for W/Z production in parton showers.
  resonancePtr = new ResonanceGmZ(93);
  setResonancePtr( 93, resonancePtr);
  resonancePtr = new ResonanceW(94);
  setResonancePtr( 94, resonancePtr);

  // Supersymmetry
  //  - Squarks
  for(int i = 1; i < 7; i++){
    resonancePtr = new ResonanceSquark(1000000 + i);
    setResonancePtr( 1000000 + i, resonancePtr);
    resonancePtr = new ResonanceSquark(2000000 + i);
    setResonancePtr( 2000000 + i, resonancePtr);
  }

  //  - Sleptons and sneutrinos
  for(int i = 1; i < 7; i++){
    resonancePtr = new ResonanceSlepton(1000010 + i);
    setResonancePtr( 1000010 + i, resonancePtr);
    resonancePtr = new ResonanceSlepton(2000010 + i);
    setResonancePtr( 2000010 + i, resonancePtr);
  }

  // - Gluino
  resonancePtr = new ResonanceGluino(1000021);
  setResonancePtr( 1000021, resonancePtr);

  // - Charginos
  resonancePtr = new ResonanceChar(1000024);
  setResonancePtr( 1000024, resonancePtr);
  resonancePtr = new ResonanceChar(1000037);
  setResonancePtr( 1000037, resonancePtr);

  // - Neutralinos
  resonancePtr = new ResonanceNeut(1000022);
  setResonancePtr( 1000022, resonancePtr);
  resonancePtr = new ResonanceNeut(1000023);
  setResonancePtr( 1000023, resonancePtr);
  resonancePtr = new ResonanceNeut(1000025);
  setResonancePtr( 1000025, resonancePtr);
  resonancePtr = new ResonanceNeut(1000035);
  setResonancePtr( 1000035, resonancePtr);
  resonancePtr = new ResonanceNeut(1000045);
  setResonancePtr( 1000045, resonancePtr);

  // Excited quarks and leptons.
  for (int i = 1; i < 7; ++i) {
    resonancePtr = new ResonanceExcited(4000000 + i);
    setResonancePtr( 4000000 + i, resonancePtr);
  }
  for (int i = 11; i < 17; ++i) {
    resonancePtr = new ResonanceExcited(4000000 + i);
    setResonancePtr( 4000000 + i, resonancePtr);
  }

  // An excited graviton/gluon in extra-dimensional scenarios.
  resonancePtr = new ResonanceGraviton(5100039);
  setResonancePtr( 5100039, resonancePtr);
  resonancePtr = new ResonanceKKgluon(5100021);
  setResonancePtr( 5100021, resonancePtr);

  // A left-right-symmetric scenario with new righthanded neutrinos,
  // righthanded gauge bosons and doubly charged Higgses.
  resonancePtr = new ResonanceNuRight(9900012);
  setResonancePtr( 9900012, resonancePtr);
  resonancePtr = new ResonanceNuRight(9900014);
  setResonancePtr( 9900014, resonancePtr);
  resonancePtr = new ResonanceNuRight(9900016);
  setResonancePtr( 9900016, resonancePtr);
  resonancePtr = new ResonanceZRight(9900023);
  setResonancePtr( 9900023, resonancePtr);
  resonancePtr = new ResonanceWRight(9900024);
  setResonancePtr( 9900024, resonancePtr);
  resonancePtr = new ResonanceHchgchgLeft(9900041);
  setResonancePtr( 9900041, resonancePtr);
  resonancePtr = new ResonanceHchgchgRight(9900042);
  setResonancePtr( 9900042, resonancePtr);

  // Attach user-defined external resonances and do basic initialization.
  for (int i = 0; i < int(resonancePtrs.size()); ++i) {
    int idNow = resonancePtrs[i]->id();
    resonancePtrs[i]->initBasic(idNow);
    setResonancePtr( idNow, resonancePtrs[i]);
  }

  // Set up lists to order resonances in ascending mass.
  vector<int>    idOrdered;
  vector<double> m0Ordered;

  // Put Z0 and W+- first, since known to be SM and often off-shell.
  idOrdered.push_back(23);
  m0Ordered.push_back(m0(23));
  idOrdered.push_back(24);
  m0Ordered.push_back(m0(24));

  // Loop through particle table to find resonances.
  for (map<int, ParticleDataEntry>::iterator pdtEntry = pdt.begin();
    pdtEntry != pdt.end(); ++pdtEntry) {
    ParticleDataEntry& pdtNow = pdtEntry->second;
    int idNow = pdtNow.id();

    // Set up a simple default object for uninitialized resonances.
    if (pdtNow.isResonance() && pdtNow.getResonancePtr() == 0) {
      resonancePtr = new ResonanceGeneric(idNow);
      setResonancePtr( idNow, resonancePtr);
    }

    // Insert resonances in ascending mass, to respect decay hierarchies.
    if (pdtNow.getResonancePtr() != 0 && idNow != 23 && idNow != 24) {
      double m0Now = pdtNow.m0();
      idOrdered.push_back(idNow);
      m0Ordered.push_back(m0Now);
      for (int i = int(idOrdered.size()) - 2; i > 1; --i) {
        if (m0Ordered[i] < m0Now) break;
        swap( idOrdered[i], idOrdered[i + 1]);
        swap( m0Ordered[i], m0Ordered[i + 1]);
      }
    }
  }

  // Initialize the resonances in ascending mass order. Reset mass generation.
  for (int i = 0; i < int(idOrdered.size()); ++i) {
    resInit( idOrdered[i]);
    ParticleDataEntry* pdtPtrNow = particleDataEntryPtr( idOrdered[i]);
    pdtPtrNow->initBWmass();
  }

}

//--------------------------------------------------------------------------

// Read in database from specific XML file (which may refer to others).

bool ParticleData::readXML(string inFile, bool reset) {

  // Normally reset whole database before beginning.
  if (reset) {pdt.clear(); isInit = false;}

  // List of files to be checked.
  vector<string> files;
  files.push_back(inFile);

  // Loop over files. Open them for read.
  for (int i = 0; i < int(files.size()); ++i) {
    const char* cstring = files[i].c_str();
    ifstream is(cstring);

    // Check that instream is OK.
    if (!is.good()) {
      infoPtr->errorMsg("Error in ParticleData::readXML:"
        " did not find file", files[i]);
      return false;
    }

    // Read in one line at a time.
    particlePtr = 0;
    string line;
    while ( getline(is, line) ) {

      // Get first word of a line.
      istringstream getfirst(line);
      string word1;
      getfirst >> word1;

      // Check for occurence of a particle. Add any continuation lines.
      if (word1 == "<particle") {
        while (line.find(">") == string::npos) {
          string addLine;
          getline(is, addLine);
          line += addLine;
        }

        // Read in particle properties.
        int idTmp          = intAttributeValue( line, "id");
        string nameTmp     = attributeValue( line, "name");
        string antiNameTmp = attributeValue( line, "antiName");
        if (antiNameTmp == "") antiNameTmp = "void";
        int spinTypeTmp    = intAttributeValue( line, "spinType");
        int chargeTypeTmp  = intAttributeValue( line, "chargeType");
        int colTypeTmp     = intAttributeValue( line, "colType");
        double m0Tmp       = doubleAttributeValue( line, "m0");
        double mWidthTmp   = doubleAttributeValue( line, "mWidth");
        double mMinTmp     = doubleAttributeValue( line, "mMin");
        double mMaxTmp     = doubleAttributeValue( line, "mMax");
        double tau0Tmp     = doubleAttributeValue( line, "tau0");

        // Erase if particle already exists.
        if (isParticle(idTmp)) pdt.erase(idTmp);

        // Store new particle. Save pointer, to be used for decay channels.
        addParticle( idTmp, nameTmp, antiNameTmp, spinTypeTmp, chargeTypeTmp,
          colTypeTmp, m0Tmp, mWidthTmp, mMinTmp, mMaxTmp, tau0Tmp);
        particlePtr = particleDataEntryPtr(idTmp);

      // Check for occurence of a decay channel. Add any continuation lines.
      } else if (word1 == "<channel") {
        while (line.find(">") == string::npos) {
          string addLine;
          getline(is, addLine);
          line += addLine;
        }

        // Read in channel properties - products so far only as a string.
        int onMode      = intAttributeValue( line, "onMode");
        double bRatio   = doubleAttributeValue( line, "bRatio");
        int meMode      = intAttributeValue( line, "meMode");
        string products = attributeValue( line, "products");

        // Read in decay products from stream. Must have at least one.
        istringstream prodStream(products);
        int prod0 = 0; int prod1 = 0; int prod2 = 0; int prod3 = 0;
        int prod4 = 0; int prod5 = 0; int prod6 = 0; int prod7 = 0;
        prodStream >> prod0 >> prod1 >> prod2 >> prod3 >> prod4 >> prod5
                   >> prod6 >> prod7;
        if (prod0 == 0) {
          infoPtr->errorMsg("Error in ParticleData::readXML:"
            " incomplete decay channel", line);
          return false;
        }

        // Store new channel (if particle already known).
        if (particlePtr == 0) {
          infoPtr->errorMsg("Error in ParticleData::readXML:"
            " orphan decay channel", line);
          return false;
        }
        particlePtr->addChannel(onMode, bRatio, meMode, prod0, prod1,
          prod2, prod3, prod4, prod5, prod6, prod7);

      // Check for occurence of a file also to be read.
      } else if (word1 == "<file") {
        string file = attributeValue(line, "name");
        if (file == "") {
          infoPtr->errorMsg("Error in ParticleData::readXML:"
            " skip unrecognized file name", line);
        } else files.push_back(file);
      }

    // End of loop over lines in input file and loop over files.
    };
  };

  // All particle data at this stage defines baseline original.
  if (reset) for (map<int, ParticleDataEntry>::iterator pdtEntry
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    particlePtr = &pdtEntry->second;
    particlePtr->setHasChanged(false);
  }

  // Done.
  isInit = true;
  return true;

}

//--------------------------------------------------------------------------

// Print out complete database in numerical order as an XML file.

void ParticleData::listXML(string outFile) {

  // Convert file name to ofstream.
    const char* cstring = outFile.c_str();
    ofstream os(cstring);

  // Iterate through the particle data table.
  for (map<int, ParticleDataEntry>::iterator pdtEntry
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    particlePtr = &pdtEntry->second;

    // Print particle properties.
    os << "<particle id=\"" << particlePtr->id() << "\""
       << " name=\"" << particlePtr->name() << "\"";
    if (particlePtr->hasAnti())
      os << " antiName=\"" << particlePtr->name(-1) << "\"";
    os << " spinType=\"" << particlePtr->spinType() << "\""
       << " chargeType=\"" << particlePtr->chargeType() << "\""
       << " colType=\"" << particlePtr->colType() << "\"\n";
    // Pick format for mass and width based on mass value.
    double m0Now = particlePtr->m0();
    if (m0Now == 0 || (m0Now > 0.1 && m0Now < 1000.))
      os << fixed << setprecision(5);
    else  os << scientific << setprecision(3);
    os << "          m0=\"" << m0Now << "\"";
    if (particlePtr->mWidth() > 0.)
      os << " mWidth=\"" << particlePtr->mWidth() << "\""
         << " mMin=\"" << particlePtr->mMin() << "\""
         << " mMax=\"" << particlePtr->mMax() << "\"";
    if (particlePtr->tau0() > 0.) os << scientific << setprecision(5)
         << " tau0=\"" << particlePtr->tau0() << "\"";
    os << ">\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->sizeChannels() > 0) {
      for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
        const DecayChannel& channel = particlePtr->channel(i);
        int mult = channel.multiplicity();

        // Print decay channel properties.
        os << " <channel onMode=\"" << channel.onMode() << "\""
           << fixed << setprecision(7)
           << " bRatio=\"" << channel.bRatio() << "\"";
        if (channel.meMode() > 0)
          os << " meMode=\"" << channel.meMode() << "\"";
        os << " products=\"";
        for (int j = 0; j < mult; ++j) {
          os << channel.product(j);
          if (j < mult - 1) os << " ";
        }

        // Finish off line and loop over allowed decay channels.
        os  << "\"/>\n";
      }
    }

    // Finish off existing particle.
    os << "</particle>\n\n";

  }

}

//--------------------------------------------------------------------------

// Read in database from specific free format file.

bool ParticleData::readFF(string inFile, bool reset) {

  // Normally reset whole database before beginning.
  if (reset) {pdt.clear(); isInit = false;}

  // Open file for read and check that instream is OK.
  const char* cstring = inFile.c_str();
  ifstream is(cstring);
  if (!is.good()) {
    infoPtr->errorMsg("Error in ParticleData::readFF:"
      " did not find file", inFile);
    return false;
  }

  // Read in one line at a time.
  particlePtr = 0;
  string line;
  bool readParticle = false;
  while ( getline(is, line) ) {

    // Empty lines begins new particle.
    if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) {
      readParticle = true;
      continue;
    }

    // Prepare to use standard read from line.
    istringstream readLine(line);

    // Read in a line with particle information.
    if (readParticle) {

      // Properties to be read.
      int    idTmp;
      string nameTmp, antiNameTmp;
      int    spinTypeTmp, chargeTypeTmp, colTypeTmp;
      double m0Tmp, mWidthTmp, mMinTmp, mMaxTmp, tau0Tmp;
      string mayTmp;

      // Do the reading.
      readLine >> idTmp >> nameTmp >> antiNameTmp >> spinTypeTmp
               >> chargeTypeTmp >> colTypeTmp >> m0Tmp >> mWidthTmp
               >> mMinTmp >> mMaxTmp >> tau0Tmp;

      // Error printout if something went wrong.
      if (!readLine) {
        infoPtr->errorMsg("Error in ParticleData::readFF:"
          " incomplete particle", line);
        return false;
      }

      // Erase if particle already exists.
      if (isParticle(idTmp)) pdt.erase(idTmp);

      // Store new particle. Save pointer, to be used for decay channels.
      addParticle( idTmp, nameTmp, antiNameTmp, spinTypeTmp, chargeTypeTmp,
        colTypeTmp, m0Tmp, mWidthTmp, mMinTmp, mMaxTmp, tau0Tmp);
      particlePtr = particleDataEntryPtr(idTmp);
      readParticle = false;

    // Read in a line with decay channel information.
    } else {

      // Properties to be read.
      int    onMode = 0;
      double bRatio = 0.;
      int    meMode = 0;
      int prod0 = 0; int prod1 = 0; int prod2 = 0; int prod3 = 0;
      int prod4 = 0; int prod5 = 0; int prod6 = 0; int prod7 = 0;

      // Read in data from stream. Need at least one decay product.
      readLine >> onMode >> bRatio >> meMode >> prod0;
      if (!readLine) {
        infoPtr->errorMsg("Error in ParticleData::readFF:"
          " incomplete decay channel", line);
        return false;
      }
      readLine >> prod1 >> prod2 >> prod3 >> prod4 >> prod5
        >> prod6  >> prod7;

      // Store new channel.
      if (particlePtr == 0) {
        infoPtr->errorMsg("Error in ParticleData::readFF:"
          " orphan decay channel", line);
        return false;
      }
      particlePtr->addChannel(onMode, bRatio, meMode, prod0, prod1,
        prod2, prod3, prod4, prod5, prod6, prod7);

    }

  // End of while loop over lines in the file.
  }


  // Done.
  isInit = true;
  return true;

}

//--------------------------------------------------------------------------

// Print out complete database in numerical order as a free format file.

void ParticleData::listFF(string outFile) {

  // Convert file name to ofstream.
    const char* cstring = outFile.c_str();
    ofstream os(cstring);

  // Iterate through the particle data table.
  for (map<int, ParticleDataEntry>::iterator pdtEntry
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    particlePtr = &pdtEntry->second;

    // Pick format for mass and width based on mass value.
    double m0Now = particlePtr->m0();
    if (m0Now == 0 || (m0Now > 0.1 && m0Now < 1000.))
      os << fixed << setprecision(5);
    else os << scientific << setprecision(3);

    // Print particle properties.
    os << "\n" << setw(8) << particlePtr->id() << "  "
       << left << setw(16) << particlePtr->name() << " "
       << setw(16) << particlePtr->name(-1) << "  "
       << right << setw(2) << particlePtr->spinType() << "  "
       << setw(2) << particlePtr->chargeType() << "  "
       << setw(2) << particlePtr->colType() << " "
       << setw(10) << particlePtr->m0() << " "
       << setw(10) << particlePtr->mWidth() << " "
       << setw(10) << particlePtr->mMin() << " "
       << setw(10) << particlePtr->mMax() << " "
       << scientific << setprecision(5)
       << setw(12) << particlePtr->tau0() << "\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->sizeChannels() > 0) {
      for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
        const DecayChannel& channel = particlePtr->channel(i);
        os << "               " << setw(6) << channel.onMode()
           << "  " << fixed << setprecision(7) << setw(10)
           << channel.bRatio() << "  "
           << setw(3) << channel.meMode() << " ";
        for (int j = 0; j < channel.multiplicity(); ++j)
          os << setw(8) << channel.product(j) << " ";
        os << "\n";
      }
    }

  }

}

//--------------------------------------------------------------------------

// Read in updates from a character string, like a line of a file.
// Is used by readString (and readFile) in Pythia.

  bool ParticleData::readString(string lineIn, bool warn, ostream& os) {

  // If empty line then done.
  if (lineIn.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;

  // Take copy that will be modified.
  string line = lineIn;

  // If first character is not a digit then taken to be a comment.
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isdigit(line[firstChar])) return true;

  // Replace colons and equal signs by blanks to make parsing simpler.
  for ( int j = 0; j < int(line.size()); ++ j)
     if (line[j] == ':' || line[j] == '=') line[j] = ' ';

  // Get particle id and property name.
  int    idTmp;
  string property;
  istringstream getWord(line);
  getWord >> idTmp >> property;
  property = toLower(property);

  // Check that valid particle.
  if ( (!isParticle(idTmp) && property  != "all" && property  != "new")
  || idTmp <= 0) {
    if (warn) os << "\n PYTHIA Error: input particle not found in Particle"
      << " Data Table:\n   " << lineIn << "\n";
    readingFailedSave = true;
    return false;
  }

  // Identify particle property and read + set its value, case by case.
  if (property == "name") {
    string nameTmp;
    getWord >> nameTmp;
    pdt[idTmp].setName(nameTmp);
    return true;
  }
  if (property == "antiname") {
    string antiNameTmp;
    getWord >> antiNameTmp;
    pdt[idTmp].setAntiName(antiNameTmp);
    return true;
  }
  if (property == "names") {
    string nameTmp, antiNameTmp;
    getWord >> nameTmp >> antiNameTmp;
    pdt[idTmp].setNames(nameTmp, antiNameTmp);
    return true;
  }
  if (property == "spintype") {
    int spinTypeTmp;
    getWord >> spinTypeTmp;
    pdt[idTmp].setSpinType(spinTypeTmp);
    return true;
  }
  if (property == "chargetype") {
    int chargeTypeTmp;
    getWord >> chargeTypeTmp;
    pdt[idTmp].setChargeType(chargeTypeTmp);
    return true;
  }
  if (property == "coltype") {
    int colTypeTmp;
    getWord >> colTypeTmp;
    pdt[idTmp].setColType(colTypeTmp);
    return true;
  }
  if (property == "m0") {
    double m0Tmp;
    getWord >> m0Tmp;
    pdt[idTmp].setM0(m0Tmp);
    return true;
  }
  if (property == "mwidth") {
    double mWidthTmp;
    getWord >> mWidthTmp;
    pdt[idTmp].setMWidth(mWidthTmp);
    return true;
  }
  if (property == "mmin") {
    double mMinTmp;
    getWord >> mMinTmp;
    pdt[idTmp].setMMin(mMinTmp);
    return true;
  }
  if (property == "mmax") {
    double mMaxTmp;
    getWord >> mMaxTmp;
    pdt[idTmp].setMMax(mMaxTmp);
    return true;
  }
  if (property == "tau0") {
    double tau0Tmp;
    getWord >> tau0Tmp;
    pdt[idTmp].setTau0(tau0Tmp);
    return true;
  }
  if (property == "isresonance") {
    string isresTmp;
    getWord >> isresTmp;
    bool isResonanceTmp = boolString(isresTmp);
    pdt[idTmp].setIsResonance(isResonanceTmp);
    return true;
  }
  if (property == "maydecay") {
    string mayTmp;
    getWord >> mayTmp;
    bool mayDecayTmp = boolString(mayTmp);
    pdt[idTmp].setMayDecay(mayDecayTmp);
    return true;
  }
  if (property == "doexternaldecay") {
    string extdecTmp;
    getWord >> extdecTmp;
    bool doExternalDecayTmp = boolString(extdecTmp);
    pdt[idTmp].setDoExternalDecay(doExternalDecayTmp);
    return true;
  }
  if (property == "isvisible") {
    string isvisTmp;
    getWord >> isvisTmp;
    bool isVisibleTmp = boolString(isvisTmp);
    pdt[idTmp].setIsVisible(isVisibleTmp);
    return true;
  }
  if (property == "doforcewidth") {
    string doforceTmp;
    getWord >> doforceTmp;
    bool doForceWidthTmp = boolString(doforceTmp);
    pdt[idTmp].setDoForceWidth(doForceWidthTmp);
    return true;
  }

  // Addition or complete replacement of a particle.
  if (property == "all" || property == "new") {

    // Default values for properties to be read.
    string nameTmp       = "void";
    string antiNameTmp   = "void";
    int    spinTypeTmp   = 0;
    int    chargeTypeTmp = 0;
    int    colTypeTmp    = 0;
    double m0Tmp         = 0.;
    double mWidthTmp     = 0.;
    double mMinTmp       = 0.;
    double mMaxTmp       = 0.;
    double tau0Tmp       = 0.;

    // Read in data from stream.
    getWord >> nameTmp >> antiNameTmp >> spinTypeTmp >> chargeTypeTmp
            >> colTypeTmp >> m0Tmp >> mWidthTmp >> mMinTmp >> mMaxTmp
            >> tau0Tmp;

    // To keep existing decay channels, only overwrite particle data.
    if (property == "all" && isParticle(idTmp)) {
      setAll( idTmp, nameTmp, antiNameTmp, spinTypeTmp, chargeTypeTmp,
        colTypeTmp, m0Tmp, mWidthTmp, mMinTmp, mMaxTmp, tau0Tmp);

    // Else start over completely from scratch.
    } else {
      if (isParticle(idTmp)) pdt.erase(idTmp);
      addParticle( idTmp, nameTmp, antiNameTmp, spinTypeTmp, chargeTypeTmp,
        colTypeTmp, m0Tmp, mWidthTmp, mMinTmp, mMaxTmp, tau0Tmp);
    }
    return true;
  }

  // Set onMode of all decay channels in one go.
  if (property == "onmode") {
      int    onMode = 0;
      string onModeIn;
      getWord >> onModeIn;
      // For onMode allow the optional possibility of Bool input.
      if (isdigit(onModeIn[0])) {
        istringstream getOnMode(onModeIn);
        getOnMode >> onMode;
      } else onMode = (boolString(onModeIn)) ? 1 : 0;
    for (int i = 0; i < pdt[idTmp].sizeChannels(); ++i)
      pdt[idTmp].channel(i).onMode(onMode);
    return true;
  }

  // Selective search for matching decay channels.
  int matchKind = 0;
  if (property == "offifany" || property == "onifany" ||
    property == "onposifany" || property == "onnegifany") matchKind = 1;
  if (property == "offifall" || property == "onifall" ||
    property == "onposifall" || property == "onnegifall") matchKind = 2;
  if (property == "offifmatch" || property == "onifmatch" ||
    property == "onposifmatch" || property == "onnegifmatch") matchKind = 3;
  if (matchKind > 0) {
    int onMode = 0;
    if (property == "onifany" || property == "onifall"
      || property == "onifmatch") onMode = 1;
    if (property == "onposifany" || property == "onposifall"
      || property == "onposifmatch") onMode = 2;
    if (property == "onnegifany" || property == "onnegifall"
      || property == "onnegifmatch") onMode = 3;

    // Read in particles to match.
    vector<int> idToMatch;
    int idRead;
    for ( ; ; ) {
      getWord >> idRead;
      if (!getWord) break;
      idToMatch.push_back(abs(idRead));
    }
    int nToMatch = idToMatch.size();

    // Loop over all decay channels.
    for (int i = 0; i < pdt[idTmp].sizeChannels(); ++i) {
      int multi = pdt[idTmp].channel(i).multiplicity();

      // Look for any match at all.
      if (matchKind == 1) {
        bool foundMatch = false;
        for (int j = 0; j < multi; ++j) {
          int idNow =  abs(pdt[idTmp].channel(i).product(j));
          for (int k = 0; k < nToMatch; ++k)
          if (idNow == idToMatch[k]) {foundMatch = true; break;}
          if (foundMatch) break;
        }
        if (foundMatch) pdt[idTmp].channel(i).onMode(onMode);

      // Look for match of all products provided.
      } else {
        int nUnmatched = nToMatch;
        if (multi < nToMatch);
        else if (multi > nToMatch && matchKind == 3);
        else {
          vector<int> idUnmatched;
          for (int k = 0; k < nToMatch; ++k)
            idUnmatched.push_back(idToMatch[k]);
          for (int j = 0; j < multi; ++j) {
            int idNow =  abs(pdt[idTmp].channel(i).product(j));
            for (int k = 0; k < nUnmatched; ++k)
            if (idNow == idUnmatched[k]) {
              idUnmatched[k] = idUnmatched[--nUnmatched];
              break;
            }
            if (nUnmatched == 0) break;
          }
        }
        if (nUnmatched == 0) pdt[idTmp].channel(i).onMode(onMode);
      }
    }
    return true;
  }

  // Rescale all branching ratios by common factor.
  if (property == "rescalebr") {
    double factor;
    getWord >> factor;
    pdt[idTmp].rescaleBR(factor);
    return true;
  }

  // Reset decay table in preparation for new input.
  if (property == "onechannel") pdt[idTmp].clearChannels();

  // Add or change a decay channel: get channel number and new property.
  if (property == "addchannel" || property == "onechannel"
    || isdigit(property[0])) {
    int channel;
    if (property == "addchannel" || property == "onechannel") {
      pdt[idTmp].addChannel();
      channel = pdt[idTmp].sizeChannels() - 1;
      property = "all";
    } else{
      istringstream getChannel(property);
      getChannel >> channel;
      getWord >> property;
      property = toLower(property);
    }

    // Check that channel exists.
    if (channel < 0 || channel >= pdt[idTmp].sizeChannels()) return false;

    // Find decay channel property and value, case by case.
    // At same time also do case where all should be replaced.
    if (property == "onmode" || property == "all") {
      int    onMode = 0;
      string onModeIn;
      getWord >> onModeIn;
      // For onMode allow the optional possibility of Bool input.
      if (isdigit(onModeIn[0])) {
        istringstream getOnMode(onModeIn);
        getOnMode >> onMode;
      } else onMode = (boolString(onModeIn)) ? 1 : 0;
      pdt[idTmp].channel(channel).onMode(onMode);
      if (property == "onmode") return true;
    }
    if (property == "bratio" || property == "all") {
      double bRatio;
      getWord >> bRatio;
      pdt[idTmp].channel(channel).bRatio(bRatio);
      if (property == "bratio") return true;
    }
    if (property == "memode" || property == "all") {
      int meMode;
      getWord >> meMode;
      pdt[idTmp].channel(channel).meMode(meMode);
      if (property == "memode") return true;
    }

    // Scan for products until end of line.
    if (property == "products" || property == "all") {
      int nProd = 0;
      for (int iProd = 0; iProd < 8; ++iProd) {
        int idProd;
        getWord >> idProd;
        if (!getWord) break;
        pdt[idTmp].channel(channel).product(iProd, idProd);
        ++nProd;
      }
      for (int iProd = nProd; iProd < 8; ++iProd)
        pdt[idTmp].channel(channel).product(iProd, 0);
      pdt[idTmp].channel(channel).multiplicity(nProd);
      return true;
    }

    // Rescale an existing branching ratio.
    if (property == "rescalebr") {
      double factor;
      getWord >> factor;
      pdt[idTmp].channel(channel).rescaleBR(factor);
      return true;
    }
  }

  // Return false if failed to recognize property.
  if (warn) os << "\n PYTHIA Error: input property not found in Particle"
    << " Data Table:\n   " << lineIn << "\n";
  readingFailedSave = true;
  return false;

}

//--------------------------------------------------------------------------

// Print out complete or changed table of database in numerical order.

void ParticleData::list(bool changedOnly, bool changedRes, ostream& os) {

  // Table header; output for bool as off/on.
  if (!changedOnly) {
    os << "\n --------  PYTHIA Particle Data Table (complete)  --------"
       << "------------------------------------------------------------"
       << "--------------\n \n";

  } else {
    os << "\n --------  PYTHIA Particle Data Table (changed only)  ----"
       << "------------------------------------------------------------"
       << "--------------\n \n";
  }
  os << "      id   name            antiName         spn chg col      m0"
     << "        mWidth      mMin       mMax       tau0    res dec ext "
     << "vis wid\n             no onMode   bRatio   meMode     products \n";

  // Iterate through the particle data table. Option to skip unchanged.
  int nList = 0;
  for (map<int, ParticleDataEntry>::iterator pdtEntry
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    particlePtr = &pdtEntry->second;
    if ( !changedOnly || particlePtr->hasChanged() ||
      ( changedRes && particlePtr->getResonancePtr() != 0 ) ) {

      // Pick format for mass and width based on mass value.
      double m0Now = particlePtr->m0();
      if (m0Now == 0 || (m0Now > 0.1 && m0Now < 1000.))
        os << fixed << setprecision(5);
      else os << scientific << setprecision(3);

      // Print particle properties.
      ++nList;
      os << "\n" << setw(8) << particlePtr->id() << "  " << left;
      if (particlePtr->name(-1) == "void")
        os << setw(33) << particlePtr->name() << "  ";
      else os << setw(16) << particlePtr->name() << " "
         << setw(16) << particlePtr->name(-1) << "  ";
      os << right << setw(2) << particlePtr->spinType() << "  "
         << setw(2) << particlePtr->chargeType() << "  "
         << setw(2) << particlePtr->colType() << " "
         << setw(10) << particlePtr->m0() << " "
         << setw(10) << particlePtr->mWidth() << " "
         << setw(10) << particlePtr->mMin() << " "
         << setw(10) << particlePtr->mMax() << " "
         << scientific << setprecision(5)
         << setw(12) << particlePtr->tau0() << "  " << setw(2)
         << particlePtr->isResonance() << "  " << setw(2)
         << (particlePtr->mayDecay() && particlePtr->canDecay())
         << "  " << setw(2) << particlePtr->doExternalDecay() << "  "
         << setw(2) << particlePtr->isVisible()<< "  "
         << setw(2) << particlePtr->doForceWidth() << "\n";

      // Loop through the decay channel table for each particle.
      if (particlePtr->sizeChannels() > 0) {
        for (int i = 0; i < int(particlePtr->sizeChannels()); ++i) {
          const DecayChannel& channel = particlePtr->channel(i);
          os << "          "  << setprecision(7)
             << setw(5) << i
             << setw(6) << channel.onMode()
             << fixed<< setw(12) << channel.bRatio()
             << setw(5) << channel.meMode() << " ";
          for (int j = 0; j < channel.multiplicity(); ++j)
            os << setw(8) << channel.product(j) << " ";
          os << "\n";
        }
      }
    }

  }

  // End of loop over database contents.
  if (changedOnly && nList == 0) os << "\n no particle data has been "
    << "changed from its default value \n";
  os << "\n --------  End PYTHIA Particle Data Table  -----------------"
     << "--------------------------------------------------------------"
     << "----------\n" << endl;

}

//--------------------------------------------------------------------------

// Print out partial table of database in input order.

void ParticleData::list(vector<int> idList, ostream& os) {

  // Table header; output for bool as off/on.
  os << "\n --------  PYTHIA Particle Data Table (partial)  ---------"
     << "------------------------------------------------------------"
     << "--------------\n \n";
  os << "      id   name            antiName         spn chg col      m0"
     << "        mWidth      mMin       mMax       tau0    res dec ext "
     << "vis wid\n             no onMode   bRatio   meMode     products \n";

  // Iterate through the given list of input particles.
  for (int i = 0; i < int(idList.size()); ++i) {
    particlePtr = particleDataEntryPtr(idList[i]);

    // Pick format for mass and width based on mass value.
    double m0Now = particlePtr->m0();
    if (m0Now == 0 || (m0Now > 0.1 && m0Now < 1000.))
      os << fixed << setprecision(5);
    else os << scientific << setprecision(3);

    // Print particle properties.
    os << "\n" << setw(8) << particlePtr->id() << "  " << left;
    if (particlePtr->name(-1) == "void")
      os << setw(33) << particlePtr->name() << "  ";
    else os << setw(16) << particlePtr->name() << " "
       << setw(16) << particlePtr->name(-1) << "  ";
    os << right << setw(2) << particlePtr->spinType() << "  "
       << setw(2) << particlePtr->chargeType() << "  "
       << setw(2) << particlePtr->colType() << " "
       << setw(10) << particlePtr->m0() << " "
       << setw(10) << particlePtr->mWidth() << " "
       << setw(10) << particlePtr->mMin() << " "
       << setw(10) << particlePtr->mMax() << " "
       << scientific << setprecision(5)
       << setw(12) << particlePtr->tau0() << "  " << setw(2)
       << particlePtr->isResonance() << "  " << setw(2)
       << (particlePtr->mayDecay() && particlePtr->canDecay())
       << "  " << setw(2) << particlePtr->doExternalDecay() << "  "
       << setw(2) << particlePtr->isVisible() << "  "
       << setw(2) << particlePtr->doForceWidth() << "\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->sizeChannels() > 0) {
      for (int j = 0; j < int(particlePtr->sizeChannels()); ++j) {
        const DecayChannel& channel = particlePtr->channel(j);
        os << "          "  << setprecision(7)
           << setw(5) << j
           << setw(6) << channel.onMode()
           << fixed<< setw(12) << channel.bRatio()
           << setw(5) << channel.meMode() << " ";
        for (int k = 0; k < channel.multiplicity(); ++k)
          os << setw(8) << channel.product(k) << " ";
        os << "\n";
      }
    }

  }

  // End of loop over database contents.
  os << "\n --------  End PYTHIA Particle Data Table  -----------------"
     << "--------------------------------------------------------------"
     << "----------\n" << endl;

}

//--------------------------------------------------------------------------

// Check that table makes sense: e.g. consistent names and mass ranges,
// that branching ratios sum to unity, that charge is conserved and
// that phase space is open in each channel.
// verbosity = 0: mimimal amount of checks, e.g. that no channels open.
//           = 1: further warning if individual channels closed
//                (except for resonances).
//           = 2:  also print branching-ratio-averaged threshold mass.
//      = 11, 12: as 1, 2, but include resonances in detailed checks.

void ParticleData::checkTable(int verbosity, ostream& os) {

  // Header.
  os << "\n --------  PYTHIA Check of Particle Data Table  ------------"
     <<"------\n\n";
  int nErr = 0;

  // Loop through all particles.
  for (map<int, ParticleDataEntry>::iterator pdtEntry
  = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    particlePtr = &pdtEntry->second;

    // Extract some particle properties. Set some flags.
    int    idNow          = particlePtr->id();
    bool   hasAntiNow     = particlePtr->hasAnti();
    int    spinTypeNow    = particlePtr->spinType();
    int    chargeTypeNow  = particlePtr->chargeType();
    int    baryonTypeNow  = particlePtr->baryonNumberType();
    double m0Now          = particlePtr->m0();
    double mMinNow        = particlePtr->mMin();
    double mMaxNow        = particlePtr->mMax();
    double mWidthNow      = particlePtr->mWidth();
    double tau0Now        = particlePtr->tau0();
    bool   isResonanceNow = particlePtr->isResonance();
    bool   hasPrinted     = false;
    bool   studyCloser    = verbosity > 10 || !isResonanceNow;

    // Check that particle name consistent with charge information.
    string particleName = particlePtr->name(1);
    if (particleName.size() > 16) {
      os << " Warning: particle " << idNow << " has name " << particleName
         << " of length " << particleName.size() << "\n";
      hasPrinted = true;
      ++nErr;
    }
    int nPos = 0;
    int nNeg = 0;
    for (int i = 0; i < int(particleName.size()); ++i) {
      if (particleName[i] == '+') ++nPos;
      if (particleName[i] == '-') ++nNeg;
    }
    if ( (nPos > 0 && nNeg > 0) || ( nPos + nNeg > 0
      && 3 * (nPos - nNeg) != chargeTypeNow )) {
      os << " Warning: particle " << idNow << " has name " << particleName
         << " inconsistent with charge type " << chargeTypeNow << "\n";
      hasPrinted = true;
      ++nErr;
    }

    // Check that antiparticle name consistent with charge information.
    if (hasAntiNow) {
      particleName = particlePtr->name(-1);
      if (particleName.size() > 16) {
        os << " Warning: particle " << idNow << " has name " << particleName
           << " of length " << particleName.size() << "\n";
        hasPrinted = true;
        ++nErr;
      }
      nPos = 0;
      nNeg = 0;
      for (int i = 0; i < int(particleName.size()); ++i) {
        if (particleName[i] == '+') ++nPos;
        if (particleName[i] == '-') ++nNeg;
      }
      if ( (nPos > 0 && nNeg > 0) || ( nPos + nNeg > 0
        && 3 * (nPos - nNeg) != -chargeTypeNow )) {
        os << " Warning: particle " << -idNow << " has name "
           << particleName << " inconsistent with charge type "
           << -chargeTypeNow << "\n";
        hasPrinted = true;
        ++nErr;
      }
    }

    // Check that mass, mass range and width are consistent.
    if (particlePtr->useBreitWigner()) {
      if (mMinNow > m0Now) {
        os << " Error: particle " << idNow << " has mMin "
           << fixed << setprecision(5) << mMinNow
           << " larger than m0 " << m0Now << "\n";
        hasPrinted = true;
        ++nErr;
      }
      if (mMaxNow > mMinNow && mMaxNow < m0Now) {
        os << " Error: particle " << idNow << " has mMax "
           << fixed << setprecision(5) << mMaxNow
           << " smaller than m0 " << m0Now << "\n";
        hasPrinted = true;
        ++nErr;
      }
      if (mMaxNow > mMinNow && mMaxNow - mMinNow < mWidthNow) {
        os << " Warning: particle " << idNow << " has mMax - mMin "
           << fixed << setprecision(5) << mMaxNow - mMinNow
           << " smaller than mWidth " << mWidthNow << "\n";
        hasPrinted = true;
        ++nErr;
      }
    }

    // Check that particle does not both have width and lifetime.
    if (mWidthNow > 0. && tau0Now > 0.) {
      os << " Warning: particle " << idNow << " has both nonvanishing width "
         << scientific << setprecision(5) << mWidthNow << " and lifetime "
         << tau0Now << "\n";
      hasPrinted = true;
      ++nErr;
    }

    // Begin study decay channels.
    if (particlePtr->sizeChannels() > 0) {

      // Loop through all decay channels.
      double bRatioSum = 0.;
      double bRatioPos = 0.;
      double bRatioNeg = 0.;
      bool hasPosBR = false;
      bool hasNegBR = false;
      double threshMass = 0.;
      bool openChannel = false;
      for (int i = 0; i < particlePtr->sizeChannels(); ++i) {

        // Extract channel properties.
        int onMode = particlePtr->channel(i).onMode();
        double bRatio = particlePtr->channel(i).bRatio();
        int meMode = particlePtr->channel(i).meMode();
        int mult = particlePtr->channel(i).multiplicity();
        int prod[8];
        for (int j = 0; j < 8; ++j)
          prod[j] = particlePtr->channel(i).product(j);

        // Sum branching ratios. Include off-channels.
        if (onMode == 0 || onMode == 1) bRatioSum += bRatio;
        else if (onMode == 2) {bRatioPos += bRatio; hasPosBR = true;}
        else if (onMode == 3) {bRatioNeg += bRatio; hasNegBR = true;}

        // Error printout when unknown decay product code.
        for (int j = 0; j < 8; ++j) {
          if ( prod[j] != 0 && !isParticle(prod[j]) ) {
            os << " Error: unknown decay product for " << idNow
               << " -> " << prod[j] << "\n";
            hasPrinted = true;
            ++nErr;
            continue;
          }
        }

        // Error printout when no decay products or 0 interspersed.
        int nLast = 0;
        for (int j = 0; j < 8; ++j)
          if (prod[j] != 0) nLast = j + 1;
        if (mult == 0 || mult != nLast) {
          os << " Error: corrupted decay product list for "
             <<  particlePtr->id() << " -> ";
          for (int j = 0; j < 8; ++j) os << prod[j] << " ";
          os << "\n";
          hasPrinted = true;
          ++nErr;
          continue;
        }

        // Check charge conservation and open phase space in decay channel.
        int chargeTypeSum = -chargeTypeNow;
        int baryonTypeSum = -baryonTypeNow;
        double avgFinalMass = 0.;
        double minFinalMass = 0.;
        bool canHandle = true;
        for (int j = 0; j < mult; ++j) {
          chargeTypeSum += chargeType( prod[j] );
          baryonTypeSum += baryonNumberType( prod[j] );
          avgFinalMass += m0( prod[j] );
          minFinalMass += m0Min( prod[j] );
          if (prod[j] == 81 || prod[j] == 82 || prod[j] == 83)
            canHandle = false;
        }
        threshMass += bRatio * avgFinalMass;

        // Error printout when charge or baryon number not conserved.
        if (chargeTypeSum != 0 && canHandle) {
          os << " Error: 3*charge changed by " << chargeTypeSum
             << " in " << idNow << " -> ";
          for (int j = 0; j < mult; ++j) os << prod[j] << " ";
          os << "\n";
          hasPrinted = true;
          ++nErr;
          continue;
        }
        if ( baryonTypeSum != 0 && canHandle && particlePtr->isHadron() ) {
          os << " Error: 3*baryon number changed by " << baryonTypeSum
             << " in " << idNow << " -> ";
          for (int j = 0; j < mult; ++j) os << prod[j] << " ";
          os << "\n";
          hasPrinted = true;
          ++nErr;
          continue;
        }

        // Begin check that some matrix elements are used correctly.
        bool correctME = true;

        // Check matrix element mode 0: recommended not into partons.
        if (meMode == 0 && !isResonanceNow) {
          bool hasPartons = false;
          for (int j = 0; j < mult; ++j) {
            int idAbs = abs(prod[j]);
            if ( idAbs < 10 || idAbs == 21 || idAbs == 81 || idAbs == 82
            || idAbs == 83 || (idAbs > 1000 && idAbs < 10000
            && (idAbs/10)%10 == 0) ) hasPartons = true;
          }
          if (hasPartons) correctME = false;
        }

        // Check matrix element mode 1: rho/omega -> pi+ pi- pi0.
        bool useME1 = ( mult == 3 && spinTypeNow == 3 && idNow > 100
          && idNow < 1000
          && particlePtr->channel(i).contains(211, -211, 111) );
        if ( meMode == 1 && !useME1 ) correctME = false;
        if ( meMode != 1 &&  useME1 ) correctME = false;

        // Check matrix element mode 2: polarization in V -> PS + PS.
        bool useME2 = ( mult == 2 && spinTypeNow == 3  && idNow > 100
          && idNow < 1000 && spinType(prod[0]) == 1
          && spinType(prod[1]) == 1 );
        if ( meMode == 2 && !useME2 ) correctME = false;
        if ( meMode != 2 &&  useME2 ) correctME = false;

        // Check matrix element mode 11, 12 and 13: Dalitz decay with
        // one or more particles in addition to the lepton pair,
        // or double Dalitz decay.
        bool useME11 = ( mult == 3 && !isResonanceNow
          && (prod[1] == 11 || prod[1] == 13 || prod[1] == 15)
          && prod[2] == -prod[1] );
        bool useME12 = ( mult > 3 && !isResonanceNow
          && (prod[mult - 2] == 11 || prod[mult - 2] == 13
          || prod[mult - 2] == 15) && prod[mult - 1] == -prod[mult - 2] );
        bool useME13 = ( mult == 4 && !isResonanceNow
          && (prod[0] == 11 || prod[0] == 13) && prod[1] == -prod[0]
          && (prod[2] == 11 || prod[2] == 13) && prod[3] == -prod[2] );
        if (useME13) useME12 = false;
        if ( meMode == 11 && !useME11 ) correctME = false;
        if ( meMode != 11 &&  useME11 ) correctME = false;
        if ( meMode == 12 && !useME12 ) correctME = false;
        if ( meMode != 12 &&  useME12 ) correctME = false;
        if ( meMode == 13 && !useME13 ) correctME = false;
        if ( meMode != 13 &&  useME13 ) correctME = false;

        // Check matrix element mode 21: tau -> nu_tau hadrons.
        bool useME21 = (idNow == 15 && mult > 2 && prod[0] == 16
          && abs(prod[1]) > 100);
        if ( meMode == 21 && !useME21 ) correctME = false;
        if ( meMode != 21 &&  useME21 ) correctME = false;

        // Check matrix element mode 22, but only for semileptonic decay.
        // For a -> b c d require types u = 2, ubar = -2, d = 1, dbar = -1.
        if ( isLepton(prod[0]) && isLepton(prod[1]) ) {
          bool useME22 = false;
          int typeA = 0;
          int typeB = 0;
          int typeC = 0;
          if (particlePtr->isLepton()) {
            typeA = (idNow > 0) ? 1 + (idNow-1)%2 : -1 - (1-idNow)%2;
          } else if (particlePtr->isHadron()) {
            int hQ = particlePtr->heaviestQuark();
            // Special case: for B_c either bbar or c decays.
            if (idNow == 541 && heaviestQuark(prod[2]) == -5) hQ = 4;
            typeA = (hQ > 0) ? 1 + (hQ-1)%2 : -1 - (1-hQ)%2;
          }
          typeB = (prod[0] > 0) ? 1 + (prod[0]-1)%2 : -1 - (1-prod[0])%2;
          typeC = (prod[1] > 0) ? 1 + (prod[1]-1)%2 : -1 - (1-prod[1])%2;
          // Special cases.
          if ( (idNow == 130 || idNow == 310) && typeC * typeA < 0)
            typeA = -typeA;
          if (mult == 3 && idNow == 2112 && prod[2] == 2212)
            typeA = 3 - typeA;
          if (mult == 3 && idNow == 3222 && prod[2] == 3122)
            typeA = 3 - typeA;
          if (mult > 2 && typeC == typeA && typeB * typeC < 0
            && (typeB + typeC)%2 != 0) useME22 = true;
          if ( meMode == 22 && !useME22 ) correctME = false;
          if ( meMode != 22 &&  useME22 ) correctME = false;
        }

        // Check for matrix element mode 31.
        if (meMode == 31) {
          int nGamma = 0;
          for (int j = 0; j < mult; ++j) if (prod[j] == 22) ++nGamma;
          if (nGamma != 1) correctME = false;
        }

        // Check for unknown mode, or resonance-only mode.
        if ( !isResonanceNow && (meMode < 0 || (meMode > 2 && meMode <= 10)
          || (meMode > 13 && meMode <= 20) || (meMode > 23 && meMode <= 30)
          || (meMode > 31 && meMode <= 41) || meMode == 51 || meMode == 61
          || meMode == 71 || (meMode > 80 && meMode <= 90)
          || (!particlePtr->isOctetHadron() && meMode > 92) ) )
          correctME = false;

        // Print if incorrect matrix element mode.
        if ( !correctME ) {
          os << " Warning: meMode " << meMode << " used for "
             << idNow << " -> ";
          for (int j = 0; j < mult; ++j) os << prod[j] << " ";
          os << "\n";
          hasPrinted = true;
          ++nErr;
        }

        // Warning printout when no phase space for decay.
        if ( studyCloser && verbosity > 0  && canHandle && onMode > 0
          && particlePtr->m0Min() - minFinalMass < 0. ) {
          if (particlePtr->m0Max() - minFinalMass < 0.)
            os << " Error: decay never possible for ";
          else  os << " Warning: decay sometimes not possible for ";
          os << idNow << " -> ";
          for (int j = 0; j < mult; ++j) os << prod[j] << " ";
          os << "\n";
          hasPrinted = true;
          ++nErr;
        }

        // End loop over decay channels.
        if (onMode > 0 && bRatio > 0.) openChannel = true;
      }

      // Optional printout of threshold.
      if (verbosity%10 > 1 && particlePtr->useBreitWigner()) {
        threshMass /= bRatioSum;
        os << " Info: particle " << idNow << fixed << setprecision(5)
           << " has average mass threshold " << threshMass
           << " while mMin is " << mMinNow << "\n";
        hasPrinted = true;
      }

      // Error printout when no acceptable decay channels found.
      if (studyCloser && !openChannel) {
        os << " Error: no acceptable decay channel found for particle "
           << idNow << "\n";
        hasPrinted = true;
        ++nErr;
      }

      // Warning printout when branching ratios do not sum to unity.
      if (studyCloser && (!hasAntiNow || (!hasPosBR && !hasNegBR))
        && abs(bRatioSum + bRatioPos - 1.) > 1e-8) {
        os << " Warning: particle " << idNow  << fixed << setprecision(8)
           << " has branching ratio sum " << bRatioSum << "\n";
        hasPrinted = true;
        ++nErr;
      } else if (studyCloser && hasAntiNow
        && (abs(bRatioSum + bRatioPos - 1.) > 1e-8
        || abs(bRatioSum + bRatioNeg - 1.) > 1e-8)) {
        os << " Warning: particle " << idNow  << fixed << setprecision(8)
           << " has branching ratio sum " << bRatioSum + bRatioPos
           << " while its antiparticle has " << bRatioSum + bRatioNeg
           << "\n";
        hasPrinted = true;
        ++nErr;
      }

    // End study of decay channels and loop over particles.
    }
    if (hasPrinted) os << "\n";
  }

  // Final output. Done.
  os << " Total number of errors and warnings is " << nErr << "\n";
  os << "\n --------  End PYTHIA Check of Particle Data Table  --------"
     << "------\n" << endl;

}

//--------------------------------------------------------------------------

// Return the id of the sequentially next particle stored in table.

int ParticleData::nextId(int idIn) {

  // Return 0 for negative or unknown codes. Return first for 0.
  if (idIn < 0 || (idIn > 0 && !isParticle(idIn))) return 0;
  if (idIn == 0) return pdt.begin()->first;

  // Find pointer to current particle and step up. Return 0 if impossible.
  map<int, ParticleDataEntry>::const_iterator pdtIn = pdt.find(idIn);
  if (pdtIn == pdt.end()) return 0;
  ++pdtIn;
  if (pdtIn == pdt.end()) return 0;
  return pdtIn->first;

}

//--------------------------------------------------------------------------

// Fractional width associated with open channels of one or two resonances.

double ParticleData::resOpenFrac(int id1In, int id2In, int id3In) {

  // Default value.
  double answer = 1.;

  // First resonance.
  if (isParticle(id1In)) answer  = pdt[abs(id1In)].resOpenFrac(id1In);

  // Possibly second resonance.
  if (isParticle(id2In)) answer *= pdt[abs(id2In)].resOpenFrac(id2In);

  // Possibly third resonance.
  if (isParticle(id3In)) answer *= pdt[abs(id2In)].resOpenFrac(id3In);

  // Done.
  return answer;

}

//--------------------------------------------------------------------------

// Extract XML value string following XML attribute.

string ParticleData::attributeValue(string line, string attribute) {

  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

}

//--------------------------------------------------------------------------

// Extract XML bool value following XML attribute.

bool ParticleData::boolAttributeValue(string line, string attribute) {

  string valString = attributeValue(line, attribute);
  if (valString == "") return false;
  return boolString(valString);

}

//--------------------------------------------------------------------------

// Extract XML int value following XML attribute.

int ParticleData::intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;

}

//--------------------------------------------------------------------------

// Extract XML double value following XML attribute.

double ParticleData::doubleAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0.;
  istringstream valStream(valString);
  double doubleVal;
  valStream >> doubleVal;
  return doubleVal;

}

//==========================================================================

} // end namespace Pythia8
