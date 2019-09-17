// FragmentationFlavZpT.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains helper classes for fragmentation.
// StringFlav is used to select quark and hadron flavours.
// StringPT is used to select transverse momenta.
// StringZ is used to sample the fragmentation function f(z).

#ifndef Pythia8_FragmentationFlavZpT_H
#define Pythia8_FragmentationFlavZpT_H

#include "Pythia8/Basics.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The FlavContainer class is a simple container for flavour,
// including the extra properties needed for popcorn baryon handling.
// id = current flavour.
// rank = current rank; 0 for endpoint flavour and then increase by 1.
// nPop = number of popcorn mesons yet to be produced (1 or 0).
// idPop = (absolute sign of) popcorn quark, shared between B and Bbar.
// idVtx = (absolute sign of) vertex (= non-shared) quark in diquark.

class FlavContainer {

public:

  // Constructor.
  FlavContainer(int idIn = 0, int rankIn = 0, int nPopIn = 0,
    int idPopIn = 0, int idVtxIn = 0) : id(idIn), rank(rankIn),
    nPop(nPopIn), idPop(idPopIn), idVtx(idVtxIn) {}

  // Overloaded equal operator.
  FlavContainer& operator=(const FlavContainer& flav) { if (this != &flav) {
    id = flav.id; rank = flav.rank; nPop = flav.nPop; idPop = flav.idPop;
    idVtx = flav.idVtx; } return *this; }

  // Invert flavour.
  FlavContainer& anti() {id = -id; return *this;}

  // Read in a container into another, without/with id sign flip.
  FlavContainer& copy(const FlavContainer& flav) { if (this != &flav) {
    id = flav.id; rank = flav.rank; nPop = flav.nPop; idPop = flav.idPop;
    idVtx = flav.idVtx; } return *this; }
  FlavContainer& anti(const FlavContainer& flav) { if (this != &flav) {
    id = -flav.id; rank = flav.rank; nPop = flav.nPop; idPop = flav.idPop;
    idVtx = flav.idVtx; } return *this; }

  // Check whether is diquark.
  bool isDiquark() {int idAbs = abs(id);
    return (idAbs > 1000 && idAbs < 10000 && (idAbs/10)%10 == 0);}

  // Stored properties.
  int id, rank, nPop, idPop, idVtx;

};

//==========================================================================

// The StringFlav class is used to select quark and hadron flavours.

class StringFlav {

public:

  // Constructor.
  StringFlav() : rndmPtr(), particleDataPtr(), infoPtr(), suppressLeadingB(),
    mT2suppression(), useWidthPre(), probQQtoQ(), probStoUD(), probSQtoQQ(),
    probQQ1toQQ0(), probQandQQ(), probQandS(), probQandSinQQ(), probQQ1corr(),
    probQQ1corrInv(), probQQ1norm(), probQQ1join(), mesonRate(),
    mesonRateSum(), mesonMix1(), mesonMix2(), etaSup(), etaPrimeSup(),
    decupletSup(), baryonCGSum(), baryonCGMax(), popcornRate(), popcornSpair(),
    popcornSmeson(), scbBM(), popFrac(), popS(), dWT(), lightLeadingBSup(),
    heavyLeadingBSup(), sigmaHad(), widthPreStrange(), widthPreDiquark(),
    thermalModel(), mesonNonetL1(), temperature(), tempPreFactor(),
    nNewQuark(), mesMixRate1(), mesMixRate2(), mesMixRate3(),
    baryonOctWeight(), baryonDecWeight(), closePacking(), exponentMPI(),
    exponentNSP(), hadronIDwin(0), idNewWin(0), hadronMassWin(-1.0) {}

  // Destructor.
  virtual ~StringFlav() {}

  // Initialize data members.
  virtual void init(Settings& settings, ParticleData* particleDataPtrIn,
    Rndm* rndmPtrIn, Info* infoPtrIn);

  // Pick a light d, u or s quark according to fixed ratios.
  int pickLightQ() { double rndmFlav = probQandS * rndmPtr->flat();
    if (rndmFlav < 1.) return 1;
    if (rndmFlav < 2.) return 2;
    return 3; }

  // Pick a new flavour (including diquarks) given an incoming one,
  // either by old standard Gaussian or new alternative exponential.
  virtual FlavContainer pick(FlavContainer& flavOld, double pT = -1.0,
    double nNSP = 0.0) { hadronIDwin = 0; idNewWin = 0; hadronMassWin = -1.0;
    if ( (thermalModel || mT2suppression) && (pT >= 0.0) )
      return pickThermal(flavOld, pT, nNSP);
    return pickGauss(flavOld); }
  virtual FlavContainer pickGauss(FlavContainer& flavOld);
  virtual FlavContainer pickThermal(FlavContainer& flavOld,
    double pT, double nNSP);

  // Combine two flavours (including diquarks) to produce a hadron.
  virtual int combine(FlavContainer& flav1, FlavContainer& flav2);

  // Ditto, simplified input argument for simple configurations.
  virtual int combineId( int id1, int id2, bool keepTrying = true) {
    FlavContainer flag1(id1); FlavContainer flag2(id2);
    for (int i = 0; i < 100; ++i) { int idNew = combine( flag1, flag2);
      if (idNew != 0 || !keepTrying) return idNew;} return 0;}

  // Return chosen hadron in case of thermal model.
  virtual int getHadronIDwin() { return hadronIDwin; }

  // Combine two flavours into hadron for last two remaining flavours
  // for thermal model.
  virtual int combineLastThermal(FlavContainer& flav1, FlavContainer& flav2,
    double pT, double nNSP);

  // General function, decides whether to just return the hadron id
  // if thermal model was use or whether to combine the two flavours.
  virtual int getHadronID(FlavContainer& flav1, FlavContainer& flav2,
    double pT = -1.0, double nNSP = 0, bool finalTwo = false) {
    if (finalTwo) return ((thermalModel || mT2suppression) ?
      combineLastThermal(flav1, flav2, pT, nNSP) : combine(flav1, flav2));
    if ((thermalModel || mT2suppression)&& (hadronIDwin != 0)
      && (idNewWin != 0)) return getHadronIDwin();
    return combine(flav1, flav2); }

  // Return hadron mass. Used one if present, pick otherwise.
  virtual double getHadronMassWin(int idHad) { return
    ((hadronMassWin < 0.0) ? particleDataPtr->mSel(idHad) : hadronMassWin); }

  // Assign popcorn quark inside an original (= rank 0) diquark.
  void assignPopQ(FlavContainer& flav);

  // Combine two quarks to produce a diquark.
  int makeDiquark(int id1, int id2, int idHad = 0);

  // Check if quark-diquark combination should be added. If so add.
  void addQuarkDiquark(vector< pair<int,int> >& quarkCombis,
    int qID, int diqID, int hadronID) {
    bool allowed = true;
    for (int iCombi = 0; iCombi < int(quarkCombis.size()); iCombi++)
      if ( (qID   == quarkCombis[iCombi].first ) &&
           (diqID == quarkCombis[iCombi].second) ) allowed = false;
    if (allowed) quarkCombis.push_back( (hadronID > 0) ?
      make_pair( qID,  diqID) : make_pair(-qID, -diqID) ); }

  // Get spin counter for mesons.
  int getMesonSpinCounter(int hadronID) { hadronID = abs(hadronID);
    int j = (hadronID % 10);
    if (hadronID <  1000) return ((j==1) ? 0 : ( (j==3) ? 1 : 5 ));
    if (hadronID < 20000) return ((j==1) ? 3 : 2);
    if (hadronID > 20000) return 4;
    return -1; }

protected:

  // Pointer to the random number generator.
  Rndm*  rndmPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to event information.
  Info* infoPtr;

  // Constants: could only be changed in the code itself.
  static const int    mesonMultipletCode[6];
  static const double baryonCGOct[6], baryonCGDec[6];

  // Settings for Gaussian model.
  bool   suppressLeadingB, mT2suppression, useWidthPre;
  double probQQtoQ, probStoUD, probSQtoQQ, probQQ1toQQ0, probQandQQ,
         probQandS, probQandSinQQ, probQQ1corr, probQQ1corrInv, probQQ1norm,
         probQQ1join[4], mesonRate[4][6], mesonRateSum[4], mesonMix1[2][6],
         mesonMix2[2][6], etaSup, etaPrimeSup, decupletSup, baryonCGSum[6],
         baryonCGMax[6], popcornRate, popcornSpair, popcornSmeson, scbBM[3],
         popFrac, popS[3], dWT[3][7], lightLeadingBSup, heavyLeadingBSup;
  double sigmaHad, widthPreStrange, widthPreDiquark;

  // Settings for thermal model.
  bool   thermalModel, mesonNonetL1;
  double temperature, tempPreFactor;
  int    nNewQuark;
  double mesMixRate1[2][6], mesMixRate2[2][6], mesMixRate3[2][6];
  double baryonOctWeight[6][6][6][2], baryonDecWeight[6][6][6][2];

  // Settings used by both models.
  bool   closePacking;
  double exponentMPI, exponentNSP;

  // Key = hadron id, value = list of constituent ids.
  map< int, vector< pair<int,int> > > hadronConstIDs;
  // Key = initial (di)quark id, value = list of possible hadron ids
  //                                     + nr in hadronConstIDs.
  map< int, vector< pair<int,int> > > possibleHadrons;
  // Key = initial (di)quark id, value = prefactor to multiply rate.
  map< int, vector<double> > possibleRatePrefacs;
  // Similar, but for combining the last two (di)quarks. Key = (di)quark pair.
  map< pair<int,int>, vector< pair<int,int> > > possibleHadronsLast;
  map< pair<int,int>, vector<double> > possibleRatePrefacsLast;

  // Selection in thermal model.
  int    hadronIDwin, idNewWin;
  double hadronMassWin;

};

//==========================================================================

// Auxiliary class to encapsulate (unnormalised) Lund FF.

class LundFFRaw : public FunctionEncapsulator {

public:

  // Constructor and destructor.
  LundFFRaw() {};
  virtual ~LundFFRaw() {};

  // f(z) = (1-z)^a / z^c * Exp ( - b * mT2 / z).
  virtual double f(vector<double> args) {
    if (args.size() < 5) return -1.;
    double z   = args[0];
    if (z <= 0. || z >= 1.) return 0.;
    double a   = args[1];
    double b   = args[2];
    double c   = args[3];
    double mT2 = args[4];
    return pow(1. - z, a) / pow(z, c) * exp(-b * mT2 / z);
  }

};

//==========================================================================

// Auxiliary class to encapsulate average, <z>, of Lund FF.

class LundFFAvg : public FunctionEncapsulator {

public:

  // Constructor and destructor.
  LundFFAvg() : check() {};
  virtual ~LundFFAvg() {};

  // <z> = integral( z * f(z) dz ) / integral( f(z) dz ),
  // argsIn[] : a, b, c, mT2, and optionally integrator tolerance parameter.
  virtual double f(vector<double> argsIn) {

    // Do sanity checks and set initial values.
    if (argsIn.size() < 4) return -1.;
    double tol = 1.e-6;
    if (argsIn.size() >= 5) tol = argsIn[4];
    double denom = 1.;
    double numerator = 0.;

    // Evaluate denominator integral. LundFF expects args: z, a, b, c, mT2.
    vector<double> args(1);
    args.insert(args.end(), argsIn.begin(), argsIn.end());
    check = lundFF.integrateGauss(denom, 0, 0., 1., args, tol);
    if ( !check || denom <= 0. ) return -1.;

    // Evaluate numerator integral. Evaluating z*f is equivalent to c -> c-1.
    args[3] -= 1.;
    check = lundFF.integrateGauss(numerator, 0, 0., 1., args, tol);
    if ( !check || numerator < 0. ) return -1.;
    return numerator/denom;
  }

private:

  LundFFRaw lundFF;
  bool check;
};

//==========================================================================

// The StringZ class is used to sample the fragmentation function f(z).

class StringZ {

public:

  // Constructor.
  StringZ() : useNonStandC(), useNonStandB(), useNonStandH(), usePetersonC(),
    usePetersonB(), usePetersonH(), mc2(), mb2(), aLund(), bLund(),
    aExtraSQuark(), aExtraDiquark(), rFactC(), rFactB(), rFactH(), aNonC(),
    aNonB(), aNonH(), bNonC(), bNonB(), bNonH(), epsilonC(), epsilonB(),
    epsilonH(), stopM(), stopNF(), stopS(), rndmPtr(), infoPtr() {}

  // Destructor.
  virtual ~StringZ() {}

  // Initialize data members.
  virtual void init(Settings& settings, ParticleData& particleData,
    Rndm* rndmPtrIn, Info* infoPtrIn);

  // Fragmentation function: top-level to determine parameters.
  virtual double zFrag( int idOld, int idNew = 0, double mT2 = 1.);

  // Parameters for stopping in the middle; overloaded for Hidden Valley.
  virtual double stopMass() {return stopM;}
  virtual double stopNewFlav() {return stopNF;}
  virtual double stopSmear() {return stopS;}

  // a and b fragmentation parameters needed in some operations.
  virtual double aAreaLund() {return aLund;}
  virtual double bAreaLund() {return bLund;}

  // Method to derive bLund from <z> (for fixed a and reference mT2).
  bool deriveBLund(Settings& settings, ParticleData& particleData);

protected:

  // Constants: could only be changed in the code itself.
  static const double CFROMUNITY, AFROMZERO, AFROMC, EXPMAX;

  // Initialization data, to be read from Settings.
  bool   useNonStandC, useNonStandB, useNonStandH,
         usePetersonC, usePetersonB, usePetersonH;
  double mc2, mb2, aLund, bLund, aExtraSQuark, aExtraDiquark, rFactC,
         rFactB, rFactH, aNonC, aNonB, aNonH, bNonC, bNonB, bNonH,
         epsilonC, epsilonB, epsilonH, stopM, stopNF, stopS;

  // Fragmentation function: select z according to provided parameters.
  double zLund( double a, double b, double c = 1.);
  double zPeterson( double epsilon);

  // Pointer to the random number generator.
  Rndm*  rndmPtr;

  // Pointer to event information.
  Info* infoPtr;

};

//==========================================================================

// The StringPT class is used to select select transverse momenta.

class StringPT {

public:

  // Constructor.
  StringPT() : useWidthPre(), sigmaQ(), enhancedFraction(), enhancedWidth(),
    sigma2Had(), widthPreStrange(), widthPreDiquark(), thermalModel(),
    temperature(), tempPreFactor(), fracSmallX(), closePacking(),
    exponentMPI(), exponentNSP(), particleDataPtr(), rndmPtr(), infoPtr() {}

  // Destructor.
  virtual ~StringPT() {}

  // Initialize data members.
  virtual void init(Settings& settings, ParticleData* particleDataPtr,
    Rndm* rndmPtrIn, Info* infoPtrIn);

  // General function, return px and py as a pair in the same call
  // in either model.
  pair<double, double>  pxy(int idIn, double nNSP = 0.0) {
    return (thermalModel ? pxyThermal(idIn, nNSP) :
    pxyGauss(idIn, nNSP)); }
  pair<double, double>  pxyGauss(int idIn = 0, double nNSP = 0.0);
  pair<double, double>  pxyThermal(int idIn, double nNSP = 0.0);

  // Gaussian suppression of given pT2; used in MiniStringFragmentation.
  double suppressPT2(double pT2) { return (thermalModel ?
    exp(-sqrt(pT2)/temperature) : exp(-pT2/sigma2Had)); }

protected:

  // Constants: could only be changed in the code itself.
  static const double SIGMAMIN;

  // Initialization data, to be read from Settings.
  // Gaussian model.
  bool   useWidthPre;
  double sigmaQ, enhancedFraction, enhancedWidth, sigma2Had,
         widthPreStrange, widthPreDiquark;
  // Thermal model.
  bool   thermalModel;
  double temperature, tempPreFactor, fracSmallX;
  // Both.
  bool   closePacking;
  double exponentMPI, exponentNSP;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*  rndmPtr;

  // Pointer to event information.
  Info* infoPtr;

private:

  // Evaluate Bessel function K_{1/4}(x).
  double BesselK14(double x);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_FragmentationFlavZpT_H
