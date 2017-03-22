// FragmentationFlavZpT.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
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
  StringFlav() {}

  // Destructor.
  virtual ~StringFlav() {}

  // Initialize data members.
  virtual void init(Settings& settings, Rndm* rndmPtrIn);

  // Pick a light d, u or s quark according to fixed ratios.
  int pickLightQ() { double rndmFlav = probQandS * rndmPtr->flat();
    if (rndmFlav < 1.) return 1; if (rndmFlav < 2.) return 2; return 3; }

  // Pick a new flavour (including diquarks) given an incoming one.
  virtual FlavContainer pick(FlavContainer& flavOld);

  // Combine two flavours (including diquarks) to produce a hadron.
  virtual int combine(FlavContainer& flav1, FlavContainer& flav2);

  // Ditto, simplified input argument for simple configurations.
  virtual int combine( int id1, int id2, bool keepTrying = true) {
    FlavContainer flag1(id1); FlavContainer flag2(id2);
    for (int i = 0; i < 100; ++i) { int idNew = combine( flag1, flag2);
      if (idNew != 0 || !keepTrying) return idNew;} return 0;}

  // Assign popcorn quark inside an original (= rank 0) diquark.
  void assignPopQ(FlavContainer& flav);

  // Combine two quarks to produce a diquark.
  int makeDiquark(int id1, int id2, int idHad = 0);

protected:

  // Pointer to the random number generator.
  Rndm*  rndmPtr;

private:

  // Constants: could only be changed in the code itself.
  static const int    mesonMultipletCode[6];
  static const double baryonCGOct[6], baryonCGDec[6];

  // Initialization data, to be read from Settings.
  bool   suppressLeadingB;
  double probQQtoQ, probStoUD, probSQtoQQ, probQQ1toQQ0, probQandQQ,
         probQandS, probQandSinQQ, probQQ1corr, probQQ1corrInv, probQQ1norm,
         probQQ1join[4], mesonRate[4][6], mesonRateSum[4], mesonMix1[2][6],
         mesonMix2[2][6], etaSup, etaPrimeSup, decupletSup, baryonCGSum[6],
         baryonCGMax[6], popcornRate, popcornSpair, popcornSmeson, scbBM[3],
         popFrac, popS[3], dWT[3][7], lightLeadingBSup, heavyLeadingBSup;

};

//==========================================================================

// The StringZ class is used to sample the fragmentation function f(z).

class StringZ {

public:

  // Constructor.
  StringZ() {}

  // Destructor.
  virtual ~StringZ() {}

  // Initialize data members.
  virtual void init(Settings& settings, ParticleData& particleData,
    Rndm* rndmPtrIn);

  // Fragmentation function: top-level to determine parameters.
  virtual double zFrag( int idOld, int idNew = 0, double mT2 = 1.);

  // Parameters for stopping in the middle; overloaded for Hidden Valley.
  virtual double stopMass() {return stopM;}
  virtual double stopNewFlav() {return stopNF;}
  virtual double stopSmear() {return stopS;}

  // b fragmentation parameter needed to weight final two solutions.
  virtual double bAreaLund() {return bLund;}

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

};

//==========================================================================

// The StringPT class is used to select select transverse momenta.

class StringPT {

public:

  // Constructor.
  StringPT() {}

  // Destructor.
  virtual ~StringPT() {}

  // Initialize data members.
  virtual void init(Settings& settings, ParticleData& particleData,
    Rndm* rndmPtrIn);

  // Return px and py as a pair in the same call.
  pair<double, double>  pxy();

  // Gaussian suppression of given pT2; used in MiniStringFragmentation.
  double suppressPT2(double pT2) { return exp( -pT2 / sigma2Had); }

protected:

  // Constants: could only be changed in the code itself.
  static const double SIGMAMIN;

  // Initialization data, to be read from Settings.
  double sigmaQ, enhancedFraction, enhancedWidth, sigma2Had;

  // Pointer to the random number generator.
  Rndm*  rndmPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_FragmentationFlavZpT_H
