// HadronScatter.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef Pythia8_HadronScatter_H
#define Pythia8_HadronScatter_H

#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/PythiaComplex.h"

namespace Pythia8 {

class SigmaPartialWave {
public:
  // Initialisation
  bool init(int, string, string, Info *, ParticleData *, Rndm *);

  // Read data file
  bool readFile(string, string);

  // Set the subprocess/incoming particles
  bool setSubprocess(int);
  bool setSubprocess(int, int);

  // Return sigma total/elastic, dSigma/dCos(theta)
  double sigmaEl(double Wcm)  { return sigma(0, Wcm); }
  double sigmaTot(double Wcm) { return sigma(1, Wcm); }
  double dSigma(double Wcm, double cTheta) { return sigma(2, Wcm, cTheta); }

  // Return a cos(theta) value
  double pickCosTheta(double);

  // Return maximum sigma elastic
  double getSigmaElMax() { return sigElMax; }

private:
  // Pointers
  Info         *infoPtr;
  ParticleData *particleDataPtr;
  Rndm         *rndmPtr;

  // Constants
  static const int     LSHIFT, ISHIFT, SUBBIN, ITER;
  static const double  CONVERT2MB, WCMBIN, CTBIN, MASSSAFETY, GRIDSAFETY;
  static const complex BINEND;

  // Settings
  int process, subprocess, subprocessMax, norm;

  // Current incoming and maximum L/I values
  int idA, idB, Lmax, Imax;

  // Masses of incoming, last bin, maximum sigma elastic
  double mA, mB, binMax, sigElMax;

  // Map subprocess to incoming and vice versa:
  //   sp2in[subprocess] = idA, idB
  //   in2sp[idA, idB]   = subprocess
  map < int, pair < int, int > > sp2in;
  map < pair < int, int >, int > in2sp;

  // Isospin coefficients isoCoeff[subprocess][2I]
  map < int, map < int, double > > isoCoeff;

  // Storage for partial wave data:
  //   pwData[L * LSHIFT + 2I * ISHIFT + 2J][eNow] = T
  map < int, map < double, complex > > pwData;

  // Values of Pl and Pl' as computed by legendreP
  vector < double > PlVec, PlpVec;

  // Integration grid [subprocess][WcmBin][cosThetaBin]
  vector < vector < vector < double > > > gridMax;
  vector < vector < double > >            gridNorm;

  // Setup subprocesses (including isospin coefficients)
  void setupSubprocesses();

  // Setup grids for integration
  void setupGrid();

  // Routine for calculating sigma total/elastic and dSigma/dCos(theta)
  double sigma(int, double, double = 0.);

  // Generate Legendre polynomials (and optionally derivatives)
  void legendreP(double, bool = false);

};


//==========================================================================

// HadronScatterPair class
//   Simple class to hold details of a pair of hadrons which will scatter.
//   Stores indices in event record and the measure used for ordering

// Store a pair of indices
typedef pair < int, int > HSIndex;

class HadronScatterPair {
public:
  // Constructor
  HadronScatterPair(const HSIndex &i1in, int yt1in, int pt1in,
                    const HSIndex &i2in, int yt2in, int pt2in,
                    double measureIn) :
      i1(i1in), yt1(yt1in), pt1(pt1in),
      i2(i2in), yt2(yt2in), pt2(pt2in),
      measure(measureIn) {}

  // Operator for sorting according to ordering measure
  bool operator<(const HadronScatterPair& in) const {
    return this->measure < in.measure;
  }

  // Indices into event record of hadrons to scatter
  HSIndex i1;
  int     yt1, pt1;
  HSIndex i2;
  int     yt2, pt2;
  // Ordering measure
  double measure;
};


//==========================================================================

// HadronScatter class

class HadronScatter {

public:

  // Constructor.
  HadronScatter() {}

  // Initialisation
  bool init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
            ParticleData *particleDataPtr);

  // Perform all hadron scatterings
  void scatter(Event&);

private:

  // Pointer to various information on the generation.
  Info* infoPtr;
  Rndm* rndmPtr;

  // Settings
  bool   doHadronScatter, afterDecay, allowDecayProd,
         scatterRepeat, doTile;
  int    hadronSelect, scatterProb;
  double Npar, kPar, pPar, jPar, rMax, rMax2;
  double pTsigma, pTsigma2, pT0MPI;

  // Tiling
  int    ytMax, ptMax;
  double yMin, yMax, ytSize, ptSize;
  vector < vector < set < HSIndex > > > tile;

  // Keep track of scattered pairs
  set < HSIndex > scattered;

  // Partial wave amplitudes
  SigmaPartialWave sigmaPW[3];

  // Maximum sigma elastic
  double sigElMax;

  // Decide if a hadron can scatter
  bool canScatter(Event &, int);

  // Probability for a pair of hadrons to scatter
  bool doesScatter(Event &, const HSIndex &, const HSIndex &);

  // Calculate measure for ordering of scatterings
  double measure(Event &, int, int);

  // Perform a single hadron scattering
  bool hadronScatter(Event &, HadronScatterPair &);

  // Tiling functions
  bool tileIntProb(vector < HadronScatterPair > &, Event &,
                   const HSIndex &, int, int, bool);
  int yTile(Event& event, int idx) {
    return (doTile) ? int((event[idx].y() - yMin) / ytSize) : 0;
  }
  int pTile(Event& event, int idx) {
    return (doTile) ? int((event[idx].phi() + M_PI) / ptSize) : 0;
  }

  // Debug
  void debugOutput();
};

//==========================================================================


} // end namespace Pythia8

#endif // Pythia8_HadronScatter_H

