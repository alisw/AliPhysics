// DeuteronProduction.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef Pythia8_DeuteronProduction_H
#define Pythia8_DeuteronProduction_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The DeuteronProduction class.

class DeuteronProduction {

public:

  // Constructor.
  DeuteronProduction() : valid(true), models(), ids(), parms(), masses(),
    norm(), mPion(), mSafety(), kMin(), kMax(), kTol(), kSteps(),
    infoPtr(), pdbPtr(), rndmPtr() {}

  // Find settings. Precalculate table used to find momentum shifts.
  bool init(Info* infoPtrIn, Settings& settings, ParticleData* pdbPtrIn,
    Rndm* rndmPtrIn);

  // Form deuterons in an event.
  bool combine(Event& event);

private:

  bool valid;                            // Flag if class has been initialized.
  vector<int> models;                    // Cross-section mode per channel.
  vector<vector<int> > ids;              // IDs and charges per channel.
  vector<vector<double> > parms, masses; // Parameters and masses per channel.
  double norm;                           // Overall normalization scale.
  double mPion;                          // Mass of the pion.
  double mSafety;                        // Safety margin for decays.
  double kMin, kMax, kTol;               // Bracketing/tolerance in k for max.
  int kSteps;                            // Number of steps for grid search.

  // Internal pointers.
  Info* infoPtr;        // Information pointer.
  ParticleData* pdbPtr; // Particle database pointer.
  Rndm* rndmPtr;        // Random number generator pointer.

  // Constants: could only be changed in the code itself.
  static const int NTRYDECAY;           // Number of times to try a decay.
  static const double WTCORRECTION[11]; // M-generator parameters.

  // Bind the nucleon-pair combinations.
  void bind(Event& event, vector<int>& prts);

  // Build the nucleon-pair combinations and shuffle.
  void combos(Event& event, vector<int>& prts, vector<pair<int, int> > &cmbs);

  // Single pion final state fit, equations 10/13/14 of arXiv:1504.07242.
  double fit(double k, vector<double>& c, int i);

  // Return the cross-section for a given channel.
  double sigma(double k, int chn);

  // N-body decay using the M-generator algorithm.
  bool decay(Event& event, int idx0, int idx1, int chn);

  // Helper methods.
  void maximum(double& k, double& s, int chn); // Find cross-section max.
  vector<int> parseIds(string line);           // Parse the ID strings.
  vector<double> parseParms(string line);      // Parse the parameter strings.

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_DeuteronProduction_H
