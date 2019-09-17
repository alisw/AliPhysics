// main63.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how you can use UserHooks to enhance rare emission rates,
// in this case q -> q gamma.
// To concentrate on the photons from the showers, MPI and hadronization
// are switched off by default.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Write own derived UserHooks class.

class EnhanceHooks : public UserHooks {

public:

  // Constructor and destructor do nothing.
  EnhanceHooks() {}
  ~EnhanceHooks() {}

  // Enhance real-emission rate. Thus no trial-emission enhancement.
  bool canEnhanceEmission() { return true;}
  bool canEnhanceTrial()    { return false;}

  // Function to return the weight enhance factor.
  double enhanceFactor(string name) {
    if (name == "isr:Q2QA") return 50.;
    return 1.0;
  }

  // Function to return the vetoing probability.
  double vetoProbability(string name) {
    if (name == "isr:Q2QA") return 0.5;
    return 0.0;
  }

};

//==========================================================================

int main() {

  // Histogram pT spectrum of photons and event weights.
  Hist gamNoEnh(   "gamma pT spectrum, no enhancement",   100, 0., 100.);
  Hist gamWithEnh( "gamma pT spectrum, with enhancement", 100, 0., 100.);
  Hist gamRatio("gamma pT spectrum, with/no enhancement", 100, 0., 100.);
  Hist gamBefWt(   "gamma pT spectrum, without weight",   100, 0., 100.);
  Hist eventWt(    "log10(event weight)",                 100, -7., 3.);

  // Compare generation without and with enhanced q -> q gamma emission.
  for (int iCase = 0; iCase < 2; ++iCase) {

    // Generator.
    Pythia pythia;
    pythia.readFile("main63.cmnd");
    int nEvent = pythia.mode("Main:numberOfEvents");

    // Set up a user hook and send it in.
    UserHooks* enhanceHooks = 0;
    if (iCase == 1) {
      enhanceHooks = new EnhanceHooks();
      pythia.setUserHooksPtr( enhanceHooks);
    }

    // LHC initialization.
    pythia.init();

    // Begin event loop.
    double sumWt = 0.;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events. Find and histogram event weight.
      pythia.next();
      double weight = (iCase == 1)
                    ? enhanceHooks->getEnhancedEventWeight() : 1.;
      if (iCase == 1) eventWt.fill( log10(weight) );
      sumWt += weight;

      // Find all final-state photons and histogram them.
      for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].id() == 22) {
        double pT = pythia.event[i].pT();
        if (iCase == 0) gamNoEnh.fill(   pT, 1.);
        if (iCase == 1) gamBefWt.fill(   pT, 1.);
        if (iCase == 1) gamWithEnh.fill( pT, weight);
      }

    // End of event loop.
    }

    // Statistics.
    pythia.stat();
    cout << "\n Average event weight = " << scientific
         << sumWt / nEvent << endl;

    // End of case loop.
    if (iCase == 1) delete enhanceHooks;
  }

  // Write histograms to output stream.
  gamRatio = gamWithEnh / gamNoEnh;
  cout << gamNoEnh << gamWithEnh << gamRatio << gamBefWt << eventWt;

  // Write histogram data to files.
  ofstream write;
  write.open("PTA_0");
  gamNoEnh.table(write);
  write.close();
  write.open("PTA_1");
  gamWithEnh.table(write);
  write.close();

  // Done.
  return 0;
}
