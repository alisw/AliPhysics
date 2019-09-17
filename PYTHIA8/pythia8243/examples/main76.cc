// main76.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to run SUSY processes in Pythia8.
// All input is specified in the main76.cmnd file.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("main76.cmnd");

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors");

  // Initialize. Set lifetime of long-lived particle.
  pythia.init();
  int iLLP = 56; // Change to 57 for X+ and 59 for X++
  cout << "Lifetime [mm] = " << scientific << pythia.particleData.tau0(iLLP)
       << endl;

  // Histograms.
  Hist life("Decay lifetime [mm]",100,0.,100.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      event.list();
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    life.fill(event[5].tau());

  // End of event loop.
  }

  // Final statistics and histogram output.
  pythia.stat();
  cout << life << endl;

  return 0;
}
