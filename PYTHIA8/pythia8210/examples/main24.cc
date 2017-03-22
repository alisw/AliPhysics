// main24.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to run SUSY processes in Pythia8.
// All input is specified in the main22.cmnd file.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("main24.cmnd");

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors");
  double eCM   = pythia.parm("Beams:eCM");

  // Initialize.
  pythia.init();

  // Histograms.
  double epTol = 1e-6 * eCM;
  Hist epCons("deviation from energy-momentum conservation",100,0.,epTol);
  Hist nFinal("final particle multiplicity",100,-0.5,799.5);
  Hist dnparticledy("dn/dy for particles",100,-10.,10.);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Loop over final particles in the event.
    int nFin = 0;
    Vec4 pSum;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
      nFin++;
      pSum += event[i].p();
      dnparticledy.fill(event[i].y());
    }

    // Check and print event with too big energy-momentum deviation.
    nFinal.fill(nFin);
    double epDev = abs(pSum.e() - eCM) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);
    if (epDev > epTol) {
      cout << " Warning! Event with epDev = " << scientific
           << setprecision(4) << epDev << " now listed:";
      event.list();
    }

  // End of event loop.
  }

  // Final statistics and histogram output.
  pythia.stat();
  cout << epCons << nFinal << dnparticledy;

  return 0;
}
