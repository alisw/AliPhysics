// main13.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how two Les Houches Event Files can be combined in PYTHIA,
// just like in main12.cc, but here with the difference that information is
// stored in main13.cmnd and read out using the subruns possibility.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // Book histogram.
  Hist nCharged("charged particle multiplicity",100,-0.5,399.5);

  // Generator.
  Pythia pythia;

  // Read in subrun-independent data from main13.cmnd.
  pythia.readFile( "main13.cmnd", 0);

  // Extract data to be used in main program. Set counters.
  int nSubrun = pythia.mode("Main:numberOfSubruns");
  int nAbort  = pythia.mode("Main:timesAllowErrors");
  int iAbort  = 0;

  // Begin loop over subruns.
  for (int iSubrun = 1; iSubrun <= nSubrun; ++iSubrun) {

    // Read in subrun-specific data from main13.cmnd.
    pythia.readFile( "main13.cmnd", iSubrun);

    // Initialize generator.
    pythia.init();

    // Begin infinite event loop - to be exited at end of file.
    for (int iEvent = 0; ; ++iEvent) {

      // Generate next event.
      if (!pythia.next()) {

        // Leave event loop if at end of file.
        if (pythia.info.atEndOfFile()) break;

        // First few failures write off as "acceptable" errors, then quit.
        if (++iAbort < nAbort) continue;
        break;
      }

      // Sum up final charged multiplicity and fill in histogram.
      int nChg = 0;
      for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) ++nChg;
      nCharged.fill(nChg);

    // End of event loop.
    }

  // End of subrun loop.
  }

  // Give statistics. Print histogram.
  pythia.stat();
  cout << nCharged;

  // Done.
  return 0;
}
