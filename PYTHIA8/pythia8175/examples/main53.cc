// main53.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. 
// It illustrates how to interface an external process with an incoming photon
// in a hadron beam, using the MRST2004QED PDF set.
// All input apart from the name of the external LHEF file is specified in the 
// main53.cmnd file.

#include "Pythia.h"

using namespace Pythia8; 

int main() {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("main53.cmnd");    

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");

  // Initialize. Either of two opions, to be picked in main53.cmnd.
  // 1) Read in external event with incoming photon in the ME,
  // from pre-generated .lhe file (thanks to SANC and R. Sadykov).
  // 2) Use internal fermion gamma -> W+- fermion' process.
  pythia.init();

  // Histograms for pT distribution in gluon production vertex.
  Hist pTprim( "pT of photon production, no ISR", 100, 0., 100.);
  Hist pTwith( "pT of photon production, with ISR", 100, 0., 100.);

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      break;
    }

    // Analyze event to find branching where photon is produced.
    int iGam = (event[3].id() == 22) ? 3 : 4;
    int iGamMother = iGam;
    for ( ; ; ) {
      iGamMother = event[iGam].mother1();
      if (iGamMother < iGam || event[iGamMother].id() != 22) break;
      iGam = iGamMother;
    }

    // Find and histogram pT in this branching.
    if (iGamMother < iGam) pTprim.fill( event[iGam].pT() );
    else {
      int iQ = iGamMother;
      int size = event.size();
      do ++iQ; 
      while (event[iQ].status() != -43 && iQ < size - 1);
      if (event[iQ].status() == -43) pTwith.fill( event[iQ].pT() );
    }

  // End of event loop.
  }

  // Final statistics and histogram output.
  pythia.stat();
  cout << pTprim << pTwith;

  return 0;
}

