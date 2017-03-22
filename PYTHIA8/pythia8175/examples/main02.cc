// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.

#include "Pythia.h"
using namespace Pythia8; 
int main() {
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:idB = -2212");    
  pythia.readString("Beams:eCM = 1960.");    
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");    
  pythia.readString("PhaseSpace:mHatMin = 80.");    
  pythia.readString("PhaseSpace:mHatMax = 120.");    
  pythia.init();
  Hist pTZ("dN/dpTZ", 100, 0., 100.);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;
    // Loop over particles in event. Find last Z0 copy. Fill its pT. 
    int iZ = 0;
    for (int i = 0; i < pythia.event.size(); ++i) 
      if (pythia.event[i].id() == 23) iZ = i;
    pTZ.fill( pythia.event[iZ].pT() );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << pTZ; 
  return 0;
}
