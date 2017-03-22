// main80.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do CKKW-L merging,
// see the Matrix Element Merging page in the online manual.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator. Input parameters.
  Pythia pythia;
  pythia.readFile("main80.cmnd");

  // Extract number of events and max number of jets in merging.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nMerge = pythia.mode("Merging:nJetMax");

  // Histograms combined over all jet multiplicities.
  Hist pTWsum("pT of W, summed over all subruns", 100, 0., 200.);

  // Merged total cross section, summed over subruns.
  double sigmaTotal = 0.;

  // Loop over subruns with varying number of jets.
  for (int iMerge = 0; iMerge <= nMerge; ++iMerge) {
    double sigmaSample = 0.;

    // Read in name of LHE file for current subrun and initialize.
    pythia.readFile("main80.cmnd", iMerge);
    pythia.init();

    // Histograms for current jet multiplicity.
    Hist weightNow("event weights, current subrun", 100, 0., 2.5);
    Hist pTWnow("pT of W, current subrun", 100, 0., 200.);

    // Start event generation loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate next event. Break out of event loop if at end of LHE file.
      if ( !pythia.next() ) {
        if ( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get CKKWL weight of current event. Histogram and accumulate it.
      double weight = pythia.info.mergingWeight();
      weightNow.fill( weight);
      sigmaSample += weight;

      // Find the final copy of the W+, which is after the full shower.
      int iW = 0;
      for (int i = 1; i < pythia.event.size(); ++i)
        if (pythia.event[i].id() == 24) iW = i;

      // Fill the pT of the W histogram, with CKKWL weight.
      double pTW = pythia.event[iW].pT();
      pTWnow.fill( pTW, weight);

    // End of event loop.
    }

    // Normalize pTW histogram, convert mb -> pb, and correct for bin width.
    pTWnow *= 1e9 * pythia.info.sigmaGen() / (2. * pythia.info.nAccepted());

    // Print cross section and histograms for current subrun.
    pythia.stat();
    cout << weightNow << pTWnow;

    // Sum up merged cross section of current run.
    sigmaSample *= pythia.info.sigmaGen() / double(pythia.info.nAccepted());
    sigmaTotal  += sigmaSample;

    // Add current histogram to the combined one. End of subrun loop.
    pTWsum += pTWnow;
  }

  // Print final histograms and info on merged cross section..
  cout << pTWsum;
  cout << "\n\n The inclusive cross section after merging is: "
       << scientific << setprecision(4) << sigmaTotal << " mb " << endl;

  // Done
  return 0;

}
