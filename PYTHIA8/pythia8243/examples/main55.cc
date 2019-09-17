// main55.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how to set up the simulation of isotropic QED gamma + gamma
// production at 750 GeV by modifying gamma + gamma -> H -> gamma + gamma.
// Contains a semirealistic analysis of associated jet production.
// Exemplifies how an lhagrid1 PDF file can be used for hard processes.
// Cross section not to be taken seriously, but event properties realistic.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 10000;

  // Generator and collision energy.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");

  // Pick gamma + gamma -> H -> gamma + gamma production.
  pythia.readString("HiggsSM:gmgm2H = on");
  pythia.readString("25:onMode = off");
  pythia.readString("25:onIfMatch = 22 22");

  // Force H to have 750 GeV mass and yet be narrow.
  pythia.readString("25:m0 = 750.");
  pythia.readString("25:mWidth = 1.");
  pythia.readString("25:doForceWidth = on");

  // Use NLO CT14qed central member for hard process, default LO for rest.
  pythia.readString("PDF:useHard = on");
  pythia.readString("PDF:pHardSet = LHAGrid1:CT14qed_proton_0000.dat");

  // Do not allow gamma -> fermion pair by showering.
  // (To simplify removal of photons for jet finding in the rest of the event.)
  pythia.readString("TimeShower:QEDshowerByGamma = off");

  // No event printout.
  pythia.readString("Next:numberShowEvent = 0");

  // Initialize. Shorthand for event.
  pythia.init();
  Event& event = pythia.event;

  // Set up anti-kT jet finder with R = 0.5, pT > 30., |eta| < 5.
  SlowJet slowJet( -1, 0.5, 30., 5.);

  // Histogram.
  Hist massH( "mass of Higgs", 100, 500., 1000.);
  Hist pTH(  "pT of Higgs", 100, 0., 500.);
  Hist pTj1( "pT of jet 1", 100, 0., 500.);
  Hist pTj2( "pT of jet 2", 100, 0., 500.);
  Hist yj1(   "y of jet 1", 100, -5., 5.);
  Hist yj2(   "y of jet 2", 100, -5., 5.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find location of (final copy of) Higgs. Histogram Higgs mass and pT.
    int iH = 0;
    for (int i = 1; i < event.size(); ++i) if (event[i].id() == 25) iH = i;
    massH.fill( event[iH].m() );
    pTH.fill( event[iH].pT() );

    // Mark the two daughter photons decayed to do jet finding without them.
    event[event[iH].daughter1()].statusNeg();
    event[event[iH].daughter2()].statusNeg();
    slowJet.analyze( pythia.event );

    // Histogram pT and rapidity of two hardest jets, if they exist.
    if (slowJet.sizeJet() > 0) {
      pTj1.fill( slowJet.pT(0) );
      yj1.fill( slowJet.y(0) );
    }
    if (slowJet.sizeJet() > 1) {
      pTj2.fill( slowJet.pT(1) );
      yj2.fill( slowJet.y(1) );
    }

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  cout << massH << pTH << pTj1 << pTj2 << yj1 << yj2;

  // Done.
  return 0;
}
