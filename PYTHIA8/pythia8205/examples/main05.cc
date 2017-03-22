// main05.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It studies jet production at the LHC, using SlowJet and CellJet.
// Note: the two finders are intended to construct approximately the same
// jet properties, but provides output in slightly different format,
// and have here not been optimized to show maximum possible agreement.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events, generated and listed ones.
  int nEvent    = 500;
  int nListJets = 5;

  // Generator. LHC process and output selection. Initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 200.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();

  // Common parameters for the two jet finders.
  double etaMax   = 4.;
  double radius   = 0.7;
  double pTjetMin = 10.;
  // Exclude neutrinos (and other invisible) from study.
  int    nSel     = 2;
  // Range and granularity of CellJet jet finder.
  int    nEta     = 80;
  int    nPhi     = 64;

  // Set up SlowJet jet finder, with anti-kT clustering
  // and pion mass assumed for non-photons..
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

  // Set up CellJet jet finder.
  CellJet cellJet( etaMax, nEta, nPhi, nSel);

  // Histograms. Note similarity in names, even when the two jet finders
  // do not calculate identically the same property (pT vs. ET, y vs. eta).
  Hist nJetsS("number of jets, SlowJet", 50, -0.5, 49.5);
  Hist nJetsC("number of jets, CellJet", 50, -0.5, 49.5);
  Hist nJetsD("number of jets, CellJet - SlowJet", 45, -22.5, 22.5);
  Hist eTjetsS("pT for jets, SlowJet", 100, 0., 500.);
  Hist eTjetsC("eT for jets, CellJet", 100, 0., 500.);
  Hist etaJetsS("y for jets, SlowJet", 100, -5., 5.);
  Hist etaJetsC("eta for jets, CellJet", 100, -5., 5.);
  Hist phiJetsS("phi for jets, SlowJwt", 100, -M_PI, M_PI);
  Hist phiJetsC("phi for jets, CellJet", 100, -M_PI, M_PI);
  Hist distJetsS("R distance between jets, SlowJet", 100, 0., 10.);
  Hist distJetsC("R distance between jets, CellJet", 100, 0., 10.);
  Hist eTdiffS("pT difference, SlowJet", 100, -100., 400.);
  Hist eTdiffC("eT difference, CellJet", 100, -100., 400.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Analyze Slowet jet properties. List first few.
    slowJet. analyze( pythia.event );
    if (iEvent < nListJets) slowJet.list();

    // Fill SlowJet inclusive jet distributions.
    nJetsS.fill( slowJet.sizeJet() );
    for (int i = 0; i < slowJet.sizeJet(); ++i) {
      eTjetsS.fill( slowJet.pT(i) );
      etaJetsS.fill( slowJet.y(i) );
      phiJetsS.fill( slowJet.phi(i) );
    }

    // Fill SlowJet distance between jets.
    for (int i = 0; i < slowJet.sizeJet() - 1; ++i)
    for (int j = i +1; j < slowJet.sizeJet(); ++j) {
      double dEta = slowJet.y(i) - slowJet.y(j);
      double dPhi = abs( slowJet.phi(i) - slowJet.phi(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dEta) + pow2(dPhi) );
      distJetsS.fill( dR );
    }

    // Fill SlowJet pT-difference between jets (to check ordering of list).
    for (int i = 1; i < slowJet.sizeJet(); ++i)
      eTdiffS.fill( slowJet.pT(i-1) - slowJet.pT(i) );

    // Analyze CellJet jet properties. List first few.
    cellJet. analyze( pythia.event, pTjetMin, radius );
    if (iEvent < nListJets) cellJet.list();

    // Fill CellJet inclusive jet distributions.
    nJetsC.fill( cellJet.size() );
    for (int i = 0; i < cellJet.size(); ++i) {
      eTjetsC.fill( cellJet.eT(i) );
      etaJetsC.fill( cellJet.etaWeighted(i) );
      phiJetsC.fill( cellJet.phiWeighted(i) );
    }

    // Fill CellJet distance between jets.
    for (int i = 0; i < cellJet.size() - 1; ++i)
    for (int j = i +1; j < cellJet.size(); ++j) {
      double dEta = cellJet.etaWeighted(i)
        - cellJet.etaWeighted(j);
      double dPhi = abs( cellJet.phiWeighted(i)
        - cellJet.phiWeighted(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dEta) + pow2(dPhi) );
      distJetsC.fill( dR );
    }

    // Fill CellJet ET-difference between jets (to check ordering of list).
    for (int i = 1; i < cellJet.size(); ++i)
      eTdiffC.fill( cellJet.eT(i-1) - cellJet.eT(i) );

    // Compare number of jets for the two finders.
    nJetsD.fill( cellJet.size() - slowJet.sizeJet() );

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  cout << nJetsS << nJetsC << nJetsD << eTjetsS << eTjetsC
       << etaJetsS << etaJetsC << phiJetsS << phiJetsC
       << distJetsS << distJetsC << eTdiffS << eTdiffC;

  // Done.
  return 0;
}
