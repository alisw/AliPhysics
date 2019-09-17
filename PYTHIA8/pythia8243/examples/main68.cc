// main68.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Christine O. Rasmussen.

// The y, pT, x_Pomeron and t distributions for forward Z bosons at the LHC,
// within the hard diffraction framework for an inclusive event sample.
// Tests the impact of successive requirements.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Create Pythia instance. Shorthand for event and info.
  Pythia pythia;
  Info&  info  = pythia.info;

  // Set it up to generate dijets at HERA.
  pythia.readString("Beams:frametype = 2");    // Beams with different energies
  pythia.readString("Beams:idA = 2212");       // p+ at 920 GeV
  pythia.readString("Beams:eA = 920");
  pythia.readString("Beams:idB = -11");        // e+ at 27.5 GeV
  pythia.readString("Beams:eB = 27.5");
  pythia.readString("PDF:lepton2gamma = on");  // Allow for photon-from lepton
  pythia.readString("Photon:ProcessType = 0"); // Allow all photon processes
  pythia.readString("Photon:Q2max = 1.");      // Maximal Q2
  pythia.readString("HardQCD:all = on");       // All dijet MEs
  pythia.readString("PhotonParton:all = on");  // All dijet MEs with photons
  pythia.readString("PhaseSpace:pThatMin = 4.");            // Minimal pT cut
  pythia.readString("MultipartonInteractions:pT0Ref = 3."); // Tuned ep value

  // Setup of diffractive framework.
  pythia.readString("Diffraction:doHard = on");
  pythia.readString("Diffraction:sampleType = 1");   // 'PDF' sample
  pythia.readString("Diffraction:hardDiffSide = 2"); // Diff. on photon side
  pythia.readString("SigmaDiffractive:PomFlux = 7"); // H1 Fit B LO
  pythia.readString("PDF:PomSet = 6");               // H1 Fit B LO

  // Simplify printout.
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Init:showMultipartonInteractions = off");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:showScaleAndVertex = off");

  // Initialize.
  pythia.init();

  // Collect information on the number of diffractive events
  int maxEvent      = 10000;
  int nListJets     = 5;
  int nDiffA        = 0;
  int nDiffB        = 0;
  int nReducedDiffA = 0;
  int nReducedDiffB = 0;

  // Histograms.
  Hist nJets0("number of jets, inclusive",            50, -0.5,   9.5);
  Hist nJets1("number of jets, after PDF selection",  50, -0.5,   9.5);
  Hist nJets2("number of jets, after MPI selection",  50, -0.5,   9.5);
  Hist eTjets0("pT for jets, inclusive",             100,  0.0,  50.0);
  Hist eTjets1("pT for jets, after PDF selection",   100,  0.0,  50.0);
  Hist eTjets2("pT for jets, after MPI selection",   100,  0.0,  50.0);
  Hist etaJets0("y for jets, inclusive",             100, -5.0,   5.0);
  Hist etaJets1("y for jets, after PDF selection",   100, -5.0,   5.0);
  Hist etaJets2("y for jets, after MPI selection",   100, -5.0,   5.0);
  Hist phiJets0("phi for jets, inclusive",           100, -M_PI, M_PI);
  Hist phiJets1("phi for jets, after PDF selection", 100, -M_PI, M_PI);
  Hist phiJets2("phi for jets, after MPI selection", 100, -M_PI, M_PI);
  Hist xP1("dN/dxPom after PDF selection",           100,  0.0,   1.0);
  Hist xP2("dN/dxPom after MPI selection",           100,  0.0,   1.0);
  Hist tP1("dN/dt after PDF selection",              100, -2.0,   0.0);
  Hist tP2("dN/dt after MPI selection",              100, -2.0,   0.0);

  // Parameters for the jet finder.
  double etaMax   = 4.;
  double radius   = 1.0;
  double pTjetMin = 3.;
  int    nSel     = 2; // Exclude neutrinos from study.

  // Set up SlowJet jet finder, with C/A clustering
  // and pion mass assumed for non-photons..
  SlowJet slowJet( 0, radius, pTjetMin, etaMax, nSel, 1);

  // Begin event loop. Generate event; skip if generation failed.
  for (int iEvent = 0; iEvent < maxEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Analyze Slowet jet properties. List first few.
    slowJet. analyze( pythia.event );
    if (iEvent < nListJets) slowJet.list();

    // Fill SlowJet inclusive jet distributions.
    nJets0.fill( slowJet.sizeJet() );
    for (int i = 0; i < slowJet.sizeJet(); ++i) {
      eTjets0.fill( slowJet.pT(i) );
      etaJets0.fill( slowJet.y(i) );
      phiJets0.fill( slowJet.phi(i) );
    }

    // Find diffractive events. Histogram y and pT.
    if ( info.isHardDiffractiveA() == 1 || info.isHardDiffractiveB() == 1) {
      nJets1.fill( slowJet.sizeJet() );
      for (int i = 1; i < slowJet.sizeJet(); ++i) {
        eTjets1.fill( slowJet.pT(i) );
        etaJets1.fill( slowJet.y(i) );
        phiJets1.fill( slowJet.phi(i) );
      }
      if (info.nMPI() == 1) {
        nJets2.fill( slowJet.sizeJet() );
        for (int i = 2; i < slowJet.sizeJet(); ++i) {
          eTjets2.fill( slowJet.pT(i) );
          etaJets2.fill( slowJet.y(i) );
          phiJets2.fill( slowJet.phi(i) );
        }
      }

      // Statistics and histogram on x_Pomeron and t.
      if ( info.isHardDiffractiveA() == 1) {
        ++nDiffA;
        xP1.fill( info.xPomeronB() );
        tP1.fill( info.tPomeronB() );
        if (info.nMPI() == 1) {
          ++nReducedDiffA;
          xP2.fill( info.xPomeronB() );
          tP2.fill( info.tPomeronB() );
        }
      }
      else if ( info.isHardDiffractiveB() == 1) {
        ++nDiffB;
        xP1.fill( info.xPomeronA() );
        tP1.fill( info.tPomeronA() );
        if (info.nMPI() == 1) {
          ++nReducedDiffB;
          xP2.fill( info.xPomeronA() );
          tP2.fill( info.tPomeronA() );
        }
      }
    }

  // End of event loop. Statistics on event generation.
  }
  pythia.stat();

  // Statistics on diffraction.
  cout << "Side A is MPI-unchecked diffractive : " << nDiffA << endl;
  cout << "Side A is MPI-checked diffractive   : " << nReducedDiffA << endl;
  cout << "Side B is MPI-unchecked diffractive : " << nDiffB << endl;
  cout << "Side B is MPI-checked diffractive   : " << nReducedDiffB << endl;
  cout << "Total MPI-unchecked diffractive events : " << fixed
       << setprecision(2) << (nDiffA + nDiffB) / double(maxEvent) * 100.
       << "%" << endl;
  cout << "Total MPI-checked diffractive events : "
       << (nReducedDiffA + nReducedDiffB) / double(maxEvent) * 100.
       << "%" << endl;

  // Histograms.
  cout << xP1 << xP2 << tP1 << tP2
       << nJets0 << nJets1 << nJets2 << eTjets0 << eTjets1 << eTjets2
       << etaJets0 << etaJets1 << etaJets2 << phiJets0 << phiJets1 << phiJets2;

  return 0;
}
