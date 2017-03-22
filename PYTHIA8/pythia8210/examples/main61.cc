// main61.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
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
  Event& event = pythia.event;
  Info&  info  = pythia.info;

  // Set it up to generate Z's at 8 TeV.
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:mMin = 70.");
  pythia.readString("23:mMax = 110.");

  // Setup of diffractive framework.
  pythia.readString("Diffraction:doHard = on");
  pythia.readString("Diffraction:sampleType = 1");
  pythia.readString("Diffraction:PomFlux = 5");
  pythia.readString("PDF:PomSet = 6");

  // Simplify printout.
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Init:showMultipartonInteractions = off");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:showScaleAndVertex = off");

  // Switch off hadronization, since not used here.
  pythia.readString("HadronLevel:all = off");

  // Initialize.
  pythia.init();

  // Collect information on the number of diffractive events
  int maxEvent      = 10000;
  int nDiffA        = 0;
  int nDiffB        = 0;
  int nReducedDiffA = 0;
  int nReducedDiffB = 0;

  // Histograms.
  Hist y0("dN/dy inclusive",               100, -5.,   5.);
  Hist y1("dN/dy after PDF selection",     100, -5.,   5.);
  Hist y2("dN/dy after MPI selection",     100, -5.,   5.);
  Hist pT0("dN/dpTZ inclusive",            100,  0., 100.);
  Hist pT1("dN/dpTZ after PDF selection",  100,  0., 100.);
  Hist pT2("dN/dpTZ after MPI selection",  100,  0., 100.);
  Hist xP1("dN/dxPom after PDF selection", 100,  0.,   1.);
  Hist xP2("dN/dxPom after MPI selection", 100,  0.,   1.);
  Hist tP1("dN/dt after PDF selection",    100, -2.,   0.);
  Hist tP2("dN/dt after MPI selection",    100, -2.,   0.);

  // Begin event loop. Generate event; skip if generation failed.
  for (int iEvent = 0; iEvent < maxEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Locate the Z0 and find its y and pT.
    int iZ = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].id() == 23) iZ = i;
    double yZ  = event[iZ].y();
    double pTZ = event[iZ].pT();
    y0.fill( yZ );
    pT0.fill( pTZ );

    // Find diffractive events. Histogram y and pT.
    if ( info.isHardDiffractiveA() == 1 || info.isHardDiffractiveB() == 1) {
      y1.fill( yZ );
      pT1.fill( pTZ );
      if (info.nMPI() == 1) {
        y2.fill( yZ );
        pT2.fill( pTZ );
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
  cout << y0 << y1 << y2 << pT0 << pT1 << pT2 << xP1 << xP2 << tP1 << tP2;

  return 0;
}
