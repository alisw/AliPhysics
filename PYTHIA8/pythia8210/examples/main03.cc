// main03.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how different processes can be selected and studied.
// All input is specified in the main03.cmnd file.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("main03.cmnd");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

  // Book histograms.
  Hist pThard("process pT scale", 100, 0., 500.);
  Hist nFinal("final particle multiplicity", 100, -0.5, 1599.5);
  Hist nCharged("charged particle multiplicity", 100, -0.5, 799.5);
  Hist dndy("dn/dy for charged particles", 100, -10., 10.);
  Hist dndeta("dn/d(eta) for charged particles", 100, -10., 10.);
  Hist dndpT("dn/dpT for charged particles", 100, 0., 10.);
  Hist epCons("deviation from energy-momentum conservation", 100, 0., 1e-6);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Fill hard scale of event.
    pThard.fill( pythia.info.pTHat() );

    // Loop over final particles in the event.
    int  nFin = 0;
    int  nChg = 0;
    Vec4 pSum;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

       // Analyze all particles.
      nFin++;
      pSum += event[i].p();

      // Analyze charged particles and fill histograms.
      if (event[i].isCharged()) {
        ++nChg;
        dndy.fill( event[i].y() );
        dndeta.fill( event[i].eta() );
        dndpT.fill( event[i].pT() );
      }

    // End of particle loop. Fill global properties.
    }
    nFinal.fill( nFin );
    nCharged.fill( nChg );
    pSum /= event[0].e();
    double epDev = abs(pSum.e() - 1.) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);

  // End of event loop.
  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();
  dndy   *=  5. / nEvent;
  dndeta *=  5. / nEvent;
  dndpT  *= 10. / nEvent;
  cout << pThard << nFinal << nCharged << dndy << dndeta << dndpT << epCons;

  // Done.
  return 0;
}
