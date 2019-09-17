// main03.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how different processes can be selected and studied.
// All input is specified in the main03.cmnd file.
// Also illustrated output to be plotted by Python/Matplotlib/pyplot.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

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
  double sigma = pythia.info.sigmaGen();
  pThard   *= sigma * 1e6 * 0.2 / nEvent;
  nFinal   *= 1. / (16. * nEvent);
  nCharged *= 1. / (8. * nEvent);
  dndy     *=  5. / nEvent;
  dndeta   *=  5. / nEvent;
  dndpT    *= 10. / nEvent;
  cout << pThard << nFinal << nCharged << dndy << dndeta << dndpT << epCons;

  // Write Python code that can generate a PDF file with the spectra.
  // Assuming you have Python installed on your platform, do as follows.
  // After the program has run, type "python main03plot.py" (without the " ")
  // in a terminal window, and open "out03plot.pdf" in a PDF viewer.
  // Colours and other choices can be omitted, but are shown as illustration.
  HistPlot hpl("main03plot");
  hpl.plotFrame("out03plot", pThard, "$p_{\\perp}$ scale of hard interaction",
    "$p_{\\perp}$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}p_{\\perp}$ (nb/GeV)",
    "h", "$p_{\\perp}$ of $2 \\to 2$ process", true);
  hpl.frame("", "Total and charged particle multiplicities",
    "$n$", "$\\mathrm{d}P/\\mathrm{d}n$");
  hpl.add( nFinal, "h,royalblue", "total");
  hpl.add( nCharged, "h,orange", "charged (even only!)");
  hpl.plot();
  hpl.frame( "", "Charged (pseudo)rapidity distribution", "$y$ or $\\eta$",
    "$\\mathrm{d}n_{\\mathrm{charged}}/\\mathrm{d}(y/\\eta)$");
  hpl.add( dndy, "-", "$\\mathrm{d}n_{\\mathrm{charged}}/\\mathrm{d}y$");
  hpl.add( dndeta, "--,magenta",
    "$\\mathrm{d}n_{\\mathrm{charged}}/\\mathrm{d}\\eta$");
  hpl.plot();
  hpl.plotFrame("", dndpT, "Charged $p_{\\perp}$ spectrum",
    "$p_{\\perp}$ (GeV)", "$\\mathrm{d}n_{\\mathrm{charged}}/\\mathrm{d}"
    "p_{\\perp}$ (GeV$^{-1}$)", "", "charged", true);

  // Done.
  return 0;
}
