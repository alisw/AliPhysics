// main09.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Generate a predetermined second hard interaction.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;
  Event& event   = pythia.event;

  // Select first hard process (just a small sample of possibilities).
  //pythia.readString("HardQCD:all = on");
  pythia.readString("Top:all = on");
  //pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  //pythia.readString("WeakSingleBoson:ffbar2W = on");

  // Select second hard process (complete list of options).
  pythia.readString("SecondHard:generate = on");
  //pythia.readString("SecondHard:TwoJets = on");
  pythia.readString("SecondHard:PhotonAndJet = on");
  //pythia.readString("SecondHard:TwoPhotons = on");
  //pythia.readString("SecondHard:SingleGmZ = on");
  //pythia.readString("SecondHard:SingleW = on");
  //pythia.readString("SecondHard:TwoBJets = on");

  // Kinematics cuts, common for the two.
  pythia.readString("PhaseSpace:mHatMin = 40.");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  // Initialize for LHC at 8 TeV.
  pythia.readString("Beams:eCM = 8000.");
  pythia.init();

  // Histogram.
  Hist pTfirst("pT first collision",    100, 0., 400.);
  Hist pTsecond("pT second collision",  100, 0., 200.);
  Hist pTdiff("pT first-second collision", 100, -100., 300.);
  Hist nMult("number of multiparton interactions", 100, -0.5, 99.5);
  Hist bMore("b enhancement factor",    100, 0., 10.);
  Hist nChg("charged multiplicity", 100, -0.5, 999.5);

  // Generate events.
  for (int iev = 0; iev < 1000; ++iev) {
    pythia.next();

    // Histogram pT.
    double pT1 = pythia.info.pTMPI(0);
    double pT2 = pythia.info.pTMPI(1);
    pTfirst.fill( pT1 );
    pTsecond.fill( pT2 );
    pTdiff.fill( pT1 - pT2 );

    // Histogram multiparton interactions
    double nMPI = pythia.info.nMPI();
    nMult.fill( nMPI );
    bMore.fill( pythia.info.enhanceMPI() );

    // Histogram charged multiplicity.
    int nCharged = 0;
    for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isCharged()) ++nCharged;
    nChg.fill( nCharged );

  }

  // Compare full statistics listing with what is set in info.
  pythia.stat();
  cout << scientific << setprecision(3) << "\n From pythia.info: sigma = "
       << pythia.info.sigmaGen() << " +- " << pythia.info.sigmaErr()
       << endl;

  // Print histograms.
  cout << pTfirst << pTsecond << pTdiff << nMult << bMore << nChg;

  // Done.
  return 0;
}
