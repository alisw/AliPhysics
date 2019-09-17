// main19.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program runs four instances of Pythia simultaneously,
// one for signal events, one for pileup background ones, and two
// For beam-gas background ones. Note that Pythia does not do nuclear
// effects, so beam-gas is represented by "fixed-target" pp collisions.
// The = and += overloaded operators are used to join several
// event records into one, but should be used with caution.

// The possibility to instantiate Pythia with Settings and ParticleData
// databases is illustrated, but not essential here. It means that the
// share/Pythia8/xmldoc/*.xml files are only read once, saving some time.

// Note that each instance of Pythia is running independently of any other,
// but with two important points to remember.
// 1) By default all generate the same random number sequence,
//    which has to be corrected if they are to generate the same
//    physics, like the two beam-gas ones below.
// 2) Interfaces to external Fortran programs are "by definition" static.
//    Thus it is not a good idea to use LHAPDF5 to set different PDF's
//    in different instances.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Method to pick a number according to a Poissonian distribution.

int poisson(double nAvg, Rndm& rndm) {

  // Set maximum to avoid overflow.
  const int NMAX = 100;

  // Random number.
  double rPoisson = rndm.flat() * exp(nAvg);

  // Initialize.
  double rSum  = 0.;
  double rTerm = 1.;

  // Add to sum and check whether done.
  for (int i = 0; i < NMAX; ) {
    rSum += rTerm;
    if (rSum > rPoisson) return i;

    // Evaluate next term.
    ++i;
    rTerm *= nAvg / i;
  }

  // Emergency return.
  return NMAX;
}

//==========================================================================

int main() {

  // Number of signal events to generate.
  int nEvent = 100;

  // Beam Energy.
  double eBeam = 7000.;

  // Average number of pileup events per signal event.
  double nPileupAvg = 2.5;

  // Average number of beam-gas events per signal ones, on two sides.
  double nBeamAGasAvg = 0.5;
  double nBeamBGasAvg = 0.5;

  // Signal generator instance.
  Pythia pythiaSignal;

  // Switch off automatic event listing (illustrates settings inheritance).
  pythiaSignal.readString("Next:numberShowInfo = 0");
  pythiaSignal.readString("Next:numberShowProcess = 0");
  pythiaSignal.readString("Next:numberShowEvent = 0");

  // Switch off K0S decay (illustrates particle data inheritance).
  pythiaSignal.readString("130:mayDecay = off");

  // Background generator instances copies settings and particle data.
  Pythia pythiaPileup(   pythiaSignal.settings, pythiaSignal.particleData);
  Pythia pythiaBeamAGas( pythiaSignal.settings, pythiaSignal.particleData);
  Pythia pythiaBeamBGas( pythiaSignal.settings, pythiaSignal.particleData);

  // Switch off Lambda decay (illustrates particle data non-inheritance).
  pythiaSignal.readString("3122:mayDecay = off");

  // One object where all individual events are to be collected.
  Event sumEvent;

  // Initialize generator for signal processes.
  pythiaSignal.readString("HardQCD:all = on");
  pythiaSignal.readString("PhaseSpace:pTHatMin = 50.");
  pythiaSignal.settings.parm("Beams:eCM", 2. * eBeam);
  pythiaSignal.init();

  // Initialize generator for pileup (background) processes.
  pythiaPileup.readString("Random:setSeed = on");
  pythiaPileup.readString("Random:seed = 10000002");
  pythiaPileup.readString("SoftQCD:all = on");
  pythiaPileup.settings.parm("Beams:eCM", 2. * eBeam);
  pythiaPileup.init();

  // Initialize generators for beam A - gas (background) processes.
  pythiaBeamAGas.readString("Random:setSeed = on");
  pythiaBeamAGas.readString("Random:seed = 10000003");
  pythiaBeamAGas.readString("SoftQCD:all = on");
  pythiaBeamAGas.readString("Beams:frameType = 2");
  pythiaBeamAGas.settings.parm("Beams:eA", eBeam);
  pythiaBeamAGas.settings.parm("Beams:eB", 0.);
  pythiaBeamAGas.init();

  // Initialize generators for beam B - gas (background) processes.
  pythiaBeamBGas.readString("Random:setSeed = on");
  pythiaBeamBGas.readString("Random:seed = 10000004");
  pythiaBeamBGas.readString("SoftQCD:all = on");
  pythiaBeamBGas.readString("Beams:frameType = 2");
  pythiaBeamBGas.settings.parm("Beams:eA", 0.);
  pythiaBeamBGas.settings.parm("Beams:eB", eBeam);
  pythiaBeamBGas.init();

  // Histograms: number of pileups, total charged multiplicity.
  Hist nPileH("number of pileup events per signal event", 100, -0.5, 99.5);
  Hist nAGH("number of beam A + gas events per signal event", 100, -0.5, 99.5);
  Hist nBGH("number of beam B + gas events per signal event", 100, -0.5, 99.5);
  Hist nChgH("number of charged multiplicity",100, -0.5, 1999.5);
  Hist sumPZH("total pZ of system",100, -100000., 100000.);

  // Loop over events.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate a signal event. Copy this event into sumEvent.
    if (!pythiaSignal.next()) continue;
    sumEvent = pythiaSignal.event;

    // Select the number of pileup events to generate.
    int nPileup = poisson(nPileupAvg, pythiaPileup.rndm);
    nPileH.fill( nPileup );

    // Generate a number of pileup events. Add them to sumEvent.
    for (int iPileup = 0; iPileup < nPileup; ++iPileup) {
      pythiaPileup.next();
      sumEvent += pythiaPileup.event;
    }

    // Select the number of beam A + gas events to generate.
    int nBeamAGas = poisson(nBeamAGasAvg, pythiaBeamAGas.rndm);
    nAGH.fill( nBeamAGas );

    // Generate a number of beam A + gas events. Add them to sumEvent.
    for (int iAG = 0; iAG < nBeamAGas; ++iAG) {
      pythiaBeamAGas.next();
      sumEvent += pythiaBeamAGas.event;
    }

    // Select the number of beam B + gas events to generate.
    int nBeamBGas = poisson(nBeamBGasAvg, pythiaBeamBGas.rndm);
    nBGH.fill( nBeamBGas );

    // Generate a number of beam B + gas events. Add them to sumEvent.
    for (int iBG = 0; iBG < nBeamBGas; ++iBG) {
      pythiaBeamBGas.next();
      sumEvent += pythiaBeamBGas.event;
    }

    // List first few events.
    if (iEvent < 1) {
      pythiaSignal.info.list();
      pythiaSignal.process.list();
      sumEvent.list();
    }

    // Find charged multiplicity.
    int nChg = 0;
    for (int i = 0; i < sumEvent.size(); ++i)
      if (sumEvent[i].isFinal() && sumEvent[i].isCharged()) ++nChg;
    nChgH.fill( nChg );

    // Fill net pZ - nonvanishing owing to beam + gas.
    sumPZH.fill( sumEvent[0].pz() );

  // End of event loop
  }

  // Statistics. Histograms.
  pythiaSignal.stat();
  pythiaPileup.stat();
  pythiaBeamAGas.stat();
  pythiaBeamBGas.stat();
  cout << nPileH << nAGH << nBGH << nChgH << sumPZH;

  return 0;
}
