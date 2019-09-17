// main41.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch, based on main01.cc.
// This program illustrates how HepMC can be interfaced to Pythia8.
// It studies the charged multiplicity distribution at the LHC.
// HepMC events are output to the hepmcout41.dat file.

// WARNING: typically one needs 25 MB/100 events at the LHC.
// Therefore large event samples may be impractical.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;

int main() {

  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io("hepmcout41.dat", std::ios::out);

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );

    // Construct new empty HepMC event and fill it.
    // Units will be as chosen for HepMC build; but can be changed
    // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( pythia, hepmcevt );

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << mult;

  // Done.
  return 0;
}
