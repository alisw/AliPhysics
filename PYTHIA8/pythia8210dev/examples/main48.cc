// main48.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Philip Ilten.

// An example where decays are performed with EvtGen. See the
// documentation in Pythia8Plugins/EvtGen.h for further details on the
// EvtGenDecays class. In this example the decay B0 -> nu_e e+ D*-[->
// gamma D-[-> nu_ebar e- pi0]] is forced. The invariant mass spectrum
// of the final state electron and pion is then plotted.

// The syntax to run this example is:
//     ./main48 <EvtGen decay file> <EvtGen particle data> <PYTHIA8DATA>
//              <flag to use EvtGen>

// This example has only been tested with EvtGen version 1.3.0. The
// EvtGen package is designed to use Pythia 8 for any decays that it
// cannot perform. For this to be possible, EvtGen must be linked
// against the Pythia 8 shared library. To build EvtGen 1.3.0 with
// Pythia 8.2 the "-llhapdfdummy" library requirement must be
// removed. Prior to running "./configure" for EvtGen this can be
// accomplished via:
//     sed -i "s/-llhapdfdummy//g" configure
// To modify how this example program is compiled (i.e to remove
// linking against the EvtGenExternal library) change the main48 rule
// in the Makefile of this directory.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/EvtGen.h"
using namespace Pythia8;
int main(int argc, char* argv[]) {

  // Check arguments.
  if (argc != 5) {
    cerr << " Unexpected number of command-line arguments. \n You are"
         << " expected to provide the arguments \n"
         << " 1. EvtGen decay file (e.g. DECAY_2010.DEC) \n"
         << " 2. EvtGen particle data (e.g. evt.pdl) \n"
         << " 3. PYTHIA8DATA path \n"
         << " 4. Flag to use EvtGen (true or false) \n"
         << " Program stopped. " << endl;
    return 1;
  }
  bool use(string(argv[4]) == "true");

  // Intialize Pythia.
  Pythia pythia;
  pythia.readString("Print:quiet = on");
  pythia.readString("HardQCD:hardbbbar = on");
  if (!use) {
    cout << "Not using EvtGen." << endl;
    pythia.readString("511:onMode = off");
    pythia.readString("511:onIfMatch = 12 -11 -413");
    pythia.readString("413:onMode = off");
    pythia.readString("413:onIfMatch = 411 22");
    pythia.readString("411:onMode = off");
    pythia.readString("411:onIfMatch = -11 12 111");
  } else cout << "Using EvtGen." << endl;
  pythia.init();

  // Initialize EvtGen.
  EvtGenDecays *evtgen = 0;
  if (use) {
    setenv("PYTHIA8DATA", argv[3], 1);
    evtgen = new EvtGenDecays(&pythia, argv[1], argv[2]);
    evtgen->readDecayFile("main48.dec");
  }

  // The event loop.
  Hist mass("m(e, pi0) [GeV]", 100, 0., 2.);
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    // Generate the event.
    if (!pythia.next()) continue;
    // Perform the decays with EvtGen.
    if (evtgen) evtgen->decay();
    // Analyze the event.
    Event &event = pythia.event;
    for (int iPrt = 0; iPrt < (int)event.size(); ++iPrt) {
      if (event[iPrt].idAbs() != 511) continue;
      int iB0(event[iPrt].iBotCopyId()), iDsm(-1), iDm(-1), iE(-1), iPi0(-1);
      for (int iDtr = event[iB0].daughter1(); iDtr <= event[iB0].daughter2();
           ++ iDtr) {
        if (event[iDtr].idAbs() == 413) {
          iDsm = event[iDtr].iBotCopyId();
          continue;
        }
      }
      if (iDsm == -1) continue;
      for (int iDtr = event[iDsm].daughter1(); iDtr <= event[iDsm].daughter2();
           ++ iDtr) {
        if (event[iDtr].idAbs() == 411) {
          iDm = event[iDtr].iBotCopyId();
          continue;
        }
      }
      if (iDm == -1) continue;
      for (int iDtr = event[iDm].daughter1(); iDtr <= event[iDm].daughter2();
           ++ iDtr) {
        if (event[iDtr].idAbs() == 11)  iE   = event[iDtr].iBotCopyId();
        if (event[iDtr].idAbs() == 111) iPi0 = event[iDtr].iBotCopyId();
      }
      if (iE == -1 || iPi0 == -1) continue;
      mass.fill((event[iE].p() + event[iPi0].p()).mCalc());
    }
  }

  // Print the statistics and histogram.
  pythia.stat();
  mass /= mass.getEntries();
  cout << mass;
  if (evtgen) delete evtgen;
  return 0;
}
