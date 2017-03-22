// main92.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// Modified by Rene Brun and Axel Naumann to put the Pythia::event
// into a TTree.

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"

using namespace Pythia8;

int main() {

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Beams:eCM = 14000.");
  pythia.init();

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open("pytree.root","recreate");
  Event *event = &pythia.event;
  TTree *T = new TTree("T","ev1 Tree");
  T->Branch("event",&event);

 // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;

    // Fill the pythia event into the TTree.
    // Warning: the files will rapidly become large if all events
    // are saved. In some cases it may be convenient to do some
    // processing of events and only save those that appear
    // interesting for future analyses.
    T->Fill();

  // End event loop.
  }

  // Statistics on event generation.
  pythia.stat();

  //  Write tree.
  T->Print();
  T->Write();
  delete file;

  // Done.
  return 0;
}
