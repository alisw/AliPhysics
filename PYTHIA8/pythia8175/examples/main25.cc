// main25.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. 
// It illustrates how Les Houches Event File input can be used in Pythia8.
// Here the very few events are generated with MadGraph, and illustrate
// more complicated colour topologies.

#include "Pythia.h"
using namespace Pythia8; 

int main() {

  // Generator           
  Pythia pythia;                            

  // Stick with default values, so do not bother with a separate file
  // for changes. However, do one change, to show readString in action.
  pythia.readString("PartonLevel:ISR = off"); 
  pythia.readString("PartonLevel:FSR = off"); 
  pythia.readString("PartonLevel:MPI = off"); 
  pythia.readString("HadronLevel:Hadronize = on"); 

  // Initialize Les Houches Event File run.
  pythia.readString("Beams:frameType = 4"); 
  pythia.readString("Beams:LHEF = main25.lhe"); 
  pythia.init();

  // Book histogram.
  Hist nCharged("charged particle multiplicity",100,-0.5,399.5); 

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.     
  for (int iEvent = 0; ; ++iEvent) {
    cout << endl << "Begin event # " << iEvent << endl;

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break; 

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    // Sum up final charged multiplicity and fill in histogram.
    int nChg = 0;                 
    for (int i = 0; i < pythia.event.size(); ++i) 
    if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) 
      ++nChg;
    nCharged.fill(nChg);               

  // End of event loop.        
  }                                           

  // Give statistics. Print histogram.
  pythia.stat();
  cout << nCharged;  

  // Done.                           
  return 0;
}
