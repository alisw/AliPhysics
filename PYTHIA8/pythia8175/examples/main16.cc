// main16.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. 
// It illustrates (a) how to collect the analysis code in a separate class
// and (b) how to provide the .cmnd filename on the command line

// Once you have linked the main program you can run it with a command line
// ./main16.exe main16.cmnd > out16

#include "Pythia.h"

using namespace Pythia8; 

//==========================================================================

// Put all your own analysis code in the myAnalysis class. 

class MyAnalysis {

public:

  // Constructor can be empty.
  MyAnalysis() {}

  // Initialization actions.
  void init();
 
  // Analysis of each new event.
  void analyze(Event& event);

  // Show final results.
  void finish();

private:

  // Declare variables and objects that span init - analyze - finish.
  int  nEvt;
  Hist yH, etaChg, mult; 

};

//--------------------------------------------------------------------------

// The initialization code. 

void MyAnalysis::init() {

  // Initialize counter for number of events.
  nEvt = 0;

  // Book histograms.
  yH.book("Higgs rapidity", 100, -10., 10.);
  etaChg.book("charged pseudorapidity", 100, -10., 10.);
  mult.book( "charged multiplicity", 100, -0.5, 799.5);

} 

//--------------------------------------------------------------------------

// The event analysis code. 

void MyAnalysis::analyze(Event& event) {

  // Increase counter.
  ++nEvt;

  // Find latest copy of Higgs and plot its rapidity.
  int iH = 0;
  for (int i = 0; i < event.size(); ++i) 
    if (event[i].id() == 25) iH = i;
  yH.fill( event[iH].y() );

  // Plot pseudorapidity distribution. Sum up charged multiplicity.
  int nChg = 0;
  for (int i = 0; i < event.size(); ++i) 
  if (event[i].isFinal() && event[i].isCharged()) {
    etaChg.fill( event[i].eta() );
    ++nChg;
  }
  mult.fill( nChg );

} 

//--------------------------------------------------------------------------

// The finishing code. 

void MyAnalysis::finish() {

  // Normalize histograms.
  double binFactor = 5. / nEvt;
  yH     *= binFactor;
  etaChg *= binFactor;

  // Print histograms.
  cout << yH << etaChg << mult;

} 

//==========================================================================

// You should not need to touch the main program: its actions are 
// determined by the .cmnd file and the rest belongs in MyAnalysis.

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 2) {
    cerr << " Unexpected number of command-line arguments. \n"
         << " You are expected to provide a file name and nothing else. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Check that the provided file name corresponds to an existing file.
  ifstream is(argv[1]);  
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Confirm that external file will be used for settings..
  cout << " PYTHIA settings will be read from file " << argv[1] << endl;

  // Declare generator. Read in commands from external file.
  Pythia pythia;
  pythia.readFile(argv[1]);
 
  // Initialization.
  pythia.init();

  // Declare user analysis class. Do initialization part of it.
  MyAnalysis myAnalysis;
  myAnalysis.init(); 

  // Read in number of event and maximal number of aborts.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Begin event loop.
  int iAbort = 0; 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    // User Analysis of current event.
    myAnalysis.analyze( pythia.event);

  // End of event loop.
  }

  // Final statistics.
  pythia.stat();

  // User finishing.
  myAnalysis.finish();

  // Done.
  return 0;
}
