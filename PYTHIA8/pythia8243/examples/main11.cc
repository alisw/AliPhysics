// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how Les Houches Event File input can be used in Pythia8.
// It uses the ttsample.lhe(.gz) input file, the latter only with 100 events.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

  // You can always read an plain LHE file,
  // but if you ran "./configure --with-gzip" before "make"
  // then you can also read a gzipped LHE file.
#ifdef GZIPSUPPORT
  bool useGzip = true;
#else
  bool useGzip = false;
#endif
  cout << " useGzip = " << useGzip << endl;

  // Generator. We here stick with default values, but changes
  // could be inserted with readString or readFile.
  Pythia pythia;

  // Initialize Les Houches Event File run. List initialization information.
  pythia.readString("Beams:frameType = 4");
  if (useGzip) pythia.readString("Beams:LHEF = ttbar.lhe.gz");
  else         pythia.readString("Beams:LHEF = ttbar.lhe");
  pythia.init();

  // Book histogram.
  Hist nCharged("charged particle multiplicity",100,-0.5,399.5);

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; ; ++iEvent) {

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
