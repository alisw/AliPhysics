// main37.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how Les Houches Event File version 3.0 input can be used
// to mix events according to several different event weights.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator
  Pythia pythia;

  // Initialize Les Houches Event File run.
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = wbj_lhef3.lhe");
  pythia.init();

  // Get number of event weights.
  int ninitrwgt = pythia.info.getInitrwgtSize();

  // Initialise as many histograms as there are event weights.
  vector<Hist> pTw;
  for (int iHist = 0; iHist < ninitrwgt; ++iHist)
    pTw.push_back( Hist("pT W",50,0.,200.) );

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent < 10; ++iEvent) {

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      else continue;
    }

    // Find the final copy of the W and its pT.
    int iW = 0;
    for (int i = pythia.event.size() - 1; i > 0; --i)
      if (pythia.event[i].idAbs() == 24) { iW = i; break;}
    double pT = pythia.event[iW].pT();

    // Loop over the event weights in the detailed format and histogram.
    unsigned int nwgt = pythia.info.getWeightsDetailedSize();
    for (unsigned int iwgt = 0; iwgt < nwgt; ++iwgt) {
      // The weights happen to have the identifiers 1001, 1002...
      string key;
      ostringstream convert;
      convert << iwgt + 1001;
      key = convert.str();
      double w = pythia.info.getWeightsDetailedValue(key);
      pTw[iwgt].fill( max(pT,0.5), w );
    }

  // End of event loop.
  }

  // Give statistics.
  pythia.stat();

  // Print histograms.
  ofstream write;
  stringstream suffix;
  for (int iHist = 0; iHist < ninitrwgt; ++iHist) {
    suffix << iHist << ".dat";
    // Write histograms to file.
    write.open( (char*)("PTW_" + suffix.str()).c_str());
    pTw[iHist].table(write);
    suffix.str("");
    write.close();
  }

  // Done.
  return 0;
}
