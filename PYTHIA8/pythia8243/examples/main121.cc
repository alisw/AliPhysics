// main121.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Illustrate how to set up automatic uncertainty band calculations.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//--------------------------------------------------------------------------

// Small helper function to get variation names.

string weightLabel(string weightString) {
  // Strip leading whitespace
  weightString.erase(0,weightString.find_first_not_of(" \t\n\r\f\v"));
  // Find first blank and use this to isolate weight label.
  int iBlank = weightString.find(" ", 0);
  return weightString.substr(0, iBlank);
}

//--------------------------------------------------------------------------

int main() {

  // Initialize Pythia.
  Pythia pythia;
  pythia.readFile("main121.cmnd");
  pythia.init();

  // Define multiple histograms, one for each variation.
  int nWeights = pythia.info.nWeights();
  vector<double> sumOfWeights;
  vector<Hist> pTtop, nCh;
  vector<string> names;
  vector<string> weightStrings = pythia.settings.wvec("UncertaintyBands:List");

  // Loop through weights to initialize the histograms.
  for (int iWeight=0; iWeight < nWeights; ++iWeight) {
    names.push_back( (iWeight==0)
      ? "baseline" : weightLabel(weightStrings[iWeight-1]));
    pTtop.push_back ( Hist("top transverse momentum",       100,  0., 200.));
    nCh.push_back   ( Hist("charged particle multiplicity", 100, -1., 399.));
    sumOfWeights.push_back(0.);
  }

  // Event generation loop.
  int nEvent = pythia.mode("Main:numberOfEvents");
  for (int iEvent=0; iEvent < nEvent; ++iEvent) {

    // Generate next event. Break out of event loop if at end of an LHE file.
    if ( !pythia.next() ) {
      if ( pythia.info.atEndOfFile() ) break;
      else continue;
    }

    // Find last top in event record and count number of charged particles.
    int iTop = 0;
    int nChg = 0;
    for (int i=0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == 6) iTop = i;
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) ++nChg;
    }
    // Get top pT.
    double pT = pythia.event[iTop].pT();

    // Fill histograms with variation weights.
    for (int iWeight = 0; iWeight < nWeights; ++iWeight) {
      // Get weight
      double w = pythia.info.weight(iWeight);
      // Add the weight of the current event to the wsum of weights.
      sumOfWeights[iWeight]  += w;
      // Fill histograms.
      pTtop[iWeight].fill(pT, w);
      nCh[iWeight].fill(nChg, w);
    }

  }

  // Print cross section information.
  pythia.stat();

  // Normalize histograms and print data tables. Also, construct a simple
  // gnuplot command to plot the histograms.
  stringstream gnuplotCommand;
  gnuplotCommand << "gnuplot -e \"plot";
  for (int iWeight=0; iWeight < nWeights; ++iWeight) {
    cout << "Normalize histogram for weight " << names[iWeight] << " to "
         << scientific << setprecision(8) << sumOfWeights[iWeight] << endl;
    pTtop[iWeight] *= pythia.info.sigmaGen() / sumOfWeights[iWeight]
                    / 2.; // ... and divide by bin width.
    nCh[iWeight]   *= pythia.info.sigmaGen() / sumOfWeights[iWeight]
                    / 4.; // ... and divide by bin width.
    // Print data tables.
    ofstream write;
    write.open( ("pTtop-" + names[iWeight] + ".dat").c_str());
    pTtop[iWeight].table(write);
    write.close();
    write.open( ("nCh-" + names[iWeight] + ".dat").c_str());
    nCh[iWeight].table(write);
    write.close();
    // Construct gnuplot command.
    gnuplotCommand << " \'pTtop-" << names[iWeight]
    << ".dat\' using 1:2 w steps title \'" << names[iWeight] << "\',";
  }
  gnuplotCommand.seekp(-1, std::ios_base::end);
  gnuplotCommand << "; pause -1\"";

  // Print suggested gnuplot command.
  cout << "You can plot the variation by e.g. using the command:" << endl;
  cout << gnuplotCommand.str() << endl << endl;

  // Done.
  return 0;

}
