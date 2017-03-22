// main81.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do CKKW-L merging, see the Matrix Element
// Merging page in the online manual. An example command is
//     ./main81 main81.cmnd w+_production_lhc_0.lhe histout81.dat
// where main81.cmnd supplies the commands, w+_production_lhc_0.lhe
// provides the input LHE events, and histout81.dat is the output
// file. This example requires FastJet.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

// Functions for histogramming
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/D0RunIIConePlugin.hh"

//==========================================================================

// Find the Durham kT separation of the clustering from
// nJetMin --> nJetMin-1 jets in the input event

double pTfirstJet( const Event& event, int nJetMin, double Rparam) {

  double yPartonMax = 4.;

  // Fastjet analysis - select algorithm and parameters
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  // For hadronic collision, use hadronic Durham kT measure
  if(event[3].colType() != 0 || event[4].colType() != 0)
    jetDef = new fastjet::JetDefinition(fastjet::kt_algorithm, Rparam,
                                      recombScheme, strategy);
  // For e+e- collision, use e+e- Durham kT measure
  else
    jetDef = new fastjet::JetDefinition(fastjet::ee_kt_algorithm,
                                      recombScheme, strategy);
  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;
  // Reset Fastjet input
  fjInputs.resize(0);

  // Loop over event record to decide what to pass to FastJet
  for (int i = 0; i < event.size(); ++i) {
    // (Final state && coloured+photons) only!
    if ( !event[i].isFinal()
      || event[i].isLepton()
      || event[i].id() == 23
      || abs(event[i].id()) == 24
      || abs(event[i].y()) > yPartonMax)
      continue;

    // Store as input to Fastjet
    fjInputs.push_back( fastjet::PseudoJet (event[i].px(),
            event[i].py(), event[i].pz(),event[i].e() ) );
  }

  // Do nothing for empty input
  if (int(fjInputs.size()) == 0) {
    delete jetDef;
    return 0.0;
  }

  // Run Fastjet algorithm
  fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
  // Extract kT of first clustering
  double pTFirst = sqrt(clustSeq.exclusive_dmerge_max(nJetMin-1));

  delete jetDef;
  // Return kT
  return pTFirst;

}

//==========================================================================

// Example main programm to illustrate merging

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Full name of the input LHE file (with path)" << endl
         << " 3. Path for output histogram files" << endl
         << " Program stopped. " << endl;
    return 1;
  }

  Pythia pythia;

  // Input parameters:
  //  1. Input file for settings
  //  2. Path to input LHE file
  //  3. Output histogram path
  pythia.readFile(argv[1]);
  string iPath = string(argv[2]);
  string oPath = string(argv[3]);
  // Number of events
  int nEvent = pythia.mode("Main:numberOfEvents");

  // For ISR regularisation off
  pythia.settings.forceParm("SpaceShower:pT0Ref",0.);

  // Declare histograms
  Hist histPTFirst("pT of first jet",100,0.,100.);
  Hist histPTSecond("pT of second jet",100,0.,100.);
  Hist histPTThird("pT of third jet",100,0.,100.);

  // Read in ME configurations
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + iPath);
  pythia.init();

  // Start generation loop
  for( int iEvent=0; iEvent<nEvent; ++iEvent ){

    // Generate next event
    if( ! pythia.next()) continue;

    // Get CKKWL weight of current event
    double weight = pythia.info.mergingWeight();

    // Fill bins with CKKWL weight
    // Functions use fastjet to get first / second jet
    double pTfirst = pTfirstJet(pythia.event,1, 0.4);
    double pTsecnd = pTfirstJet(pythia.event,2, 0.4);
    double pTthird = pTfirstJet(pythia.event,3, 0.4);
    histPTFirst.fill( pTfirst, weight);
    histPTSecond.fill( pTsecnd, weight);
    histPTThird.fill( pTthird, weight);

    if(iEvent%1000 == 0) cout << iEvent << endl;

  } // end loop over events to generate

  // print cross section, errors
  pythia.stat();

  // Normalise histograms
  double norm = 1.
              * pythia.info.sigmaGen()
              * 1./ double(nEvent);

  histPTFirst           *= norm;
  histPTSecond          *= norm;
  histPTThird           *= norm;

  // Get the number of jets in the LHE file from the file name
  string jetsInLHEF = iPath.substr(iPath.size()-5, iPath.size());
  jetsInLHEF = jetsInLHEF.substr(0, jetsInLHEF.size()-4);

  // Write histograms to dat file. Use "jetsInLHEF" to label the files
  // Once all the samples up to the maximal desired jet multiplicity from the
  // matrix element are run, add all histograms to produce a
  // matrix-element-merged prediction

  ofstream write;
  stringstream suffix;
  suffix << jetsInLHEF << "_wv.dat";

  // Write histograms to file
  write.open( (char*)(oPath + "PTjet1_" + suffix.str()).c_str());
  histPTFirst.table(write);
  write.close();

  write.open( (char*)(oPath + "PTjet2_" + suffix.str()).c_str());
  histPTSecond.table(write);
  write.close();

  write.open( (char*)(oPath + "PTjet3_" + suffix.str()).c_str());
  histPTThird.table(write);
  write.close();

  // Done
  return 0;

}
