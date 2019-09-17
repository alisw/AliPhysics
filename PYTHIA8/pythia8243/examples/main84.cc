// main84.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do CKKW-L merging, see the Matrix Element
// Merging page in the online manual. An example command is
//     ./main84 main84.cmnd hepmcout84.dat 2 w+_production_lhc histout84.dat
// where main84.cmnd supplies the commands, hepmcout84.dat is the
// HepMC output, 2 is the maximial number of jets, w+_production_lhc
// provides the input LHE events, and histout84.dat is the output
// histogram file. This example requires FastJet and HepMC.

#include <time.h>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;

// Functions for histogramming
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/D0RunIIConePlugin.hh"

//==========================================================================

// Find the Durham kT separation of the clustering from
// nJetMin --> nJetMin-1 jets in te input event

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

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 6) {
    cerr << " Unexpected number of command-line arguments. \n You are"
         << " expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Name of output HepMC file" << endl
         << " 3. Maximal number of additional jets"
         << " (not used internally in Pythia, only used to construct the full"
         << " name of lhe files with additional jets, and to label output"
         << " histograms)" << endl
         << " 4. Full name of the input LHE file (with path"
         << " , without any _0.lhe suffix)" << endl
         << " 5. Path for output histogram files" << endl
         << " Program stopped. " << endl;
    return 1;
  }

  Pythia pythia;

  // First argument: Get input from an input file
  pythia.readFile(argv[1]);
  int nEvent = pythia.mode("Main:numberOfEvents");

  // Interface for conversion from Pythia8::Event to HepMC event.
  // Will fill cross section and event weight directly in this program,
  // so switch it off for normal conversion routine.
  HepMC::Pythia8ToHepMC ToHepMC;
  ToHepMC.set_store_xsec(false);

  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);

  // Third argument: Maximal number of additional jets
  int njet = atoi(argv[3]);

  // Read input and output paths
  string iPath = string(argv[4]);
  string oPath = string(argv[5]);

  // To write correctly normalized events to hepmc file, first get
  // a reasonable accurate of the cross section
  int njetCounterEstimate = njet;
  vector<double> xsecEstimate;

  vector<double> nTrialEstimate;
  vector<double> nAcceptEstimate;

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 42390964");

  while(njetCounterEstimate >= 0) {

    // Number of runs
    int nRun = 1;
    double nTrial = 0.;
    double nAccept = 0.;

    int countEvents = 0;

    // Run pythia nRun times with the same lhe file to get nRun times
    // higher statistics in the histograms
    for(int n = 1; n <= nRun ; ++n ) {

      // Get process and events from LHE file, initialize only the
      // first time
      if(n > 1) pythia.readString("Main:LHEFskipInit = on");

      // From njet, choose LHE file
      stringstream in;
      in   << "_" << njetCounterEstimate << ".lhe";

      string LHEfile = iPath + in.str();

      pythia.readString("HadronLevel:all = off");

      // Read in ME configurations
      pythia.readString("Beams:frameType = 4");
      pythia.readString("Beams:LHEF = " + LHEfile);
      pythia.init();

      for( int iEvent=0; iEvent<nEvent; ++iEvent ){
        countEvents++;

        nTrial += 1.;
        if(iEvent == 0) pythia.stat();

        // Generate next event
        if(pythia.next()) nAccept += 1.;

        if(countEvents == nEvent*nRun-1){
          xsecEstimate.push_back(pythia.info.sigmaGen());
          nTrialEstimate.push_back(nTrial+1.);
          nAcceptEstimate.push_back(nAccept+1.);
        }


      } // end loop over events to generate

    } // end outer loop to rerun pythia with the same lhe file

    // Restart with ME of a reduced the number of jets
    if( njetCounterEstimate > 0 )
      njetCounterEstimate--;
    else
      break;

  } // end loop over different jet multiplicities

  cout << endl << "Finished estimating cross section"
    << endl;

  for(int i=0; i < int(xsecEstimate.size()); ++i)
    cout << "  Cross section estimate for " << njet-i << " jets :"
      << scientific << setprecision(8) << xsecEstimate[i]
      << endl;
  for(int i=0; i < int(nTrialEstimate.size()); ++i)
    cout << "  Trial events for " << njet-i << " jets :"
      << scientific << setprecision(3) << nTrialEstimate[i]
      << "  Accepted events for " << njet-i << " jets :"
      << scientific << setprecision(3) << nAcceptEstimate[i]
      << endl;
  cout << endl;

  // Now start merging procedure
  int njetCounter = njet;

  Hist histPTFirstSum("pT of first jet",100,0.,100.);
  Hist histPTSecondSum("pT of second jet",100,0.,100.);

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 42390964");

  // Sum of event weights
  double sigma = 0.0;
  double sigma2 = 0.0;

  while(njetCounter >= 0) {

    cout << "   Path to lhe files: " << iPath << "_*" << endl;
    cout << "   Output written to: " << oPath << "'name'.dat" << endl;

    // Set up histograms of pT of the first jet
    Hist histPTFirst("pT of first jet",100,0.,200.);
    Hist histPTSecond("pT of second jet",100,0.,200.);
    Hist histPTThird("pT of third jet",100,0.,200.);
    Hist histPTFourth("pT of fourth jet",50,0.,100.);
    Hist histPTFifth("pT of fifth jet",30,0.,50.);
    Hist histPTSixth("pT of sixth jet",30,0.,50.);

    // Number of runs
    int nRun = 1;
    // Number of tried events
    int nTriedEvents = 0;
    // Number of accepted events
    int nAccepEvents = 0;

    // Run pythia nRun times with the same lhe file to get nRun times
    // higher statistics in the histograms
    for(int n = 1; n <= nRun ; ++n ) {

      // Get process and events from LHE file, initialize only the
      // first time
      if(n > 1) pythia.readString("Main:LHEFskipInit = on");

      // From njet, choose LHE file
      stringstream in;
      in   << "_" << njetCounter << ".lhe";

      string LHEfile = iPath + in.str();

      cout << endl << endl
        << "\t LHE FILE FOR + " << njetCounter
        << " JET SAMPLE READ FROM " << LHEfile
        << endl << endl;

      cout << "Normalise with xsection " << xsecEstimate[njet-njetCounter]
        << endl << endl;

      pythia.readString("HadronLevel:all = on");

      // Read in ME configurations
      pythia.readString("Beams:frameType = 4");
      pythia.readString("Beams:LHEF = " + LHEfile);
      pythia.init();

      for( int iEvent=0; iEvent<nEvent; ++iEvent ){

        nTriedEvents++;
        if(iEvent == 0) pythia.stat();

        // Generate next event
        if( pythia.next()) {

          double weight = pythia.info.mergingWeight();
          nAccepEvents++;

          // Jet pT's
          double D = 0.4;
          double pTfirst = pTfirstJet(pythia.event,1, D);
          double pTsecnd = pTfirstJet(pythia.event,2, D);
          double pTthird = pTfirstJet(pythia.event,3, D);
          double pTfourt = pTfirstJet(pythia.event,4, D);
          double pTfifth = pTfirstJet(pythia.event,5, D);
          double pTsixth = pTfirstJet(pythia.event,6, D);
          histPTFirst.fill( pTfirst, weight);
          histPTSecond.fill( pTsecnd, weight);
          histPTThird.fill( pTthird, weight);
          histPTFourth.fill( pTfourt, weight);
          histPTFifth.fill( pTfifth, weight);
          histPTSixth.fill( pTsixth, weight);

          if(weight > 0.){
            // Construct new empty HepMC event and fill it.
            // Units will be as chosen for HepMC build, but can be changed
            // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
            HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();

            double normhepmc = 1.* xsecEstimate[njet-njetCounter]
                * nTrialEstimate[njet-njetCounter]
                / nAcceptEstimate[njet-njetCounter]
                * 1./ (double(nRun)*double(nEvent));

            sigma += weight*normhepmc;
            sigma2 += pow(weight*normhepmc,2);
            // Set event weight
            hepmcevt->weights().push_back(weight*normhepmc);

            // Fill summed histograms
            histPTFirstSum.fill( pTfirst, weight*normhepmc);
            histPTSecondSum.fill( pTsecnd, weight*normhepmc);

            // Fill HepMC event, with PDF info.
            ToHepMC.fill_next_event( pythia, hepmcevt );

            // Report cross section to hepmc
            HepMC::GenCrossSection xsec;
            xsec.set_cross_section( sigma*1e9,
              pythia.info.sigmaErr()*1e9 );
            hepmcevt->set_cross_section( xsec );

            // Write the HepMC event to file. Done with it.
            ascii_io << hepmcevt;
            delete hepmcevt;
          }

        } // if( pythia.next() )

        if(nTriedEvents%10000 == 0)
          cout << nTriedEvents << endl;

      } // end loop over events to generate

      // print cross section, errors
      pythia.stat();

    } // end outer loop to rerun pythia with the same lhe file

    // Normalise histograms for this particular multiplicity
    double norm = 1.
                * pythia.info.sigmaGen()
                * double(nTriedEvents)/double(nAccepEvents)
                * 1./ (double(nRun)*double(nEvent));

    histPTFirst           *= norm;
    histPTSecond          *= norm;
    histPTThird           *= norm;
    histPTFourth          *= norm;
    histPTFifth           *= norm;
    histPTSixth           *= norm;

    // Write histograms for this particular multiplicity to file
    ofstream write;
    stringstream suffix;
    suffix << njet << "_" << njetCounter;
    suffix << "_wv.dat";

    write.open( (char*)(oPath + "PTjet1_" + suffix.str()).c_str());
    histPTFirst.table(write);
    write.close();

    write.open( (char*)(oPath + "PTjet2_" + suffix.str()).c_str());
    histPTSecond.table(write);
    write.close();

    write.open( (char*)(oPath + "PTjet3_" + suffix.str()).c_str());
    histPTThird.table(write);
    write.close();

    write.open( (char*)(oPath + "PTjet4_" + suffix.str()).c_str());
    histPTFourth.table(write);
    write.close();

    write.open( (char*)(oPath + "PTjet5_" + suffix.str()).c_str());
    histPTFifth.table(write);
    write.close();

    write.open( (char*)(oPath + "PTjet6_" + suffix.str()).c_str());
    histPTSixth.table(write);
    write.close();

    histPTFirst.null();
    histPTSecond.null();
    histPTThird.null();
    histPTFourth.null();
    histPTFifth.null();
    histPTSixth.null();

    // Restart with ME of a reduced the number of jets
    if( njetCounter > 0 )
      njetCounter--;
    else
      break;

  } // end loop over different jet multiplicities

  // Since the histograms have been filled with the correct weight for
  // each jet multiplicity, no normalisation is needed.
  // Write summed histograms to file.
  ofstream writeSum;
  stringstream suffixSum;
  suffixSum << njet << "_wv.dat";

  writeSum.open( (char*)(oPath + "PTjet1Sum_" + suffixSum.str()).c_str());
  histPTFirstSum.table(writeSum);
  writeSum.close();

  writeSum.open( (char*)(oPath + "PTjet2Sum_" + suffixSum.str()).c_str());
  histPTSecondSum.table(writeSum);
  writeSum.close();

  for(int i=0; i < int(xsecEstimate.size()); ++i)
    cout << "  Cross section estimate for " << njet-i << " jets :"
      << scientific << setprecision(8) << xsecEstimate[i]
      << endl;
  for(int i=0; i < int(nTrialEstimate.size()); ++i)
    cout << "  Trial events for " << njet-i << " jets :"
      << scientific << setprecision(3) << nTrialEstimate[i]
      << "  Accepted events for " << njet-i << " jets :"
      << scientific << setprecision(3) << nAcceptEstimate[i]
      << endl;
  cout << endl;

  cout << "Histogrammed cross section for "
     << iPath << " with " << njet << " additional jets is "
     << scientific << setprecision(8) << sigma
     << " error " << sqrt(sigma2) << endl;

  return 0;
}
