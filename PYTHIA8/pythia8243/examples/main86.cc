// main86.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do UMEPS merging, see the Matrix Element
// Merging page in the online manual. An example command is
//     ./main86 main86.cmnd w_production hepmcout86.dat
// where main86.cmnd supplies the commands, w_production provides the
// input LHE events, and hepmcout86.dat is the output file. This
// example requires HepMC.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include <unistd.h>

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate merging

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Name of the input LHE file (with path), up to the '_tree'"
         << " identifier" << endl
         << " 3. Path for output HepMC" << endl
         << " Program stopped. " << endl;
    return 1;
  }

  Pythia pythia;

  // Input parameters:
  pythia.readFile(argv[1]);
  // Interface for conversion from Pythia8::Event to HepMC one.
  HepMC::Pythia8ToHepMC ToHepMC;
  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(argv[3], std::ios::out);
  // Switch off warnings for parton-level events.
  ToHepMC.set_print_inconsistency(false);
  ToHepMC.set_free_parton_exception(false);
  // Do not store cross section information, as this will be done manually.
  ToHepMC.set_store_pdf(false);
  ToHepMC.set_store_proc(false);
  ToHepMC.set_store_xsec(false);

  // Path to input events, with name up to the "_tree" identifier included.
  string iPath = string(argv[2]);

  // Number of events
  int nEvent = pythia.mode("Main:numberOfEvents");
  // Maximal number of additional LO jets.
  int nMaxLO =  pythia.mode("Merging:nJetMax");

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  // Switch off all showering and MPI when extimating the cross section after
  // the merging scale cut.
  bool fsr = pythia.flag("PartonLevel:FSR");
  bool isr = pythia.flag("PartonLevel:ISR");
  bool mpi = pythia.flag("PartonLevel:MPI");
  bool had = pythia.flag("HadronLevel:all");
  pythia.settings.flag("PartonLevel:FSR",false);
  pythia.settings.flag("PartonLevel:ISR",false);
  pythia.settings.flag("HadronLevel:all",false);
  pythia.settings.flag("PartonLevel:MPI",false);

  // Switch on cross section estimation procedure.
  pythia.settings.flag("Merging:doXSectionEstimate", true);
  pythia.settings.flag("Merging:doUMEPSTree",true);

  int njetcounterLO = nMaxLO;
  string iPathTree  = iPath + "_tree";

  // Save estimates in vectors.
  vector<double> xsecLO;
  vector<double> nAcceptLO;

  cout << endl << endl << endl;
  cout << "Start estimating umeps tree level cross section" << endl;

  while(njetcounterLO >= 0) {

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterLO << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathTree+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathTree + in.str();

    pythia.settings.mode("Merging:nRequested", njetcounterLO);
    pythia.settings.mode("Beams:frameType", 4);
    pythia.settings.word("Beams:LHEF", LHEfile);
    pythia.init();

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){
      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ){
          break;
        }
        else continue;
      }
    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();

    xsecLO.push_back(pythia.info.sigmaGen());
    nAcceptLO.push_back(pythia.info.nAccepted());

    // Restart with ME of a reduced the number of jets
    if( njetcounterLO > 0 )
      njetcounterLO--;
    else
      break;

  } // end loop over different jet multiplicities


  // Switch off cross section estimation.
  pythia.settings.flag("Merging:doXSectionEstimate", false);

  // Switch showering and multiple interaction back on.
  pythia.settings.flag("PartonLevel:FSR",fsr);
  pythia.settings.flag("PartonLevel:ISR",isr);
  pythia.settings.flag("HadronLevel:all",had);
  pythia.settings.flag("PartonLevel:MPI",mpi);

  // Declare sample cross section for output.
  double sigmaTemp  = 0.;
  vector<double> sampleXStree;
  vector<double> sampleXSsubtTree;
  // Cross section an error.
  double sigmaTotal  = 0.;
  double errorTotal  = 0.;

  int sizeLO    = int(xsecLO.size());
  njetcounterLO = nMaxLO;
  iPathTree     = iPath + "_tree";

  while(njetcounterLO >= 0){

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterLO << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathTree+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathTree + in.str();

    pythia.settings.flag("Merging:doUMEPSTree",true);
    pythia.settings.flag("Merging:doUMEPSSubt",false);
    pythia.settings.mode("Merging:nRecluster",0);

    cout << endl << endl << endl
         << "Start tree level treatment for " << njetcounterLO << " jets"
         << endl;

    pythia.settings.mode("Merging:nRequested", njetcounterLO);
    pythia.settings.mode("Beams:frameType", 4);
    pythia.settings.word("Beams:LHEF", LHEfile);
    pythia.init();

    // Remember position in vector of cross section estimates.
    int iNow = sizeLO-1-njetcounterLO;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weight = pythia.info.mergingWeight();
      double evtweight = pythia.info.weight();
      weight *= evtweight;
      // Do not print zero-weight events.
      if ( weight == 0. ) continue;

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = xsecLO[iNow] / nAcceptLO[iNow];
      // Set event weight
      hepmcevt->weights().push_back(weight*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weight*normhepmc;
      sigmaTemp  += weight*normhepmc;
      errorTotal += pow2(weight*normhepmc);
      // Report cross section to hepmc
      HepMC::GenCrossSection xsec;
      xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt->set_cross_section( xsec );
      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;

    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();
    // Save sample cross section for output.
    sampleXStree.push_back(sigmaTemp);
    sigmaTemp = 0.;

    // Restart with ME of a reduced the number of jets
    if( njetcounterLO > 0 )
      njetcounterLO--;
    else
      break;

  }

  cout << endl << endl << endl;
  cout << "Do UMEPS subtraction" << endl;

  int njetcounterLS   = nMaxLO;
  string iPathSubt    = iPath + "_tree";

  while(njetcounterLS >= 1){

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterLS << ".lhe";
#ifdef GZIPSUPPORT
    if(access( (iPathSubt+in.str()+".gz").c_str(), F_OK) != -1) in << ".gz";
#endif
    string LHEfile = iPathSubt + in.str();

    pythia.settings.flag("Merging:doUMEPSTree",false);
    pythia.settings.flag("Merging:doUMEPSSubt",true);
    pythia.settings.mode("Merging:nRecluster",1);

    cout << endl << endl << endl
         << "Start subtractive treatment for " << njetcounterLS << " jets"
         << endl;

    pythia.settings.mode("Merging:nRequested", njetcounterLS);
    pythia.settings.mode("Beams:frameType", 4);
    pythia.settings.word("Beams:LHEF", LHEfile);
    pythia.init();

    // Remember position in vector of cross section estimates.
    int iNow = sizeLO-1-njetcounterLS;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weight    = pythia.info.mergingWeight();
      double evtweight = pythia.info.weight();
      weight *= evtweight;
      // Do not print zero-weight events.
      if ( weight == 0. ) continue;

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = -1*xsecLO[iNow] / nAcceptLO[iNow];
      // Set event weight
      hepmcevt->weights().push_back(weight*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weight*normhepmc;
      sigmaTemp  += weight*normhepmc;
      errorTotal += pow2(weight*normhepmc);
      // Report cross section to hepmc.
      HepMC::GenCrossSection xsec;
      xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt->set_cross_section( xsec );
      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;

    } // end loop over events to generate

    // print cross section, errors
    pythia.stat();
    // Save sample cross section for output.
    sampleXSsubtTree.push_back(sigmaTemp);
    sigmaTemp = 0.;

    // Restart with ME of a reduced the number of jets
    if( njetcounterLS > 1 )
      njetcounterLS--;
    else
      break;
  }

  // Print cross section information.
  cout << endl << endl;
  cout << " *---------------------------------------------------*" << endl;
  cout << " |                                                   |" << endl;
  cout << " | Sample cross sections after UMEPS merging         |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | Leading order cross sections (mb):                |" << endl;
  for (int i = 0; i < int(sampleXStree.size()); ++i)
    cout << " |     " << sampleXStree.size()-1-i << "-jet:  "
         << setw(17) << scientific << setprecision(6)
         << sampleXStree[i] << "                     |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | Leading-order subtractive cross sections (mb):    |" << endl;
  for (int i = 0; i < int(sampleXSsubtTree.size()); ++i)
    cout << " |     " << sampleXSsubtTree.size()-1-i+1 << "-jet:  "
         << setw(17) << scientific << setprecision(6)
         << sampleXSsubtTree[i] << "                     |" << endl;
  cout << " |                                                   |" << endl;
  cout << " |---------------------------------------------------|" << endl;
  cout << " |---------------------------------------------------|" << endl;
  cout << " | Inclusive cross sections:                         |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | UMEPS merged inclusive cross section:             |" << endl;
  cout << " |    " << setw(17) << scientific << setprecision(6)
       << sigmaTotal << "  +-  " << setw(17) << sqrt(errorTotal) << " mb "
       << "   |" << endl;
  cout << " |                                                   |" << endl;
  cout << " | LO inclusive cross section:                       |" << endl;
  cout << " |    " << setw(17) << scientific << setprecision(6)
       << xsecLO.back() << " mb                           |" << endl;
  cout << " |                                                   |" << endl;
  cout << " *---------------------------------------------------*" << endl;
  cout << endl << endl;

  // Done
  return 0;

}
