// main88.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do UNLOPS merging, 
// see the NLO Merging page in the online manual. 

#include "Pythia.h"

#include "HepMCInterface.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
// Following line to be used with HepMC 2.04 onwards.
#include "HepMC/Units.h"

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate UNLOPS merging

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Name of the input LHE file (with path), up to the '_tree'"
         << " or '_powheg' identifiers" << endl
         << " 3. Path for output histogram files" << endl
         << " Program stopped. " << endl;
    return 1;
  }

  Pythia pythia;

  // Input parameters:
  //  1. Input file for settings
  //  2. Path to input LHE file
  //  3. OUtput histogram path
  pythia.readFile(argv[1]);

  // Interface for conversion from Pythia8::Event to HepMC one. 
  HepMC::I_Pythia8 ToHepMC;
  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(argv[3], std::ios::out);
  // Switch off warnings for parton-level events.
  ToHepMC.set_print_inconsistency(false);
  ToHepMC.set_free_parton_warnings(false);
  // Do not store cross section information, as this will be done manually.
  ToHepMC.set_store_pdf(false);
  ToHepMC.set_store_proc(false);
  ToHepMC.set_store_xsec(false);

  // Path to input events, with name up to the "_tree", "_powheg" identifier
  // included.
  string iPath = string(argv[2]);

  // Number of events
  int nEvent   = pythia.mode("Main:numberOfEvents");
  // Maximal number of additional LO jets.
  int nMaxLO   = pythia.mode("Merging:nJetMax");
  // maximal number of additional NLO jets.
  int nMaxNLO  = pythia.mode("Merging:nJetMaxNLO");

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

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
  pythia.settings.flag("Merging:doUNLOPSTree", true);

  int njetcounterLO  = nMaxLO;
  string iPathTree   = iPath + "_tree";

  // Save estimates in vectors.
  vector<double> xsecLO;
  vector<double> nSelectedLO;
  vector<double> nAcceptLO;
  vector<int> strategyLO;

  cout << endl << endl << endl;
  cout << "Start estimating unlops tree level cross section" << endl;

  while(njetcounterLO >= 0){

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterLO << ".lhe";
    string LHEfile = iPathTree + in.str();

    pythia.readString("Beams:frameType = 4"); 
    pythia.settings.word("Beams:LHEF", LHEfile);  
    pythia.settings.mode("Merging:nRequested", njetcounterLO);
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
    nSelectedLO.push_back(pythia.info.nSelected());
    nAcceptLO.push_back(pythia.info.nAccepted());
    strategyLO.push_back(pythia.info.lhaStrategy());

    // Restart with ME of a reduced the number of jets
    if( njetcounterLO > 0 )
      njetcounterLO--;
    else
      break;

  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  cout << endl << endl << endl;
  cout << "Start estimating unlops virtual corrections cross section" << endl;

  pythia.settings.flag("Merging:doUNLOPSTree",false);
  pythia.settings.flag("Merging:doUNLOPSLoop", true);

  int njetcounterNLO  = nMaxNLO;
  string iPathLoop    = iPath + "_powheg";

  // Save estimates in vectors.
  vector<double> xsecNLO;
  vector<double> nSelectedNLO;
  vector<double> nAcceptNLO;
  vector<int> strategyNLO;

  while(njetcounterNLO >= 0){

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterNLO << ".lhe";
    string LHEfile = iPathLoop + in.str();

    pythia.readString("Beams:frameType = 4"); 
    pythia.settings.word("Beams:LHEF", LHEfile);  
    pythia.settings.mode("Merging:nRequested", njetcounterNLO);
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

    xsecNLO.push_back(pythia.info.sigmaGen());
    nSelectedNLO.push_back(pythia.info.nSelected());
    nAcceptNLO.push_back(pythia.info.nAccepted());
    strategyNLO.push_back(pythia.info.lhaStrategy());

    // Restart with ME of a reduced the number of jets
    if( njetcounterNLO > 0 )
      njetcounterNLO--;
    else
      break;

  }

  // Set k-factors
  int sizeLO   = int(xsecLO.size());
  int sizeNLO  = int(xsecNLO.size());
  double k0    = 1.;
  double k1    = 1.;
  double k2    = 1.;
  // Lowest order k-factor only
  if ( false ) k1 = k2 = k0 = xsecNLO.back() / xsecLO.back();
  // No k-factors
  if ( true ) k0 = k1 = k2 = 1.;

  cout << " K-Factors :" << endl;
  cout << "k0 = " << k0 << endl;
  cout << "k1 = " << k1 << endl;
  cout << "k2 = " << k2 << endl;

  // Switch off cross section estimation.
  pythia.settings.flag("Merging:doXSectionEstimate", false);

  // Switch showering and multiple interaction back on.
  pythia.settings.flag("PartonLevel:FSR",fsr);
  pythia.settings.flag("PartonLevel:ISR",isr);
  pythia.settings.flag("HadronLevel:all",had);
  pythia.settings.flag("PartonLevel:MPI",mpi);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  // Cross section an error.
  double sigmaTotal  = 0.;
  double errorTotal  = 0.;

  // Switch on tree-level processing.
  pythia.settings.flag("Merging:doUNLOPSTree",true);
  pythia.settings.flag("Merging:doUNLOPSLoop",false);
  pythia.settings.flag("Merging:doUNLOPSSubt",false);
  pythia.settings.flag("Merging:doUNLOPSSubtNLO",false);
  pythia.settings.mode("Merging:nRecluster",0);

  // Start looping through input event files.
  njetcounterLO = nMaxLO;
  iPathTree     = iPath + "_tree";

  while(njetcounterLO >= 0){

    // Set k factors
    pythia.settings.parm("Merging:kFactor0j", k0);
    pythia.settings.parm("Merging:kFactor1j", k1);
    pythia.settings.parm("Merging:kFactor2j", k2);

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterLO << ".lhe";
    string LHEfile = iPathTree + in.str();

    cout << endl << endl << endl
         << "Start tree level treatment for " << njetcounterLO << " jets"
         << endl;

    // UNLOPS does not contain a zero-jet tree-level sample.
    if ( njetcounterLO == 0 ) break;

    pythia.readString("Beams:frameType = 4"); 
    pythia.settings.word("Beams:LHEF", LHEfile);  
    pythia.settings.mode("Merging:nRequested", njetcounterLO);
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
      double weightNLO  = pythia.info.mergingWeightNLO();
      double evtweight  = pythia.info.weight();
      weightNLO        *= evtweight;
      // Do not print zero-weight events.
      if ( weightNLO == 0. ) continue; 

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = xsecLO[iNow] / nAcceptLO[iNow];
      // Set hepmc event weight.
      hepmcevt->weights().push_back(weightNLO*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weightNLO*normhepmc;
      errorTotal += pow2(weightNLO*normhepmc);
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

    // Restart with ME of a reduced the number of jets
    if( njetcounterLO > 0 )
      njetcounterLO--;
    else
      break;

  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  cout << endl << endl << endl;
  cout << "Start unlops virtual corrections part" << endl;

  // Switch on loop-level processing.
  pythia.settings.flag("Merging:doUNLOPSTree",false);
  pythia.settings.flag("Merging:doUNLOPSLoop",true);
  pythia.settings.flag("Merging:doUNLOPSSubt",false);
  pythia.settings.flag("Merging:doUNLOPSSubtNLO",false);
  pythia.settings.mode("Merging:nRecluster",0);

  njetcounterNLO = nMaxNLO;
  iPathLoop= iPath + "_powheg";

  while(njetcounterNLO >= 0){

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterNLO << ".lhe";
    string LHEfile = iPathLoop + in.str();

    cout << endl << endl << endl
         << "Start loop level treatment for " << njetcounterNLO << " jets"
         << endl;

    pythia.readString("Beams:frameType = 4"); 
    pythia.settings.word("Beams:LHEF", LHEfile);  
    pythia.settings.mode("Merging:nRequested", njetcounterNLO);
    pythia.init();
    // Remember position in vector of cross section estimates.
    int iNow = sizeNLO-1-njetcounterNLO;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weightNLO  = pythia.info.mergingWeightNLO();
      double evtweight  = pythia.info.weight();
      weightNLO        *= evtweight;
      // Do not print zero-weight events.
      if ( weightNLO == 0. ) continue;

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = xsecNLO[iNow] / nAcceptNLO[iNow];
      // powheg weighted events
      if( abs(strategyNLO[iNow]) == 4) 
        normhepmc = 1. / (1e9*nSelectedNLO[iNow]);
      // Set hepmc event weight.
      hepmcevt->weights().push_back(weightNLO*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weightNLO*normhepmc;
      errorTotal += pow2(weightNLO*normhepmc);
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

    // Restart with ME of a reduced the number of jets
    if( njetcounterNLO > 0)
      njetcounterNLO--;
    else
      break;

  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  cout << endl << endl << endl;
  cout << "Shower subtractive events" << endl;

  // Switch on processing of counter-events.
  pythia.settings.flag("Merging:doUNLOPSTree",false);
  pythia.settings.flag("Merging:doUNLOPSLoop",false);
  pythia.settings.flag("Merging:doUNLOPSSubt",true);
  pythia.settings.flag("Merging:doUNLOPSSubtNLO",false);
  pythia.settings.mode("Merging:nRecluster",1);

  int nMaxCT = nMaxLO;
  int njetcounterCT = nMaxCT;
  string iPathSubt= iPath + "_tree";

  while(njetcounterCT >= 1){

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterCT << ".lhe";
    string LHEfile = iPathSubt + in.str();

    cout << endl << endl << endl
         << "Start subtractive treatment for " << njetcounterCT << " jets"
         << endl;

    pythia.readString("Beams:frameType = 4"); 
    pythia.settings.word("Beams:LHEF", LHEfile);  
    pythia.settings.mode("Merging:nRequested", njetcounterCT);
    pythia.init();
    // Remember position in vector of cross section estimates.
    int iNow = sizeLO-1-njetcounterCT;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weightNLO  = pythia.info.mergingWeightNLO();
      double evtweight  = pythia.info.weight();
      weightNLO        *= evtweight;
      // Do not print zero-weight events.
      if ( weightNLO == 0. ) continue; 

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = -1*xsecLO[iNow] / nAcceptLO[iNow];
      // Set hepmc event weight.
      hepmcevt->weights().push_back(weightNLO*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weightNLO*normhepmc;
      errorTotal += pow2(weightNLO*normhepmc);
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

    // Restart with ME of a reduced the number of jets
    if( njetcounterCT > 1 )
      njetcounterCT--;
    else
      break;

  }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  cout << endl << endl << endl;
  cout << "Shower subtractive events" << endl;

  pythia.settings.flag("Merging:doUNLOPSTree",false);
  pythia.settings.flag("Merging:doUNLOPSLoop",false);
  pythia.settings.flag("Merging:doUNLOPSSubt",false);
  pythia.settings.flag("Merging:doUNLOPSSubtNLO",true);
  pythia.settings.mode("Merging:nRecluster",1);

  nMaxCT = nMaxNLO;
  njetcounterCT = nMaxCT;
  iPathSubt= iPath + "_powheg";

  while(njetcounterCT >= 1){

    // From njet, choose LHE file
    stringstream in;
    in   << "_" << njetcounterCT << ".lhe";
    string LHEfile = iPathSubt + in.str();

    cout << endl << endl << endl
         << "Start subtractive treatment for " << njetcounterCT << " nlo jets"
         << endl;

    pythia.readString("Beams:frameType = 4"); 
    pythia.settings.word("Beams:LHEF", LHEfile);  
    pythia.settings.mode("Merging:nRequested", njetcounterCT);
    pythia.init();
    // Remember position in vector of cross section estimates.
    int iNow = sizeNLO-1-njetcounterCT;

    // Start generation loop
    for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if( pythia.info.atEndOfFile() ) break;
        else continue;
      }

      // Get event weight(s).
      double weightNLO  = pythia.info.mergingWeightNLO();
      double evtweight  = pythia.info.weight();
      weightNLO        *= evtweight;
      // Do not print zero-weight events.
      if ( weightNLO == 0. ) continue; 

      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      // Get correct cross section from previous estimate.
      double normhepmc = -1*xsecNLO[iNow] / nAcceptNLO[iNow];
      // powheg weighted events
      if( abs(strategyNLO[iNow]) == 4) 
        normhepmc = -1. / (1e9*nSelectedNLO[iNow]);
      // Set hepmc event weight.
      hepmcevt->weights().push_back(weightNLO*normhepmc);
      // Fill HepMC event.
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal += weightNLO*normhepmc;
      errorTotal += pow2(weightNLO*normhepmc);
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

    // Restart with ME of a reduced the number of jets
    if( njetcounterCT > 1 )
      njetcounterCT--;
    else
      break;

  }

  cout << "UNLOPS merged cross section: " << scientific << setprecision(8)
       << sigmaTotal << "  +-  " << sqrt(errorTotal) << " mb " << endl;
  cout << "NLO inclusive cross section: " << scientific << setprecision(8)
       << xsecNLO.back() << " mb "
       << endl;
  cout << "LO inclusive cross section:  " << scientific << setprecision(8)
       << xsecLO.back() << " mb "
       << endl;
  cout << endl << endl << endl;

  // Done
  return 0;

}
