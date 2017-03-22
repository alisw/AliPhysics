// main82.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do CKKW-L merging, 
// see the Matrix Element Merging page in the online manual. 

#include "Pythia.h"

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

// Class for user interaction with the merging

class MyMergingHooks : public MergingHooks {

private:

public:

  // Default constructor
  MyMergingHooks();
  // Destructor
  ~MyMergingHooks();

  // Functional definition of the merging scale
  virtual double tmsDefinition( const Event& event);

  // Helper function for tms definition
  double myKTdurham(const Particle& RadAfterBranch, 
           const Particle& EmtAfterBranch, int Type, double D );

};

//--------------------------------------------------------------------------

// Constructor
MyMergingHooks::MyMergingHooks() {}

// Destructor
MyMergingHooks::~MyMergingHooks() {}

//--------------------------------------------------------------------------

// Definition of the merging scale

double MyMergingHooks::tmsDefinition( const Event& event){

  // Cut only on QCD partons!
  // Count particle types
  int nFinalColoured = 0;
  int nFinalNow =0;
  for( int i=0; i < event.size(); ++i) {
    if(event[i].isFinal()){
      if(event[i].id() != 23 && abs(event[i].id()) != 24)
        nFinalNow++;
      if( event[i].colType() != 0)
        nFinalColoured++;
    }
  }

  // Use MergingHooks in-built functions to get information on the hard process
  int nLeptons = nHardOutLeptons();
  int nQuarks  = nHardOutPartons();
  int nResNow  = nResInCurrent();

  // Check if photons, electrons etc. have been produced. If so, do not veto
  if(nFinalNow - ( (nLeptons+nQuarks)/2 - nResNow)*2 != nFinalColoured){
    // Sometimes, Pythia detaches the decay products even though no
    // resonance was put into the LHE file, to catch this, add another
    // if statement
    if(nFinalNow != nFinalColoured) return 0.;
  }

  // Check that one parton has been produced. If not (e.g. in MPI), do not veto
  int nMPI = infoPtr->nMPI();
  if(nMPI > 1) return 0.;

  // Declare kT algorithm parameters
  double Dparam = 0.4;
  int kTtype = -1;
  // Declare final parton vector
  vector <int> FinalPartPos;
  FinalPartPos.clear();
  // Search event record for final state partons
  for (int i=0; i < event.size(); ++i)
    if(event[i].isFinal() && event[i].colType() != 0)
      FinalPartPos.push_back(i);

  // Find minimal Durham kT in event, using own function: Check
  // definition of separation
  int type = (event[3].colType() == 0 && event[4].colType() == 0) ? 1 : kTtype;
  // Find minimal kT
  double ktmin = event[0].e();
  for(int i=0; i < int(FinalPartPos.size()); ++i){
    double kt12  = ktmin;
    // Compute separation to the beam axis for hadronic collisions
    if(type == -1 || type == -2) {
      double temp = event[FinalPartPos[i]].pT();
      kt12 = min(kt12, temp);
    }
    // Compute separation to other final state jets
    for(int j=i+1; j < int(FinalPartPos.size()); ++j) {
      double temp = kTdurham( event[FinalPartPos[i]], event[FinalPartPos[j]],
                      type, Dparam);
      kt12 = min(kt12, temp);
    }
    // Keep the minimal Durham separation
    ktmin = min(ktmin,kt12);
  }

  // Return minimal Durham kT
  return ktmin;

}

//--------------------------------------------------------------------------

// Function to compute durham y separation from Particle input

double MyMergingHooks::myKTdurham(const Particle& RadAfterBranch,
  const Particle& EmtAfterBranch, int Type, double D ){

  // Declare return variable
  double ktdur;
  // Save 4-momenta of final state particles
  Vec4 jet1 = RadAfterBranch.p();
  Vec4 jet2 = EmtAfterBranch.p();

  if( Type == 1) {
    // Get angle between jets for e+e- collisions, make sure that
    // -1 <= cos(theta) <= 1
    double costh;
    if (jet1.pAbs()*jet2.pAbs() <=0.) costh = 1.;
    else {
      costh = costheta(jet1,jet2);
    }
    // Calculate kt durham separation between jets for e+e- collisions
    ktdur = 2.0*min( pow(jet1.e(),2) , (pow(jet2.e(),2)) )*(1.0 - costh);
  } else if( Type == -1 ){
    // Get delta_eta and cosh(Delta_eta) for hadronic collisions
    double eta1 = 0.5*log( (jet1.e() + jet1.pz()) / (jet1.e() - jet1.pz()) );
    double eta2 = 0.5*log( (jet2.e() + jet2.pz()) / (jet2.e() - jet2.pz()) );
    // Get delta_phi and cos(Delta_phi) for hadronic collisions
    double pt1 = sqrt( pow(jet1.px(),2) + pow(jet1.py(),2) );
    double pt2 = sqrt( pow(jet2.px(),2) + pow(jet2.py(),2) );  
    double cosdPhi = ( jet1.px()*jet2.px() + jet1.py()*jet2.py() ) / (pt1*pt2);
    double dPhi = acos( cosdPhi );
    // Calculate kT durham like fastjet
     ktdur = min( pow(pt1,2),pow(pt2,2) )
           * ( pow(eta1-eta2,2) + pow(dPhi,2) ) / pow(D,2);
  } else if( Type == -2 ){
    // Get delta_eta and cosh(Delta_eta) for hadronic collisions
    double eta1 = 0.5*log( (jet1.e() + jet1.pz()) / (jet1.e() - jet1.pz()) );
    double eta2 = 0.5*log( (jet2.e() + jet2.pz()) / (jet2.e() - jet2.pz()) );
     double coshdEta = cosh( eta1 - eta2 );
    // Get delta_phi and cos(Delta_phi) for hadronic collisions
    double pt1 = sqrt( pow(jet1.px(),2) + pow(jet1.py(),2) );
    double pt2 = sqrt( pow(jet2.px(),2) + pow(jet2.py(),2) );  
    double cosdPhi = ( jet1.px()*jet2.px() + jet1.py()*jet2.py() ) / (pt1*pt2);
    // Calculate kT durham separation "SHERPA-like"
     ktdur = 2.0*min( pow(pt1,2),pow(pt2,2) )
           * ( coshdEta - cosdPhi ) / pow(D,2);
  } else {
    ktdur = 0.0;
  }
  // Return kT
  return sqrt(ktdur);
}

//==========================================================================

// Example main programm to illustrate merging

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments. \n You are"
         << " expected to provide the arguments \n"
         << " 1. Input file for settings \n"
         << " 2. Full name of the input LHE file (with path) \n"
         << " 3. Path for output histogram files \n"
         << " Program stopped. " << endl;
    return 1;
  }

  Pythia pythia;

  // Input parameters:
  //  1. Input file for settings
  //  2. Path to input LHE file
  //  3. OUtput histogram path
  pythia.readFile(argv[1]);
  string iPath = string(argv[2]);
  string oPath = string(argv[3]);

  // Number of events
  int nEvent = pythia.mode("Main:numberOfEvents");

  // Construct user inut for merging
  MergingHooks* myMergingHooks = new MyMergingHooks();
  pythia.setMergingHooksPtr( myMergingHooks );

  // For ISR regularisation off
  pythia.settings.forceParm("SpaceShower:pT0Ref",0.);

  // Declare histograms
  Hist histPTFirst("pT of first jet",100,0.,100.);
  Hist histPTSecond("pT of second jet",100,0.,100.);

  // Read in ME configurations
  pythia.init(iPath,false);

  if(pythia.flag("Main:showChangedSettings")) {
    pythia.settings.listChanged();
  }

  // Start generation loop
  for( int iEvent=0; iEvent<nEvent; ++iEvent ){

    // Generate next event
    if( ! pythia.next()) continue;

    // Get CKKWL weight of current event
    double weight = pythia.info.mergingWeight();

    // Fill bins with CKKWL weight
    double pTfirst = pTfirstJet(pythia.event,1, 0.4);
    double pTsecnd = pTfirstJet(pythia.event,2, 0.4);
    histPTFirst.fill( pTfirst, weight);
    histPTSecond.fill( pTsecnd, weight);

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


  delete myMergingHooks;
  return 0;

  // Done
}
