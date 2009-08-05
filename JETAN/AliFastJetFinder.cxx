/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 

//---------------------------------------------------------------------
// FastJet v2.3.4 finder algorithm interface
// Last modification: Neutral cell energy included in the jet reconstruction
//
// Authors: Rafael.Diaz.Valdes@cern.ch
//          Magali.estienne@subatech.in2p3.fr (neutral part + bg subtraction option)
//
//---------------------------------------------------------------------


#include <Riostream.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TArrayF.h>
#include <TRandom.h>
#include <TClonesArray.h>

#include "AliFastJetFinder.h"
#include "AliFastJetHeaderV1.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJet.h"
#include "AliJetUnitArray.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
// get info on how fastjet was configured
#include "fastjet/config.h"

#ifdef ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#endif

#include<sstream>  // needed for internal io
#include<vector> 
#include <cmath> 

using namespace std;


ClassImp(AliFastJetFinder)


//____________________________________________________________________________

AliFastJetFinder::AliFastJetFinder():
    AliJetFinder()
{
  // Constructor
}

//____________________________________________________________________________

AliFastJetFinder::~AliFastJetFinder()
{
  // destructor
}

//______________________________________________________________________________
void AliFastJetFinder::FindJets()
{
  
  //pick up fastjet header
  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Bool_t debug  = header->GetDebug();     // debug option
  Bool_t bgMode = header->GetBGMode();    // choose to subtract BG or not
  Int_t fOpt    = fReader->GetReaderHeader()->GetDetector();

  // check if we are reading AOD jets
  TRefArray *refs = 0;
  Bool_t fromAod = !strcmp(fReader->ClassName(),"AliJetAODReader");
  if (fromAod) { refs = fReader->GetReferences(); }
  
  // RUN ALGORITHM  
  // read input particles -----------------------------
  vector<fastjet::PseudoJet> inputParticles;
  if(fOpt==0)
    {
      TClonesArray *lvArray = fReader->GetMomentumArray();
      if(lvArray == 0) { cout << "Could not get the momentum array" << endl; return; }
      Int_t nIn =  lvArray->GetEntries();
      if(nIn == 0) { if (debug) cout << "entries = 0 ; Event empty !!!" << endl ; return; }
      fJets->SetNinput(nIn) ; // number of input objects
      Float_t px,py,pz,en;
      // load input vectors
      for(Int_t i = 0; i < nIn; i++){ // loop for all input particles
	TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
	px = lv->Px();
	py = lv->Py();
	pz = lv->Pz();
	en = lv->Energy();
    
	fastjet::PseudoJet inputPart(px,py,pz,en); // create PseudoJet object
	inputPart.set_user_index(i); //label the particle into Fastjet algortihm
	inputParticles.push_back(inputPart);  // back of the inputParticles vector  
      } // end loop 
    }
  else {
    TClonesArray* fUnit = fReader->GetUnitArray();
    if(fUnit == 0) { cout << "Could not get the momentum array" << endl; return; }
    Int_t         nCandidate = fReader->GetNumCandidate();
    Int_t         nIn = fUnit->GetEntries();
    if(nIn == 0) { if (debug) cout << "entries = 0 ; Event empty !!!" << endl ; return; }
    fJets->SetNinput(nCandidate); // number of input objects // ME
    // Information extracted from fUnitArray
    // load input vectors and calculate total energy in array
    Float_t pt,eta,phi,theta,px,py,pz,en;
    Int_t ipart = 0;
    for(Int_t i=0; i<nIn; i++) 
      {
	AliJetUnitArray *uArray = (AliJetUnitArray*)fUnit->At(i);
	
	if(uArray->GetUnitEnergy()>0.){
	  
	  // It is not necessary anymore to cut on particle pt
	  pt    = uArray->GetUnitEnergy();
	  eta   = uArray->GetUnitEta();
	  phi   = uArray->GetUnitPhi();
	  theta = EtaToTheta(eta);
	  en    = (TMath::Abs(TMath::Sin(theta)) == 0) ? pt : pt/TMath::Abs(TMath::Sin(theta));
	  px    = TMath::Cos(phi)*pt;
	  py    = TMath::Sin(phi)*pt;
	  pz    = en*TMath::TanH(eta);
	  if(debug) cout << "pt: " << pt << ", eta: " << eta << ", phi: " << phi << ", en: " << en << ", px: " << px << ", py: " << py << ", pz: " << pz << endl;

	  fastjet::PseudoJet inputPart(px,py,pz,en); // create PseudoJet object
	  inputPart.set_user_index(ipart); //label the particle into Fastjet algortihm
	  inputParticles.push_back(inputPart);  // back of the inputParticles vector 
	  ipart++;
	}
      } // End loop on UnitArray 
  }
  
  // create an object that represents your choice of jet algorithm, and 
  // the associated parameters
  double rParam = header->GetRparam();
  fastjet::Strategy strategy = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();
  fastjet::JetAlgorithm algorithm = header->GetAlgorithm(); 
  fastjet::JetDefinition jet_def(algorithm, rParam, recombScheme, strategy);

  // create an object that specifies how we to define the area
  fastjet::AreaDefinition areaDef;
  double ghostEtamax = header->GetGhostEtaMax(); 
  double ghostArea   = header->GetGhostArea(); 
  int    activeAreaRepeats = header->GetActiveAreaRepeats(); 
  
  // now create the object that holds info about ghosts
  fastjet::GhostedAreaSpec ghost_spec(ghostEtamax, activeAreaRepeats, ghostArea);
  // and from that get an area definition
  fastjet::AreaType areaType = header->GetAreaType();
  areaDef = fastjet::AreaDefinition(areaType,ghost_spec);
  
  if(bgMode) // BG subtraction
    {
      //***************************** JETS FINDING AND EXTRACTION
      // run the jet clustering with the above jet definition
      fastjet::ClusterSequenceArea clust_seq(inputParticles, jet_def, areaDef);

      // save a comment in the header
      
      TString comment = "Running FastJet algorithm with the following setup. ";
      comment+= "Jet definition: ";
      comment+= TString(jet_def.description());
      comment+= ". Area definition: ";
      comment+= TString(areaDef.description());
      comment+= ". Strategy adopted by FastJet: ";
      comment+= TString(clust_seq.strategy_string());
      header->SetComment(comment);
      if(debug){
	cout << "--------------------------------------------------------" << endl;
	cout << comment << endl;
	cout << "--------------------------------------------------------" << endl;
      }
      //header->PrintParameters();
      
      
      // extract the inclusive jets with pt > ptmin, sorted by pt
      double ptmin = header->GetPtMin(); 
      vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(ptmin);
      
      //cout << "Number of unclustered particles: " << clust_seq.unclustered_particles().size() << endl;
      
      
      //subtract background // ===========================================
      // set the rapididty , phi range within which to study the background 
      double rapMax = header->GetRapMax(); 
      double rapMin = header->GetRapMin();
      double phiMax = header->GetPhiMax();
      double phiMin = header->GetPhiMin();
      fastjet::RangeDefinition range(rapMin, rapMax, phiMin, phiMax);
      
      // subtract background
      vector<fastjet::PseudoJet> subJets =  clust_seq.subtracted_jets(range,ptmin);  
      
      // print out
      //cout << "Printing inclusive sub jets with pt > "<< ptmin<<" GeV\n";
      //cout << "---------------------------------------\n";
      //cout << endl;
      //printf(" ijet   rap      phi        Pt         area  +-   err\n");
      
      // sort jets into increasing pt
      vector<fastjet::PseudoJet> jets = sorted_by_pt(subJets);  
      for (size_t j = 0; j < jets.size(); j++) { // loop for jets
	
	double area     = clust_seq.area(jets[j]);
	double areaError = clust_seq.area_error(jets[j]);
	
	printf("Jet found %5d %9.5f %8.5f %10.3f %8.3f +- %6.3f\n", (Int_t)j,jets[j].rap(),jets[j].phi(),jets[j].perp(), area, areaError);
	
	// go to write AOD  info
	AliAODJet aodjet (jets[j].px(), jets[j].py(), jets[j].pz(), jets[j].E());
	//cout << "Printing jet " << endl;
	if(debug) aodjet.Print("");
	//cout << "Adding jet ... " ;
	AddJet(aodjet);
	//cout << "added \n" << endl;
	
      }
    }
  else { // No BG subtraction

    fastjet::ClusterSequence clust_seq(inputParticles, jet_def); 

    // save a comment in the header
    
    TString comment = "Running FastJet algorithm with the following setup. ";
    comment+= "Jet definition: ";
    comment+= TString(jet_def.description());
    comment+= ". Strategy adopted by FastJet: ";
    comment+= TString(clust_seq.strategy_string());
    header->SetComment(comment);
    if(debug){
      cout << "--------------------------------------------------------" << endl;
      cout << comment << endl;
      cout << "--------------------------------------------------------" << endl;
    }
    //header->PrintParameters();
  
      // extract the inclusive jets with pt > ptmin, sorted by pt
    double ptmin = header->GetPtMin(); 
    vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(ptmin);
    
    //cout << "Number of unclustered particles: " << clust_seq.unclustered_particles().size() << endl;
    
    vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets); // Added by me
    for (size_t j = 0; j < jets.size(); j++) { // loop for jets     // Added by me
      
      printf("Jet found %5d %9.5f %8.5f %10.3f \n",(Int_t)j,jets[j].rap(),jets[j].phi(),jets[j].perp());
      
      // go to write AOD  info
      AliAODJet aodjet (jets[j].px(), jets[j].py(), jets[j].pz(), jets[j].E());
      //cout << "Printing jet " << endl;
      if(debug) aodjet.Print("");
      //cout << "Adding jet ... " ;
      AddJet(aodjet);
      //cout << "added \n" << endl;
      
    } // end loop for jets
  } 

}

//____________________________________________________________________________
void AliFastJetFinder::RunTest(const char* datafile)

{

   // This simple test run the kt algorithm for an ascii file testdata.dat
   // read input particles -----------------------------
  vector<fastjet::PseudoJet> inputParticles;
  Float_t px,py,pz,en;
  ifstream in;
  Int_t nlines = 0;
  // we assume a file basic.dat in the current directory
  // this file has 3 columns of float data
  in.open(datafile);
  while (1) {
      in >> px >> py >> pz >> en;
      if (!in.good()) break;
      //printf("px=%8f, py=%8f, pz=%8fn",px,py,pz);
      nlines++;
      inputParticles.push_back(fastjet::PseudoJet(px,py,pz,en)); 
   }
   //printf(" found %d pointsn",nlines);
   in.close();
   //////////////////////////////////////////////////
 
  // create an object that represents your choice of jet algorithm, and 
  // the associated parameters
  double rParam = 1.0;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recombScheme = fastjet::BIpt_scheme;
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, rParam, recombScheme, strategy);
  
  
 
  // create an object that specifies how we to define the area
  fastjet::AreaDefinition areaDef;
  double ghostEtamax = 7.0;
  double ghostArea    = 0.05;
  int    activeAreaRepeats = 1;
  

  // now create the object that holds info about ghosts
  fastjet::GhostedAreaSpec ghost_spec(ghostEtamax, activeAreaRepeats, ghostArea);
  // and from that get an area definition
  areaDef = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  

  // run the jet clustering with the above jet definition
  fastjet::ClusterSequenceArea clust_seq(inputParticles, jet_def, areaDef);
  
  
  // tell the user what was done
  cout << "--------------------------------------------------------" << endl;
  cout << "Jet definition was: " << jet_def.description() << endl;
  cout << "Area definition was: " << areaDef.description() << endl;
  cout << "Strategy adopted by FastJet was "<< clust_seq.strategy_string()<<endl<<endl;
  cout << "--------------------------------------------------------" << endl;
 
  
  // extract the inclusive jets with pt > 5 GeV, sorted by pt
  double ptmin = 5.0;
  vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(ptmin);
  
  cout << "Number of unclustered particles: " << clust_seq.unclustered_particles().size() << endl;
 
 
  //subtract background // ===========================================
  // set the rapididty range within which to study the background 
  double rapMax = ghostEtamax - rParam;
  fastjet::RangeDefinition range(rapMax);
  // subtract background
  vector<fastjet::PseudoJet> subJets =  clust_seq.subtracted_jets(range,ptmin);  
  
  // print them out //================================================
  cout << "Printing inclusive jets  after background subtraction \n";
  cout << "------------------------------------------------------\n";
  // sort jets into increasing pt
  vector<fastjet::PseudoJet> jets = sorted_by_pt(subJets);  

  printf(" ijet   rap      phi        Pt         area  +-   err\n");
  for (size_t j = 0; j < jets.size(); j++) {

    double area     = clust_seq.area(jets[j]);
    double areaError = clust_seq.area_error(jets[j]);

    printf("%5d %9.5f %8.5f %10.3f %8.3f +- %6.3f\n",(Int_t)j,jets[j].rap(),
	   jets[j].phi(),jets[j].perp(), area, areaError);
  }
  cout << endl;
  // ================================================================

  
 
}

//____________________________________________________________________________

void AliFastJetFinder::WriteJHeaderToFile()
{
  fHeader->Write();
}

//____________________________________________________________________________

Float_t  AliFastJetFinder::EtaToTheta(Float_t arg)
{
  //  return (180./TMath::Pi())*2.*atan(exp(-arg));
  return 2.*atan(exp(-arg));


}

//____________________________________________________________________________

void AliFastJetFinder::InitTask(TChain *tree)
{

  printf("Fast jet finder initialization ******************");
  fReader->CreateTasks(tree);

}
