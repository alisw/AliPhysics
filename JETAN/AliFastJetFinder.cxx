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
 
/* $Id$ */

//---------------------------------------------------------------------
// FastJet v2.3.4 finder algorithm interface
// Last modification: Neutral cell energy included in the jet reconstruction
//
// Authors: Rafael.Diaz.Valdes@cern.ch
//          Magali.estienne@subatech.in2p3.fr (neutral part + bg subtraction option)
//
// ** 2011 magali.estienne@subatech.in2p3.fr &  alexandre.shabetai@cern.ch
// new implementation of background subtraction
// allowing to subtract bkg using a different algo than the one used for signal jets
//---------------------------------------------------------------------

#include <Riostream.h>

#include "AliFastJetFinder.h"
#include "AliFastJetHeaderV1.h"
#include "AliFastJetInput.h"
#include "AliFastJetBkg.h"
#include "AliAODJetEventBackground.h"
#include "AliAODJet.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"

#include<vector> 

using namespace std;

ClassImp(AliFastJetFinder)

////////////////////////////////////////////////////////////////////////

AliFastJetFinder::AliFastJetFinder():
  AliJetFinder(),
  fInputFJ(new AliFastJetInput()),
  fJetBkg(new  AliFastJetBkg())
{
  // Constructor
}

//____________________________________________________________________________
AliFastJetFinder::~AliFastJetFinder()
{
  // destructor
  delete  fInputFJ;
  delete  fJetBkg;

}

//______________________________________________________________________________
void AliFastJetFinder::FindJets()
{
  // runs a FASTJET based algo

  //pick up fastjet header
  AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;
  Int_t  debug  = header->GetDebug();     // debug option
  Bool_t bgMode = header->GetBGMode();    // choose to subtract BG or not
  if(debug>0) cout<<"----------in AliFastJetFinder::FindJets() ------------------"<<endl;

  // RUN ALGORITHM  
  // read input particles -----------------------------

  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  if(inputParticles.size()==0){
    if(debug>0) Printf("%s:%d No input particles found, skipping event",(char*)__FILE__,__LINE__);
    return;
  }

  // create an object that represents your choice of jet algorithm, and 
  // the associated parameters
  double rParam = header->GetRparam();
  double rBkgParam = header->GetRparamBkg();
  fastjet::Strategy strategy = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();
  fastjet::JetAlgorithm algorithm = header->GetAlgorithm(); 
  fastjet::JetDefinition jetDef(algorithm, rParam, recombScheme, strategy);

  // create an object that specifies how we to define the area
  fastjet::AreaDefinition areaDef;
  double ghostEtamax       = header->GetGhostEtaMax(); 
  double ghostArea         = header->GetGhostArea(); 
  int    activeAreaRepeats = header->GetActiveAreaRepeats(); 
  
  // now create the object that holds info about ghosts
  fastjet::GhostedAreaSpec ghostSpec(ghostEtamax, activeAreaRepeats, ghostArea);
  // and from that get an area definition
  fastjet::AreaType areaType = header->GetAreaType();
  areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
  
  //***************************** JETS FINDING
  // run the jet clustering with the above jet definition
  fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef, areaDef);

  vector<fastjet::PseudoJet> jets;

  if(bgMode) // Do BG subtraction directly with the same algorithm (cambridge or kt) for jet signal and background
    {
      //***************************** JETS FINDING FOR RHO ESTIMATION
      // run the jet clustering with the above jet definition
      fastjet::JetAlgorithm algorithmBkg = header->GetBGAlgorithm();
      fastjet::JetDefinition jetDefBkg(algorithmBkg, rBkgParam, recombScheme, strategy);
      fastjet::ClusterSequenceArea clust_seq_bkg(inputParticles, jetDefBkg, areaDef);

      // save a comment in the header
      TString comment = "Running FastJet algorithm with the following setup. ";
      comment+= "Jet definition: ";
      comment+= TString(jetDef.description());
      comment+= "Jet bckg definition: ";
      comment+= TString(jetDefBkg.description());
      comment+= ". Area definition: ";
      comment+= TString(areaDef.description());
      comment+= ". Strategy adopted by FastJet: ";
      comment+= TString(clust_seq.strategy_string());
      header->SetComment(comment);
      if(debug>0){
	cout << "--------------------------------------------------------" << endl;
	cout << comment << endl;
	cout << "--------------------------------------------------------" << endl;
      }
      
      // extract the inclusive jets sorted by pt
      double ptmin = header->GetPtMin(); 
      vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets();
      
      //subtract background // ===========================================
      // set the rapididty , phi range within which to study the background 
      double rapMax = header->GetRapMax(); 
      double rapMin = header->GetRapMin();
      double phiMax = header->GetPhiMax();
      double phiMin = header->GetPhiMin();
      fastjet::RangeDefinition range(rapMin, rapMax, phiMin, phiMax);
     
      // Extract rho and sigma
      Double_t rho = 0.;
      Double_t sigma = 0.;
      Double_t meanarea = 0.;
      Bool_t kUse4VectorArea = header->Use4VectorArea();
      vector<fastjet::PseudoJet> bkgJets = clust_seq_bkg.inclusive_jets();
      clust_seq_bkg.get_median_rho_and_sigma(bkgJets,range, kUse4VectorArea, rho, sigma, meanarea, false);

      // subtract background and extract jets bkg subtracted
      vector<fastjet::PseudoJet> subJets = clust_seq.subtracted_jets(rho,ptmin);
      
      // sort jets into increasing pt
      jets = sorted_by_pt(subJets);  

    }

  else { // No BG subtraction!!!!!!!! Default header is bgmode=0.
    
    // save a comment in the header
    TString comment = "Running FastJet algorithm with the following setup. ";
    comment+= "Jet definition: ";
    comment+= TString(jetDef.description());
    comment+= ". Strategy adopted by FastJet: ";
    comment+= TString(clust_seq.strategy_string());
    header->SetComment(comment);
    if(debug>0){
      cout << "--------------------------------------------------------" << endl;
      cout << comment << endl;
      cout << "--------------------------------------------------------" << endl;
    }
  
    // extract the inclusive jets with pt > ptmin, sorted by pt
    double ptmin = header->GetPtMin(); 
    vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(ptmin);
    
    jets = sorted_by_pt(inclusiveJets); 
  
  }
   
  for (size_t j = 0; j < jets.size(); j++) { // loop for jets     

    double area      = clust_seq.area(jets[j]);
    double areaError = clust_seq.area_error(jets[j]);

    if(debug>0) printf("Jet found %5d %9.5f %8.5f %10.3f \n",(Int_t)j,jets[j].rap(),jets[j].phi(),jets[j].perp());

    vector<fastjet::PseudoJet> constituents = clust_seq.constituents(jets[j]);
    int nCon= constituents.size();
    TArrayI ind(nCon);
      
    if ((jets[j].eta() > (header->GetJetEtaMax())) ||
	(jets[j].eta() < (header->GetJetEtaMin())) ||
	(jets[j].phi() > (header->GetJetPhiMax())) ||
	(jets[j].phi() < (header->GetJetPhiMin())) ||
	(jets[j].perp() < header->GetPtMin())) continue; // acceptance eta range and etmin

    // go to write AOD  info
    AliAODJet aodjet (jets[j].px(), jets[j].py(), jets[j].pz(), jets[j].E());
    aodjet.SetEffArea(area,areaError);
    //cout << "Printing jet " << endl;
    if(debug>0) aodjet.Print("");
      
    for (int i=0; i < nCon; i++)
      {
	fastjet::PseudoJet mPart=constituents[i];
	ind[i]=mPart.user_index();

	// Jet constituents (charged tracks) added to the AliAODJet
	AliJetCalTrkEvent* calEvt  = GetCalTrkEvent();
	for(Int_t itrack=0; itrack<calEvt->GetNCalTrkTracks(); itrack++)
	  {
	    if(itrack==ind[i])
	      {
		TObject *track = calEvt->GetCalTrkTrack(itrack)->GetTrackObject();
		aodjet.AddTrack(track);
	      }
	  }
      } // End loop on Constituents

    AddJet(aodjet);
           
  } // end loop for jets

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
  fastjet::JetDefinition jetDef(fastjet::kt_algorithm, rParam, recombScheme, strategy);
  
  // create an object that specifies how we to define the area
  fastjet::AreaDefinition areaDef;
  double ghostEtamax = 7.0;
  double ghostArea    = 0.05;
  int    activeAreaRepeats = 1;
  
  // now create the object that holds info about ghosts
  fastjet::GhostedAreaSpec ghostSpec(ghostEtamax, activeAreaRepeats, ghostArea);
  // and from that get an area definition
  areaDef = fastjet::AreaDefinition(fastjet::active_area,ghostSpec);

  // run the jet clustering with the above jet definition
  fastjet::ClusterSequenceArea clust_seq(inputParticles, jetDef, areaDef);
  
  // tell the user what was done
  cout << "--------------------------------------------------------" << endl;
  cout << "Jet definition was: " << jetDef.description() << endl;
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
    double area      = clust_seq.area(jets[j]);
    double areaError = clust_seq.area_error(jets[j]);
    printf("%5d %9.5f %8.5f %10.3f %8.3f +- %6.3f\n",(Int_t)j,jets[j].rap(),
	   jets[j].phi(),jets[j].perp(), area, areaError);
  }
  cout << endl;
  // ================================================================

}

//____________________________________________________________________________
void AliFastJetFinder::WriteJHeaderToFile() const
{
  // Write Jet Header to file
  fHeader->Write();

}

//____________________________________________________________________________
Bool_t AliFastJetFinder::ProcessEvent()
{
  // Process one event
  // Charged only or charged+neutral jets

  fInputFJ->SetHeader(fHeader);
  fInputFJ->SetCalTrkEvent(GetCalTrkEvent());
  fInputFJ->FillInput();
  
  // Find Jets
  FindJets(); 

  // Background
  if(fAODEvBkg){
    AliFastJetHeaderV1 *header = (AliFastJetHeaderV1*)fHeader;

    fJetBkg->SetHeader(fHeader);
    fJetBkg->SetFastJetInput(fInputFJ);
    
     Int_t count = 0;  
     if(header->GetBkgFastJetb()){
        Double_t sigma1 = 0 , meanarea1= 0, bkg1 = 0;
        fJetBkg->BkgFastJetb(bkg1,sigma1,meanarea1);
        fAODEvBkg->SetBackground(count,bkg1,sigma1,meanarea1);
        count++;
     }

     if(header->GetBkgFastJetWoHardest()){
        Double_t sigma2 = 0 , meanarea2 = 0, bkg2 = 0;
        fJetBkg->BkgFastJetWoHardest(bkg2,sigma2,meanarea2);
        fAODEvBkg->SetBackground(count,bkg2,sigma2,meanarea2);
        count++;
     }

  }  

  Reset();  
  return kTRUE;

}
