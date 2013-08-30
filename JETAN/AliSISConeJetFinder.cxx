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
//
// Author: swensy.jangal@ires.in2p3.fr 
//
// ** 2011 magali.estienne@subatech.in2p3.fr &  alexandre.shabetai@cern.ch
// Modified accordingly to reader/finder splitting and new handling of neutral information (via FastJetInput) 
//---------------------------------------------------------------------

#include <Riostream.h>

#include "AliFastJetInput.h"
#include "AliFastJetBkg.h"
#include "AliAODJetEventBackground.h"
#include "AliAODJet.h"
#include "AliSISConeJetFinder.h"
#include "AliSISConeJetHeader.h"

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
// get info on how fastjet was configured
#include "fastjet/config.h"

#if defined(ENABLE_PLUGIN_SISCONE) || defined(FASTJET_ENABLE_PLUGIN_SISCONE)
#include "fastjet/SISConePlugin.hh"
#endif

#include<vector> 


using namespace std;

ClassImp(AliSISConeJetFinder)

////////////////////////////////////////////////////////////////////////

AliSISConeJetFinder::AliSISConeJetFinder():
  AliJetFinder(),
  fInputFJ(new AliFastJetInput()),
  fJetBkg(new  AliFastJetBkg())
{
  // Constructor
}

//____________________________________________________________________________
AliSISConeJetFinder::~AliSISConeJetFinder()
{
  // destructor
  delete  fInputFJ;
  delete  fJetBkg;

}

//______________________________________________________________________________
void AliSISConeJetFinder::FindJets()
{
  // run the SISCone Jet finder

  // Pick up siscone header  
  AliSISConeJetHeader *header = (AliSISConeJetHeader*)fHeader;
  Int_t debug                 = header->GetDebug();     // debug option
  Bool_t bgMode               = header->GetBGMode();// Here one choose to subtract BG or not

  // Read input particles 
  vector<fastjet::PseudoJet> inputParticles=fInputFJ->GetInputParticles();
  if(inputParticles.size()==0){
    if(debug>0) Printf("%s:%d No input particles found, skipping event",(char*)__FILE__,__LINE__);
    return;
  }

  //------------------- SISCONE PLUGIN CONFIGURATION ----------------------
  // Look for SISCone parameters in the header and definition of the plugin.  

  Double_t coneRadius                       = header->GetConeRadius();                 // cone radius
  Double_t overlapThreshold                 = header->GetOverlapThreshold();           // overlap parameter
  Int_t    nPassMax                         = header->GetNPassMax();                   // maximum number of passes
  Double_t ptProtoJetMin                    = header->GetPtProtojetMin();              // pT min of protojets
  Double_t caching                          = header->GetCaching();                    // do we record found cones for this set of data?
  // For bckg
  double rBkgParam                          = header->GetRparamBkg();
  fastjet::Strategy strategy                = header->GetStrategy();
  fastjet::RecombinationScheme recombScheme = header->GetRecombScheme();

  //  if (header->GetSplitMergeScale() == 0) fastjet::SISConePlugin::SplitMergeScale splitMergeScale = fastjet::SISConePlugin::SM_pttilde; // There's only one split merge scale
  //  Double_t splitMergeStoppingScale = header->GetSplitMergeStoppingScale(); // Additional cut on pt_tilde of protojets

  fastjet::JetDefinition::Plugin * plugin;
  plugin = new fastjet::SISConePlugin(coneRadius, overlapThreshold, nPassMax, ptProtoJetMin, caching);

  //------------------- CHOICE OF JET AREA ---------------------- 
  // Definition of jet areas for background subtraction
  // For more informations about jet areas see : The Catchment Area of Jets M. Cacciari, G. Salam and G. Soyez

  Double_t ghostEtamax        = header->GetGhostEtaMax();       // maximum eta in which a ghost can be generated
  Double_t ghostArea          = header->GetGhostArea();         // area of a ghost
  Int_t    activeAreaRepeats  = header->GetActiveAreaRepeats(); // do we repeat area calculation?
  Double_t gridScatter        = header->GetGridScatter();       // fractional random fluctuations of the position of the ghosts on the y-phi grid
  Double_t ktScatter          = header->GetKtScatter();         // fractional random fluctuations of the tranverse momentum of the ghosts on the y-phi grid
  Double_t meanGhostKt        = header->GetMeanGhostKt();       // average transverse momentum of the ghosts.

  Double_t areaTypeNumber  = header->GetAreaTypeNumber();       // the number determines jet area type 
  fastjet::AreaType areaType = fastjet::active_area;
  if (areaTypeNumber == 1)  areaType = fastjet::active_area;                 
  if (areaTypeNumber == 2)  areaType = fastjet::active_area_explicit_ghosts; 
  if (areaTypeNumber == 3)  areaType = fastjet::one_ghost_passive_area;      
  if (areaTypeNumber == 4)  areaType = fastjet::passive_area;                
  if (areaTypeNumber == 5)  areaType = fastjet::voronoi_area;                

  fastjet::AreaDefinition areaDef;

  if (areaTypeNumber < 5) 
    {
      fastjet::GhostedAreaSpec ghostSpec(ghostEtamax, activeAreaRepeats, ghostArea, gridScatter, ktScatter, meanGhostKt);
      areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
    }

  if (areaTypeNumber == 5)
    {
      Double_t effectiveRFact = header->GetEffectiveRFact();
      fastjet::VoronoiAreaSpec ghostSpec(effectiveRFact);
      areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
    }

  //------------------- JETS FINDING AND EXTRACTION ----------------------
  fastjet::ClusterSequenceArea clust_seq(inputParticles, plugin, areaDef);

  vector<fastjet::PseudoJet> jets;

  if (bgMode == 1)// BG subtraction mode
    {
      //------------------- CLUSTER JETS FINDING FOR RHO ESTIMATION ----------------------
      // run the jet clustering with the above jet definition
      fastjet::JetAlgorithm algorithmBkg = fastjet::kt_algorithm;
      Int_t algo = header->GetBGAlgorithm();
      if (algo == 0) algorithmBkg = fastjet::kt_algorithm;
      if (algo == 1) algorithmBkg = fastjet::cambridge_algorithm;
      fastjet::JetDefinition jetDefBkg(algorithmBkg, rBkgParam, recombScheme, strategy);
      fastjet::ClusterSequenceArea clust_seq_bkg(inputParticles, jetDefBkg, areaDef);
      
      // save a comment in the header
      TString comment = "Running Siscone algorithm with the following setup. ";
      comment+= "Jet definition: ";
      //      comment+= TString(plugin.description());
      comment+= "Jet bckg definition: ";
      comment+= TString(jetDefBkg.description());
      comment+= ". Area definition: ";
      comment+= TString(areaDef.description());
      comment+= ". Strategy adopted by FastJet and bkg: ";
      comment+= TString(clust_seq.strategy_string());
      header->SetComment(comment);
      if(debug>0){
	cout << "--------------------------------------------------------" << endl;
	cout << comment << endl;
	cout << "--------------------------------------------------------" << endl;
      }

      // Here we extract inclusive jets with pt > ptmin, sorted by pt 
      Double_t ptMin = header->GetMinJetPt(); 
      vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(); // ptMin removed
      //      vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusiveJets);

      //------------------- BACKGROUND SUBTRACTION ----------------------

      // Set the rapidity-azimuth range within which to study background
      Double_t rapMin = header->GetRapMin();
      Double_t rapMax = header->GetRapMax();
      Double_t phiMin = header->GetPhiMin();
      Double_t phiMax = header->GetPhiMax();
      fastjet::RangeDefinition range(rapMin, rapMax, phiMin, phiMax);

      // Extract rho and sigma
      Double_t rho = 0.;
      Double_t sigma = 0.;
      Double_t meanarea = 0.;
      Bool_t kUse4VectorArea = header->Use4VectorArea();
      vector<fastjet::PseudoJet> bkgJets = clust_seq_bkg.inclusive_jets();
      clust_seq_bkg.get_median_rho_and_sigma(bkgJets,range, kUse4VectorArea, rho, sigma, meanarea, false);

      // subtract background and extract jets bkg subtracted
      vector<fastjet::PseudoJet> subJets = clust_seq.subtracted_jets(rho,ptMin);

      // sort jets into increasing pt
      jets = sorted_by_pt(subJets);  

    }
  else // No BG subtraction
    {

      // save a comment in the header
      TString comment = "Running Siscone algorithm with the following setup. ";
      comment+= "Jet definition: ";
      //    comment+= TString(plugin.description());
      comment+= ". Strategy adopted by FastJet: ";
      comment+= TString(clust_seq.strategy_string());
      header->SetComment(comment);
      if(debug>0){
	cout << "--------------------------------------------------------" << endl;
	cout << comment << endl;
	cout << "--------------------------------------------------------" << endl;
      }
      //header->PrintParameters();

      // Here we extract inclusive jets with pt > ptmin, sorted by pt 
      Double_t ptMin = header->GetMinJetPt(); 
      vector<fastjet::PseudoJet> inclusiveJets = clust_seq.inclusive_jets(ptMin);
      jets = sorted_by_pt(inclusiveJets);

    }

  //------------------- JET AND TRACK STORAGE ----------------------
  for (size_t j = 0; j < jets.size(); j++) { // loop for jets

    double area      = clust_seq.area(jets[j]);
    double areaError = clust_seq.area_error(jets[j]);

    if(debug>0) printf("Jet found %5d %9.5f %8.5f %10.3f %8.3f +- %6.3f\n", (Int_t)j,jets[j].rap(),jets[j].phi(),jets[j].perp(), area, areaError);

    vector<fastjet::PseudoJet> constituents = clust_seq.constituents(jets[j]);
    int nCon= constituents.size();
    TArrayI ind(nCon);
      
    if ((jets[j].eta() > (header->GetJetEtaMax())) ||
        (jets[j].eta() < (header->GetJetEtaMin())) ||
        (jets[j].phi() > (header->GetJetPhiMax())) ||
        (jets[j].phi() < (header->GetJetPhiMin())) ||
        (jets[j].perp() < header->GetMinJetPt())) continue; // acceptance eta range and etmin

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
  
  delete plugin;

}
 
//____________________________________________________________________________
void AliSISConeJetFinder::WriteJHeaderToFile() const
{
  // Write Jet Header To File(
  fHeader->Write();
}

//____________________________________________________________________________
Bool_t AliSISConeJetFinder::ProcessEvent()
{
  // Process one event
  // Charged only or charged+neutral jets

  fInputFJ->SetHeader(fHeader);
  fInputFJ->SetCalTrkEvent(GetCalTrkEvent());
  fInputFJ->FillInput();
  
  // Jets
  FindJets();

  // Background  
  if( fAODEvBkg){
    fJetBkg->SetHeader(fHeader);
    Double_t sigma1 = 0,meanarea1= 0,sigma2 = 0,meanarea2 = 0;
    Double_t bkg1 = 0,bkg2 = 0;

    fJetBkg->SetFastJetInput(fInputFJ);
    fJetBkg->BkgFastJetb(bkg1,sigma1,meanarea1);
    fJetBkg->BkgFastJetWoHardest(bkg2,sigma2,meanarea2);
    fAODEvBkg->SetBackground(0,bkg1,sigma1,meanarea1);
    fAODEvBkg->SetBackground(1,bkg2,sigma2,meanarea2);
  }

  Reset();
  return kTRUE;

}
