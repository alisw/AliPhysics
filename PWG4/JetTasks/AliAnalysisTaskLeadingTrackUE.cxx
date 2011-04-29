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

/* $Id:$ */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TTree.h>
#include <TH2F.h>
#include <TRandom.h>
//#include <TVector3.h>

#include "AliAnalysisTaskLeadingTrackUE.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliUEHistograms.h"
#include "AliUEHist.h"

#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliVParticle.h"


////////////////////////////////////////////////////////////////////////
//
// Analysis class for Underlying Event studies w.r.t. leading track
//
// Look for correlations on the tranverse regions w.r.t
// the leading track in the event
//
// This class needs input AODs.
// The output is a list of analysis-specific containers.
//
// The AOD can be either connected to the InputEventHandler  
// for a chain of AOD files 
// or 
// to the OutputEventHandler
// for a chain of ESD files,
// in this case the class should be in the train after the jet-finder
//
//    Authors:
//    Arian Abrahantes Quintana 
//    Jan Fiete Grosse-Oetringhaus
//    Ernesto Lopez Torres
//    Sara Vallero
// 
////////////////////////////////////////////////////////////////////////


ClassImp( AliAnalysisTaskLeadingTrackUE )

// Define global pointer
AliAnalysisTaskLeadingTrackUE* AliAnalysisTaskLeadingTrackUE::fgTaskLeadingTrackUE=NULL;

//____________________________________________________________________
AliAnalysisTaskLeadingTrackUE:: AliAnalysisTaskLeadingTrackUE(const char* name):
AliAnalysisTask(name,""),
// general configuration
fDebug(0),
fMode(0),
fReduceMemoryFootprint(kFALSE),
// pointers to UE classes
fAnalyseUE(0x0),
fHistosUE(0x0),
fkTrackingEfficiency(0x0),
// handlers and events
fAOD(0x0),           
fArrayMC(0x0),
fInputHandler(0x0),
fMcEvent(0x0),
fMcHandler(0x0),
// histogram settings
fListOfHistos(0x0), 
fBinsPtInHist(30),     
fMinJetPtInHist(0.),
fMaxJetPtInHist(300.), 
// event QA
fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default) 
fZVertex(10.),
// track cuts
fTrackEtaCut(0.8),
fLeadingTrackEtaCut(0.8),
fFilterBit(0xFF),
fSelectBit(0),
fUseChargeHadrons(kFALSE),
//For MC weighting
fAvgTrials(1)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());

}

AliAnalysisTaskLeadingTrackUE::~AliAnalysisTaskLeadingTrackUE() 
{ 
  // destructor
  
  if (fListOfHistos  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) 
    delete fListOfHistos;
}

/************** INTERFACE METHODS *****************************/

//______________________________________________________________
Bool_t AliAnalysisTaskLeadingTrackUE::Notify()
{
  
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root.
  // This will be used when merging different MC samples.
  // (Copied from AliAnalysisTaskJFSystematics)
  
  fAvgTrials = 1;
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t trials  = 1;
  if(tree){
  	TFile *curfile = tree->GetCurrentFile();
  	if (!curfile) {
  		Error("Notify","No current file");
  		return kFALSE;
  		}
  	
	AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,trials); 
        
	//TO-DO
        //fHistosUE->GetXsec()->Fill("<#sigma>",xsection);

  	// construct average trials
  	Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
  	if(trials>=nEntries && nEntries>0.)fAvgTrials = trials/nEntries;
  	}
  
  return kTRUE;
}

//____________________________________________________________________
void AliAnalysisTaskLeadingTrackUE::ConnectInputData(Option_t* /*option*/)
{
  
  // Connect the input data  
  if (fDebug > 1) AliInfo("ConnectInputData() ");
  
  // Since AODs can either be connected to the InputEventHandler
  // or to the OutputEventHandler ( the AOD is created by a previus task in the train )
  // we need to get the pointer to the AODEvent correctly.
  
  // Delta AODs are also accepted.
  
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) { // input AOD
  	fAOD = ((AliAODInputHandler*)handler)->GetEvent();
	if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODInputHandler");
  } else {  //output AOD
  	handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
  	if( handler && handler->InheritsFrom("AliAODHandler") ) {
  		fAOD = ((AliAODHandler*)handler)->GetAOD();
  		if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODHandler");
 	} else {  // no AOD
		AliFatal("I can't get any AOD Event Handler");
  		return;
  		}
  	}	
  
  // Initialize common pointers
  Initialize();
   
}

//____________________________________________________________________
void  AliAnalysisTaskLeadingTrackUE::CreateOutputObjects()
{
  // Create the output container
  
  if (fDebug > 1) AliInfo("CreateOutputObjects()");
   
  // Initialize class with main algorithms, event and track selection. 
  fAnalyseUE = new AliAnalyseLeadingTrackUE();
  fAnalyseUE->SetParticleSelectionCriteria(fFilterBit, fUseChargeHadrons, fLeadingTrackEtaCut);
  fAnalyseUE->SetDebug(fDebug); 

  // Initialize output list of containers
  if (fListOfHistos != NULL){
	delete fListOfHistos;
        fListOfHistos = NULL;
  	}
  if (!fListOfHistos){
  	fListOfHistos = new TList();
  	fListOfHistos->SetOwner(kTRUE); 
  	}

  // Initialize class to handle histograms 
  fHistosUE = new AliUEHistograms;
  
  // add histograms to list
  fListOfHistos->Add(fHistosUE);
  
  //fListOfHistos->Add(new TH2F("multVsLeadStep5", ";multiplicity;leading pT", 100, -0.5, 99.5, 20, 0, 10));
  
  // Add task configuration to output list 
  AddSettingsTree();
  

  PostData(0,fListOfHistos);
}

//____________________________________________________________________
void  AliAnalysisTaskLeadingTrackUE::Exec(Option_t */*option*/)
{
  // array of MC particles
  if (fMcHandler){
  	fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  	if (!fArrayMC)AliFatal("No array of MC particles found !!!");
	}

  // Get number of trials from MC header
  Float_t nTrials = 1;
  if (fMcHandler) {
  	fMcEvent = fMcHandler->MCEvent();
  	if (fMcEvent) {
                // TO-DO: extend to PHOJET
  		AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(fMcEvent);
  		if(pythiaGenHeader){
  			nTrials = pythiaGenHeader->Trials();
  			}
  		}
  	}

  // TO-DO	
  //fHistosUE->GetTrials()->Fill("#sum{ntrials}",fAvgTrials);
  
  // Analyse the event
  if (fMode) AnalyseCorrectionMode();
  else AnalyseDataMode();

  PostData(0,fListOfHistos);
}

//____________________________________________________________________
void  AliAnalysisTaskLeadingTrackUE::Terminate(Option_t */*option*/)
{
  
  // Terminate analysis
  if( fDebug > 1 ) AliInfo("End analysis");

  if (!gROOT->IsBatch()){
  	fListOfHistos = dynamic_cast<TList*> (GetOutputData(0));
  	if (!fListOfHistos){
		AliError("Histogram List is not available");
		return;
  	}else{
	// Draw something
	}
  } else {
        AliInfo(" Batch mode, not histograms will be shown...");
  	}

}


/******************** ANALYSIS METHODS *****************************/

//____________________________________________________________________
void  AliAnalysisTaskLeadingTrackUE::AddSettingsTree()
{
  //Write settings to output list
  TTree *settingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fSelectBit", &fSelectBit,"EventSelectionBit/I");
  settingsTree->Branch("fLeadingTrackEtaCut", &fLeadingTrackEtaCut, "LeadingTrackEtaCut/D");
  settingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  settingsTree->Branch("fkTrackingEfficiency", "TH1D", &fkTrackingEfficiency);
  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
}  

//____________________________________________________________________
void  AliAnalysisTaskLeadingTrackUE::AnalyseCorrectionMode()
{
  // Run the analysis on MC to get the correction maps
  //
  // if fReduceMemoryFootprint step 3,4,5,7,9 are not filled

  PostData(0,fListOfHistos);
  
  if ( fDebug > 3 ) AliInfo( " Processing event in Corrections mode ..." );
  
  //PROCESS TYPE (ND,SD,DD)
  AliAnalysisHelperJetTasks::MCProcessType eventId = AliAnalysisHelperJetTasks::kInvalidProcess;
  AliGenEventHeader* genHeader = fMcEvent->GenEventHeader();
  eventId = AliAnalysisHelperJetTasks::GetPythiaEventProcessType(genHeader,kFALSE);
  if (eventId<0)
    eventId = AliAnalysisHelperJetTasks::GetDPMjetEventProcessType(genHeader,kFALSE);
  if (eventId<0 && fDebug>1)
    AliInfo("No Pythia or Phojet Header retrived!");
    
  Int_t fillId=-1;
  if (eventId == AliAnalysisHelperJetTasks::kND)fillId = 0; 
  if (eventId == AliAnalysisHelperJetTasks::kSD)fillId = 1; 
  if (eventId == AliAnalysisHelperJetTasks::kDD)fillId = 2; 
  
  // count all events
  fHistosUE->FillEvent(fillId, -1);
  
  // Only consider MC events within the vtx-z region used also as cut on the reconstructed vertex
  if (!fAnalyseUE->VertexSelection(fMcEvent, 0, fZVertex)) 
    return;
  
  // Get MC-true leading particle (but do not cut out events!)
  TObjArray *ltMC = (TObjArray*)fAnalyseUE->FindLeadingObjects(fArrayMC);
  AliVParticle* leadingMC = 0;
  if (ltMC)
    leadingMC = (AliVParticle*) ltMC->At(0);
    
  // it can happen that there is no MC leading particle in the acceptance required (|eta|<0.8)
  // and we do not want to base the event slection on MC information
  
  // Sort MC-true charged particles
  // as output you get an array of 3 lists of  particles belonging to different regions:
  // - at 0: towards
  // - at 1: away
  // - at 2: transverse MIN
  // - at 3: transverse MAX
  TObjArray *regionSortedParticlesMC = (TObjArray*)fAnalyseUE->SortRegions(leadingMC, fArrayMC, 0x0); 
  TObjArray *regionsMinMaxMC = (TObjArray*)fAnalyseUE->GetMinMaxRegion((TList*)regionSortedParticlesMC->At(2),(TList*)regionSortedParticlesMC->At(3));
  // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
  // (MC-true leading particle and MC-true all particles)
  // STEP 0
  fHistosUE->Fill(fillId,0,AliUEHist::kCFStepAll,leadingMC,(TList*)regionSortedParticlesMC->At(0),(TList*)regionSortedParticlesMC->At(1),(TList*)regionsMinMaxMC->At(0),(TList*)regionsMinMaxMC->At(1));
  
  // Trigger selection ************************************************
  if (fAnalyseUE->TriggerSelection(fInputHandler))
  {
    // PILEUP-CUT 
    Bool_t select = kFALSE;
    if (fSelectBit) select = AliAnalysisHelperJetTasks::TestSelectInfo(fSelectBit);
    if (select)
      fHistosUE->FillEvent(fillId, -2);
    else
    {
    
	    // Count events that pass AliPhysicsSelection
    
  	    // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
    	    // (MC-true leading particle and MC-true all particles)
    	    // STEP 1
	    fHistosUE->Fill(fillId,0,AliUEHist::kCFStepTriggered,leadingMC,(TList*)regionSortedParticlesMC->At(0),(TList*)regionSortedParticlesMC->At(1),(TList*)regionsMinMaxMC->At(0),(TList*)regionsMinMaxMC->At(1));
  
	    // count number of MC tracks above 150 MeV/c
	    Int_t nMCTracks = 0; 
	    if (leadingMC && leadingMC->Pt() > 0.15)
	      nMCTracks++;
	    for (Int_t i=0; i<4; i++)
	      for (Int_t j=0; j<((TList*)regionSortedParticlesMC->At(i))->GetEntries(); j++)
	        if (((AliVParticle*) ((TList*)regionSortedParticlesMC->At(i))->At(j))->Pt() > 0.15)
	          nMCTracks++;
      
	    //((TH2F*)fListOfHistos->FindObject("multVsLeadStep5"))->Fill(nMCTracks, leadingMC->Pt());
    
	    // Vertex selection *************************************************
    	    if (fAnalyseUE->VertexSelection(fAOD, fnTracksVertex, fZVertex))
	    {
	      // Count events that pass Vertex selection
            
	      // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
	      // (MC-true leading particle and MC-true all particles)
	      // STEP 2
	      fHistosUE->Fill(fillId,0,AliUEHist::kCFStepVertex,leadingMC,(TList*)regionSortedParticlesMC->At(0),(TList*)regionSortedParticlesMC->At(1),(TList*)regionsMinMaxMC->At(0),(TList*)regionsMinMaxMC->At(1));
    
              // fill here for tracking efficiency
              // loop over particle species
	      for (Int_t particleSpecies = 0; particleSpecies < 4; particleSpecies++)
              {
              	TObjArray* primMCParticles = fAnalyseUE->GetAcceptedParticles(fArrayMC, 0x0, kTRUE, particleSpecies);
              	TObjArray* primRecoTracksMatched = fAnalyseUE->GetAcceptedParticles(fAOD, fArrayMC, kTRUE, particleSpecies);
              	TObjArray* allRecoTracksMatched = fAnalyseUE->GetAcceptedParticles(fAOD, fArrayMC, kFALSE, particleSpecies);
              
              	fHistosUE->FillTrackingEfficiency(primMCParticles, primRecoTracksMatched, allRecoTracksMatched, particleSpecies);
              
              	delete primMCParticles;
              	delete primRecoTracksMatched;
              	delete allRecoTracksMatched;
	      }
      
	      // Get Reconstructed leading particle *******************************
	      TObjArray *ltRECO = fAnalyseUE->FindLeadingObjects(fAOD);
	      if (ltRECO)
	      {
	        // Count events where a reconstructed track was found in |eta|<0.8
	        // the pT cut will be set when projecting output containers
	        // for leading particle correlation plots
	        if (leadingMC) {
	              fHistosUE->Fill(leadingMC, (AliVParticle*)ltRECO->At(0));
	              }
              
	        // If there is no MC leading track the container is not filled, so the number of entries in the container might be different  
	        // from the number of events after the selection, since the selection is based on RECO tracks
	        // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
	        // (MC-true leading particle and MC-true all particles)
	        // STEP 3
	        if (!fReduceMemoryFootprint)
	          fHistosUE->Fill(fillId,0,AliUEHist::kCFStepAnaTopology,leadingMC,(TList*)regionSortedParticlesMC->At(0),(TList*)regionSortedParticlesMC->At(1),(TList*)regionsMinMaxMC->At(0),(TList*)regionsMinMaxMC->At(1));
      
	        //Sort RECO particles w.r.t. MC-leading and return matched (primary) MC particle 
	        // (you cannot sort tracks w.r.t. RECO-leading and plot it vs. MC-leading ...)
	        TObjArray *regionSortedParticlesRECOLTMC = (TObjArray*)fAnalyseUE->SortRegions(leadingMC, fAOD, fArrayMC, kTRUE); 
	        TObjArray *regionsMinMaxRECOLTMC = (TObjArray*)fAnalyseUE->GetMinMaxRegion((TList*)regionSortedParticlesRECOLTMC->At(2),(TList*)regionSortedParticlesRECOLTMC->At(3));
	        // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
	        // (MC leading particle and RECO-matched (quantities from MC particle)  all particles)
	        // STEP 4
	        //if (!fReduceMemoryFootprint)
  	          fHistosUE->Fill(fillId,0,AliUEHist::kCFStepTrackedOnlyPrim,leadingMC,(TList*)regionSortedParticlesRECOLTMC->At(0),(TList*)regionSortedParticlesRECOLTMC->At(1),(TList*)regionsMinMaxRECOLTMC->At(0),(TList*)regionsMinMaxRECOLTMC->At(1));
	        // comparing this step with step 3 (for all-tracks observables) you get the tracking efficiency
        
	        //Sort RECO particles w.r.t. MC-leading and return matched (primary+secondary) MC particle 
	        // (you cannot sort tracks w.r.t. RECO-leading and plot it vs. MC-leading ...)
	        TObjArray *regionSortedParticlesRECOLTMC2 = (TObjArray*)fAnalyseUE->SortRegions(leadingMC, fAOD, fArrayMC,kFALSE); 
	        TObjArray *regionsMinMaxRECOLTMC2 = (TObjArray*)fAnalyseUE->GetMinMaxRegion((TList*)regionSortedParticlesRECOLTMC2->At(2),(TList*)regionSortedParticlesRECOLTMC2->At(3));
        
	        // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
	        // (MC leading particle and RECO-matched (quantities from MC particle)  all particles)
	        // STEP 5
	        //if (!fReduceMemoryFootprint)
                  fHistosUE->Fill(fillId,0,AliUEHist::kCFStepTracked,leadingMC,(TList*)regionSortedParticlesRECOLTMC2->At(0),(TList*)regionSortedParticlesRECOLTMC2->At(1),(TList*)regionsMinMaxRECOLTMC2->At(0),(TList*)regionsMinMaxRECOLTMC2->At(1));
	        // comparing this step with step 3 (for all-tracks observables) you get the tracking efficiency
          
	        // SWITCH TO RECONSTRUCTED TRACKS  ************************************
	        // The next steps correspond to track selections
	        // Sort RECO particles w.r.t. RECO-leading and return RECO particle  
	        TObjArray *regionSortedParticlesRECO = (TObjArray*)fAnalyseUE->SortRegions((AliVParticle*)ltRECO->At(0), fAOD,0); 
	        TObjArray *regionsMinMaxRECO = (TObjArray*)fAnalyseUE->GetMinMaxRegion((TList*)regionSortedParticlesRECO->At(2),(TList*)regionSortedParticlesRECO->At(3));
	        // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
	        // (RECO leading particle and RECO  all particles)
	        // STEP 6
	        fHistosUE->Fill(fillId,0,AliUEHist::kCFStepReconstructed,(AliVParticle*)ltRECO->At(0),(TList*)regionSortedParticlesRECO->At(0),(TList*)regionSortedParticlesRECO->At(1),(TList*)regionsMinMaxRECO->At(0),(TList*)regionsMinMaxRECO->At(1));
        
	        // STEP 8 for reduced efficiency study
	        FillReducedEfficiency(fillId, AliUEHist::kCFStepBiasStudy, ltRECO, kFALSE);
      	        if (!fReduceMemoryFootprint)
                  FillReducedEfficiency(fillId, AliUEHist::kCFStepBiasStudy2, ltRECO, kTRUE);
        
	        // count number of reco tracks above 150 MeV/c
	        Int_t nRecoTracks = 0; 
	        if (((AliVParticle*) ltRECO->At(0))->Pt() > 0.15)
	          nRecoTracks++;
	        for (Int_t i=0; i<4; i++)
	          for (Int_t j=0; j<((TList*)regionSortedParticlesRECO->At(i))->GetEntries(); j++)
	            if (((AliVParticle*) ((TList*)regionSortedParticlesRECO->At(i))->At(j))->Pt() > 0.15)
	              nRecoTracks++;
	        
	        if (leadingMC && leadingMC->Pt() > 0.5)
	          fHistosUE->GetCorrelationMultiplicity()->Fill(nMCTracks, nRecoTracks);
        
	        if (leadingMC)
	        {
	          // Match reco leading track with true *********************************
	          Int_t recoLabel = ((AliAODTrack*)ltRECO->At(0))->GetLabel();
	          Int_t mcLabel   = ((AliAODMCParticle*)leadingMC)->GetLabel();
	          if (recoLabel != mcLabel) 
	            return;
          
	          // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
	          // (RECO-MATCHED leading particle and RECO all particles)
	          // STEP 7
	          fHistosUE->Fill(fillId,0,AliUEHist::kCFStepRealLeading,(AliVParticle*)ltRECO->At(0),(TList*)regionSortedParticlesRECO->At(0),(TList*)regionSortedParticlesRECO->At(1),(TList*)regionsMinMaxRECO->At(0),(TList*)regionsMinMaxRECO->At(1));
	          // comparing this step with step 6 (for leading-track observables) you get the efficiency to reconstruct the leading track
          	  // comparing this step with step 6 (for all-tracks observables) you see how leading-track misidentification affects the final distributions
        	}
        
	        delete regionSortedParticlesRECOLTMC;
	        delete regionsMinMaxRECOLTMC;
	        delete regionSortedParticlesRECOLTMC2;
	        delete regionsMinMaxRECOLTMC2;
	        delete regionSortedParticlesRECO;
	        delete regionsMinMaxRECO;
	        delete ltRECO;
	      } // lt reco
    	} // vertex
    } // pileup
  } //phyiscs selection
  
  if (ltMC)
    delete ltMC;
  delete regionSortedParticlesMC;
  delete regionsMinMaxMC;
}

//____________________________________________________________________
void  AliAnalysisTaskLeadingTrackUE::AnalyseDataMode()
{

  // Run the analysis on DATA or MC to get raw distributions
 
  PostData(0,fListOfHistos);
  
  if ( fDebug > 3 ) AliInfo( " Processing event in Data mode ..." );
  Int_t eventId = 0;

  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistosUE->FillEvent(eventId, AliUEHist::kCFStepAll);

  // Trigger selection ************************************************
  if (!fAnalyseUE->TriggerSelection(fInputHandler)) return;
  // PILEUP-CUT 
  Bool_t select = kFALSE;
  if (fSelectBit) select = AliAnalysisHelperJetTasks::TestSelectInfo(fSelectBit);
  if (select) 
  {
    fHistosUE->FillEvent(eventId, -2);
    return;
  }
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistosUE->FillEvent(eventId, AliUEHist::kCFStepTriggered);
  
  // Vertex selection *************************************************
  if(!fAnalyseUE->VertexSelection(fAOD, fnTracksVertex, fZVertex)) return;
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistosUE->FillEvent(eventId, AliUEHist::kCFStepVertex);
  // comparing this step with previous one you get the vertex selection efficiency from data (?)
 

  // Get Reconstructed leading particle *******************************
  TObjArray *ltRECO = fAnalyseUE->FindLeadingObjects(fAOD);
  if (!ltRECO) return;
  
  // fill control distributions
  fHistosUE->Fill(0, (AliVParticle*)ltRECO->At(0));
  
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistosUE->FillEvent(eventId, AliUEHist::kCFStepAnaTopology);

  // Switch to reconstructed tracks ************************************
  // Sort RECO particles w.r.t. RECO-leading and return RECO particle  
  TObjArray *regionSortedParticlesRECO = (TObjArray*)fAnalyseUE->SortRegions((AliVParticle*)ltRECO->At(0), fAOD,0); 
  TObjArray *regionsMinMaxRECO = (TObjArray*)fAnalyseUE->GetMinMaxRegion((TList*)regionSortedParticlesRECO->At(2),(TList*)regionSortedParticlesRECO->At(3));
  // Fill UE containers (step, leading track, towards particles, away particles, transverse MIN and MAX particles)
  // (RECO leading particle and RECO  all particles)
  // STEP 6
  fHistosUE->Fill(eventId,0,AliUEHist::kCFStepReconstructed,(AliVParticle*)ltRECO->At(0),(TList*)regionSortedParticlesRECO->At(0),(TList*)regionSortedParticlesRECO->At(1),(TList*)regionsMinMaxRECO->At(0),(TList*)regionsMinMaxRECO->At(1));
  
  // STEP 8 and 9 for reduced efficiency study
  FillReducedEfficiency(eventId, AliUEHist::kCFStepBiasStudy, ltRECO, kFALSE);
  FillReducedEfficiency(eventId, AliUEHist::kCFStepBiasStudy2, ltRECO, kTRUE);
  
  delete ltRECO;
  delete regionSortedParticlesRECO;
  delete regionsMinMaxRECO;
}

//____________________________________________________________________
void AliAnalysisTaskLeadingTrackUE::FillReducedEfficiency(Int_t eventId, AliUEHist::CFStep step, const TObjArray* ltRECO, Bool_t twoStep)
{
  // remove leading particle using fkTrackingEfficiency and use subleading particle to fill the histograms
  //
  // if twoStep is kTRUE, do a two step procedure where in each step only 50% of the loss due to the tracking efficiency is applied 
  
  if (!fkTrackingEfficiency)
    return;
    
  TObjArray* particleList =  new TObjArray(*ltRECO);
  AliVParticle* leading = (AliVParticle*) particleList->At(0);
  if (!leading)
  {
    delete particleList;
    return;
  }
  
  // remove particles depending on tracking efficiency
  Int_t count = (twoStep) ? 2 : 1;
  
  for (Int_t i=0; i<count; i++)
  {
    Float_t trackingEff = fkTrackingEfficiency->GetBinContent(fkTrackingEfficiency->GetXaxis()->FindBin(leading->Pt()));
    if (twoStep)
      trackingEff = 0.5 * (trackingEff + 1);
      
    if (gRandom->Uniform() > trackingEff)
    {
      //Printf("LOWEFF: Removing leading particle");
      particleList->RemoveAt(0);
      particleList->Compress();
    }
      
    if (particleList->GetEntries() == 0)
    {
      delete particleList;
      return;
    }
    
    leading = (AliVParticle*) particleList->At(0);
  }
  
  TObjArray *regionSortedParticlesRECOLowEff = fAnalyseUE->SortRegions(leading, particleList, 0); 
  TObjArray *regionsMinMaxRECOLowEff = fAnalyseUE->GetMinMaxRegion((TList*)regionSortedParticlesRECOLowEff->At(2), (TList*)regionSortedParticlesRECOLowEff->At(3));
    
  fHistosUE->Fill(eventId,0,step,leading,(TList*)regionSortedParticlesRECOLowEff->At(0), (TList*)regionSortedParticlesRECOLowEff->At(1), (TList*)regionsMinMaxRECOLowEff->At(0), (TList*)regionsMinMaxRECOLowEff->At(1));

  delete regionSortedParticlesRECOLowEff;
  delete regionsMinMaxRECOLowEff;
  delete particleList;
}

//____________________________________________________________________
void  AliAnalysisTaskLeadingTrackUE::Initialize()
{
  // input handler
  fInputHandler = (AliInputEventHandler*)
         ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  // MC handler
  fMcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}






