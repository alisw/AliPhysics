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

#include "AliAnalysisTaskPhiCorrelations.h"
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

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"

#include "EventMixing/AliMixEventInputHandler.h"


////////////////////////////////////////////////////////////////////////
//
// Analysis class for azimuthal correlation studies
// Based on the UE task from Sara Vallero and Jan Fiete
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
//    Jan Fiete Grosse-Oetringhaus
// 
////////////////////////////////////////////////////////////////////////


ClassImp( AliAnalysisTaskPhiCorrelations )

//____________________________________________________________________
AliAnalysisTaskPhiCorrelations:: AliAnalysisTaskPhiCorrelations(const char* name):
AliAnalysisTask(name,""),
// general configuration
fDebug(0),
fMode(0),
fReduceMemoryFootprint(kFALSE),
// pointers to UE classes
fAnalyseUE(0x0),
fHistos(0x0),
fHistosMixed(0),
fkTrackingEfficiency(0x0),
// handlers and events
fAOD(0x0),           
fArrayMC(0x0),
fInputHandler(0x0),
fMcEvent(0x0),
fMcHandler(0x0),
// histogram settings
fListOfHistos(0x0), 
// event QA
fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default) 
fZVertex(10.),
// track cuts
fTrackEtaCut(0.8),
fPtMin(0.5),
fFilterBit(0xFF),
fSelectBit(0),
fUseChargeHadrons(kFALSE)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());

}

AliAnalysisTaskPhiCorrelations::~AliAnalysisTaskPhiCorrelations() 
{ 
  // destructor
  
  if (fListOfHistos  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) 
    delete fListOfHistos;
}

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::ConnectInputData(Option_t* /*option*/)
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
void  AliAnalysisTaskPhiCorrelations::CreateOutputObjects()
{
  // Create the output container
  
  if (fDebug > 1) AliInfo("CreateOutputObjects()");
   
  // Initialize class with main algorithms, event and track selection. 
  fAnalyseUE = new AliAnalyseLeadingTrackUE();
  fAnalyseUE->SetParticleSelectionCriteria(fFilterBit, fUseChargeHadrons, fTrackEtaCut, 1.0);
  fAnalyseUE->SetDebug(fDebug); 
  fAnalyseUE->DefineESDCuts(0);

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
  fHistos = new AliUEHistograms("AliUEHistogramsSame", "4");
  fHistosMixed = new AliUEHistograms("AliUEHistogramsMixed", "4");
  
  // add histograms to list
  fListOfHistos->Add(fHistos);
  fListOfHistos->Add(fHistosMixed);
  
  //fListOfHistos->Add(new TH2F("multVsLeadStep5", ";multiplicity;leading pT", 100, -0.5, 99.5, 20, 0, 10));
  
  // Add task configuration to output list 
  AddSettingsTree();

  PostData(0,fListOfHistos);
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::Exec(Option_t */*option*/)
{
  // array of MC particles
  if (fMcHandler){
  	fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  	if (!fArrayMC)AliFatal("No array of MC particles found !!!");
	}

  // Analyse the event
  if (fMode) AnalyseCorrectionMode();
  else AnalyseDataMode();

  PostData(0,fListOfHistos);
}

/******************** ANALYSIS METHODS *****************************/

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AddSettingsTree()
{
  //Write settings to output list
  TTree *settingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fSelectBit", &fSelectBit,"EventSelectionBit/I");
  settingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  settingsTree->Branch("fPtMin", &fPtMin, "PtMin/D");
  settingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  settingsTree->Branch("fkTrackingEfficiency", "TH1D", &fkTrackingEfficiency);
  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
}  

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AnalyseCorrectionMode()
{
  // Run the analysis on MC to get the correction maps
  //
 
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AnalyseDataMode()
{

  // Run the analysis on DATA or MC to get raw distributions
 
  PostData(0,fListOfHistos);
  
  if ( fDebug > 3 ) AliInfo( " Processing event in Data mode ..." );
  Int_t eventId = 0;

  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(eventId, AliUEHist::kCFStepAll);

  // Trigger selection ************************************************
  if (!fAnalyseUE->TriggerSelection(fInputHandler)) return;
  
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(eventId, AliUEHist::kCFStepTriggered);
  
  // Vertex selection *************************************************
  if(!fAnalyseUE->VertexSelection(fAOD, fnTracksVertex, fZVertex)) return;
  
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(eventId, AliUEHist::kCFStepVertex);
 
  // TODO HACK centrality should be retrieved from AOD, this needs the ESD
  AliESDInputHandler* handler = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliESDEvent* esd = handler->GetEvent();
  Int_t centrality = esd->GetMultiplicity()->GetNumberOfITSClusters(1);

  TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(fAOD, 0, kTRUE, -1, kTRUE);
  
  // reduce to pt and eta range
  /*
  TObjArray* reduced = new TObjArray;
  for (Int_t i=0; i<tracks->GetEntries(); i++)
  {
    AliVParticle* particle = (AliVParticle*) tracks->At(i);
    if (TMath::Abs(particle->Eta()) < fTrackEtaCut && particle->Pt() > fPtMin)
      reduced->Add(particle);
  }
  */
  
  // Fill containers at STEP 6 (reconstructed)
  fHistos->FillCorrelations(eventId, centrality, AliUEHist::kCFStepReconstructed, tracks);
    
  AliInfo("**Main Event");
  AliInfo(Form("  Tracklets %d", fAOD->GetTracklets()->GetNumberOfTracklets()));
  AliInfo(Form("  Vz %f", fAOD->GetPrimaryVertex()->GetZ()));  

  AliMixEventInputHandler* mixEH = (AliMixEventInputHandler*) fInputHandler->MixingHandler();
  if (mixEH)
  {
    mixEH->GetEntry();
    AliInfo(Form("num mixed %d", mixEH->MixedEventNumber()));
          
    if (mixEH->MixedEventNumber() > 0) 
    {
      // TODO in principle need to apply same quality cuts as above
      // or somewhere before on the level of the pool (would be more efficient)
      for (Int_t i = 0; i < mixEH->BufferSize(); i++) 
      {
        AliESDEvent *eventMix = (AliESDEvent*) mixEH->InputEventHandler(i)->GetEvent();
        if (!eventMix) {
          AliError(Form("Could not retrieve event %d", i));
        }
        else 
        {
          AliInfo(Form("**Mixed Event %d", i));
          AliInfo(Form("  Tracklets %d", eventMix->GetMultiplicity()->GetNumberOfTracklets()));
          AliInfo(Form("  Vz %f", eventMix->GetPrimaryVertex()->GetZ()));
          
          TObjArray* tracksMixed = fAnalyseUE->GetAcceptedParticles(eventMix, 0, kTRUE, -1, kTRUE);
          
          fHistosMixed->FillCorrelations(eventId, centrality, AliUEHist::kCFStepReconstructed, tracks, tracksMixed);
          
          tracksMixed->SetOwner(); // contains tpc only tracks, that is why we have to delete the content as well
          delete tracksMixed;
          
        }
      }
    }
  }
  
  delete tracks;
  //delete reduced;
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::Initialize()
{
  // input handler
  fInputHandler = (AliInputEventHandler*)
         ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  // MC handler
  fMcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}






