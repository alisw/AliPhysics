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
#include <TMCProcess.h>

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
#include "AliCFContainer.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliStack.h"

#include "AliEventPoolManager.h"


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
fFillMixed(kTRUE),
fCompareCentralities(kFALSE),
fTwoTrackEfficiencyStudy(kFALSE),
// pointers to UE classes
fAnalyseUE(0x0),
fHistos(0x0),
fHistosMixed(0),
fkTrackingEfficiency(0x0),
// handlers and events
fAOD(0x0),
fESD(0x0),
fArrayMC(0x0),
fInputHandler(0x0),
fMcEvent(0x0),
fMcHandler(0x0),
fPoolMgr(0x0),
// histogram settings
fListOfHistos(0x0), 
// event QA
fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default) 
fZVertex(7.),
fCentralityMethod("V0M"),
// track cuts
fTrackEtaCut(0.8),
fPtMin(0.5),
fFilterBit(0xFF),
fSelectBit(0),
fUseChargeHadrons(kFALSE),
fSelectCharge(0),
fFillpT(kFALSE)
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
  
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) 
  { // input AOD
    fAOD = ((AliAODInputHandler*)handler)->GetEvent();
    if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODInputHandler");
  } 
  else 
  {  //output AOD
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if (handler && handler->InheritsFrom("AliAODHandler") ) 
    {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODHandler");
    } 
    else 
    {  // no AOD
      AliWarning("I can't get any AOD Event Handler");
    }
  }
  
  if (handler && handler->InheritsFrom("AliESDInputHandler") ) 
  { // input ESD
    // pointer received per event in ::Exec
    if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliESDInputHandler");
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
  fAnalyseUE->SetParticleSelectionCriteria(fFilterBit, fUseChargeHadrons, fTrackEtaCut, fPtMin);
  fAnalyseUE->SetDebug(fDebug); 
  fAnalyseUE->DefineESDCuts(fFilterBit);

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
  
  fHistos->SetSelectCharge(fSelectCharge);
  fHistosMixed->SetSelectCharge(fSelectCharge);
  
  // add histograms to list
  fListOfHistos->Add(fHistos);
  fListOfHistos->Add(fHistosMixed);
  
  fListOfHistos->Add(new TH2F("trackletsVsV0Cent", ";L1 clusters;v0 centrality", 100, -0.5, 9999.5, 101, 0, 101));
  fListOfHistos->Add(new TH2F("processIDs", ";#Delta#phi;process id", 100, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(), kPNoProcess + 1, -0.5, kPNoProcess + 0.5));
  
  PostData(0,fListOfHistos);
  
  // Add task configuration to output list 
  AddSettingsTree();

  // event mixing
  //  Int_t trackDepth = 100; // Require e.g. 20 5-track events, or 2 50-track events
  Int_t trackDepth = 50000; 
  Int_t poolsize   = 1000;  // Maximum number of events
  
  Int_t nCentralityBins  = fHistos->GetUEHist(2)->GetEventHist()->GetNBins(1);
  Double_t* centralityBins = (Double_t*) fHistos->GetUEHist(2)->GetEventHist()->GetAxis(1, 0)->GetXbins()->GetArray();
  
  Int_t nZvtxBins  = 7;
  Double_t vertexBins[] = { -7, -5, -3, -1, 1, 3, 5, 7 };
  Double_t* zvtxbin = vertexBins;
  
  if (fHistos->GetUEHist(2)->GetEventHist()->GetNVar() > 2)
  {
    nZvtxBins = fHistos->GetUEHist(2)->GetEventHist()->GetNBins(2);
    zvtxbin = (Double_t*) fHistos->GetUEHist(2)->GetEventHist()->GetAxis(2, 0)->GetXbins()->GetArray();
  }

  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nZvtxBins, zvtxbin);
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::Exec(Option_t */*option*/)
{
  // array of MC particles
  if (fMcHandler) {
    if (fAOD)
    {
      fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
	AliFatal("No array of MC particles found !!!");
    }
    fMcEvent = fMcHandler->MCEvent();
  }

  // receive ESD pointer if we are not running AOD analysis
  if (!fAOD)
  {
    AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    if (handler && handler->InheritsFrom("AliESDInputHandler"))
      fESD = ((AliESDInputHandler*)handler)->GetEvent();
  }

  // Analyse the event
  if (fMode) AnalyseCorrectionMode();
  else AnalyseDataMode();
}

/******************** ANALYSIS METHODS *****************************/

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AddSettingsTree()
{
  //Write settings to output list
  TTree *settingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  settingsTree->Branch("fnTracksVertex", &fnTracksVertex,"nTracksVertex/I");
  settingsTree->Branch("fZVertex", &fZVertex,"ZVertex/D");
  //settingsTree->Branch("fCentralityMethod", fCentralityMethod.Data(),"CentralityMethod/C");
  settingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  settingsTree->Branch("fPtMin", &fPtMin, "PtMin/D");
  settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fSelectBit", &fSelectBit,"EventSelectionBit/I");
  settingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  settingsTree->Branch("fSelectCharge", &fSelectCharge,"SelectCharge/I");
  settingsTree->Branch("fFillpT", &fFillpT,"FillpT/O");
  settingsTree->Branch("fkTrackingEfficiency", "TH1D", &fkTrackingEfficiency);
  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
}  

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AnalyseCorrectionMode()
{
  // Run the analysis on MC to get the correction maps
  //
 
  if ( fDebug > 3 ) AliInfo( " Processing event in Corrections mode ..." );
  
  Double_t centrality = 0;
  
  if (fCentralityMethod.Length() > 0)
  {
    AliCentrality *centralityObj = 0;
    if (fAOD)
      centralityObj = fAOD->GetHeader()->GetCentralityP();
    else if (fESD)
      centralityObj = fESD->GetCentrality();
    
    if (centralityObj)
    {
      centrality = centralityObj->GetCentralityPercentileUnchecked(fCentralityMethod);
      AliInfo(Form("Centrality is %f", centrality));
    }
    else
    {
      Printf("WARNING: Centrality object is 0");
      centrality = -1;
     }
  }
  
  // Support for ESD and AOD based analysis
  AliVEvent* inputEvent = fAOD;
  if (!inputEvent)
    inputEvent = fESD;
  
  TObject* mc = fArrayMC;
  if (!mc)
    mc = fMcEvent;
  
  // count all events
  fHistos->FillEvent(centrality, -1);
  
  // Only consider MC events within the vtx-z region used also as cut on the reconstructed vertex
  if (!fAnalyseUE->VertexSelection(fMcEvent, 0, fZVertex)) 
    return;
  
  Float_t zVtx = fMcEvent->GetPrimaryVertex()->GetZ();
  Float_t weight = 1;
  if (fFillpT)
    weight = -1;
    
  // Get MC primaries
  TObjArray* tracksMC = fAnalyseUE->GetAcceptedParticles(mc, 0, kTRUE, -1, kTRUE);
  
  // (MC-true all particles)
  // STEP 0
  fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepAll, tracksMC, 0, weight);
  
  // Trigger selection ************************************************
  if (fAnalyseUE->TriggerSelection(fInputHandler))
  {  
    // (MC-true all particles)
    // STEP 1
    if (!fReduceMemoryFootprint)
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTriggered, tracksMC, 0, weight);
    else
      fHistos->FillEvent(centrality, AliUEHist::kCFStepTriggered);
      
    if (!inputEvent) {
      AliFatal("UNEXPECTED: inputEvent is 0. Trigger selection should have failed");
      return;
    }
    
    // Vertex selection *************************************************
    if (fAnalyseUE->VertexSelection(inputEvent, fnTracksVertex, fZVertex))
    {
      // fill here for tracking efficiency
      // loop over particle species
      
      for (Int_t particleSpecies = 0; particleSpecies < 4; particleSpecies++)
      {
        TObjArray* primMCParticles = fAnalyseUE->GetAcceptedParticles(mc, 0x0, kTRUE, particleSpecies, kTRUE);
        TObjArray* primRecoTracksMatched = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kTRUE, particleSpecies, kTRUE);
        TObjArray* allRecoTracksMatched = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kFALSE, particleSpecies, kTRUE);
      
        fHistos->FillTrackingEfficiency(primMCParticles, primRecoTracksMatched, allRecoTracksMatched, particleSpecies, centrality);
        
// 	Printf("%d --> %d %d %d", particleSpecies, primMCParticles->GetEntries(), primRecoTracksMatched->GetEntries(), allRecoTracksMatched->GetEntries());

	delete primMCParticles;
        delete primRecoTracksMatched;
        delete allRecoTracksMatched;
      }
    
      // (MC-true all particles)
      // STEP 2
      if (!fReduceMemoryFootprint)
	fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepVertex, tracksMC, 0, weight);
      else
	fHistos->FillEvent(centrality, AliUEHist::kCFStepVertex);
      
      // Get MC primaries that match reconstructed track
      TObjArray* tracksRecoMatchedPrim = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kTRUE, -1, kTRUE);
      
      // (RECO-matched (quantities from MC particle) primary particles)
      // STEP 4
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTrackedOnlyPrim, tracksRecoMatchedPrim, 0, weight);
      
      // Get MC primaries + secondaries that match reconstructed track
      TObjArray* tracksRecoMatchedAll = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kFALSE, -1, kTRUE);
      
      // (RECO-matched (quantities from MC particle) all particles)
      // STEP 5
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTracked, tracksRecoMatchedAll, 0, weight);
      
      // Get RECO tracks
      TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, -1, kTRUE);
      
      // (RECO all tracks)
      // STEP 6
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, 0, weight);
      
      if (0 && !fReduceMemoryFootprint)
      {
        // make list of secondaries (matched with MC)
        TObjArray* tracksRecoMatchedSecondaries = new TObjArray;
        for (Int_t i=0; i<tracksRecoMatchedAll->GetEntries(); i++)
          if (((AliAODMCParticle*)tracksRecoMatchedAll->At(i))->IsPhysicalPrimary() == kFALSE)
            tracksRecoMatchedSecondaries->Add(tracksRecoMatchedAll->At(i));
      
        // Study: Use only secondaries as trigger particles and plot the correlation vs. all particles; store in step 9
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy2, tracksRecoMatchedSecondaries, tracksRecoMatchedAll, weight);
        
        // Study: Use only primaries as trigger particles and plot the correlation vs. secondaries; store in step 8
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracksRecoMatchedPrim, tracksRecoMatchedSecondaries, weight);
      
        // plot delta phi vs process id of secondaries
        // trigger particles: primaries in 4 < pT < 10
        // associated particles: secondaries in 1 < pT < 10
        
        for (Int_t i=0; i<tracksRecoMatchedPrim->GetEntries(); i++)
        {
          AliVParticle* triggerParticle = (AliVParticle*) tracksRecoMatchedPrim->At(i);
          
          if (triggerParticle->Pt() < 4 || triggerParticle->Pt() > 10)
            continue;
          
          for (Int_t j=0; j<tracksRecoMatchedSecondaries->GetEntries(); j++)
          {
            AliAODMCParticle* particle = (AliAODMCParticle*) tracksRecoMatchedSecondaries->At(j);
            
            if (particle->Pt() < 1 || particle->Pt() > 10)
              continue;
            
            if (particle->Pt() > triggerParticle->Pt())
              continue;
              
            Double_t deltaPhi = triggerParticle->Phi() - particle->Phi();
            if (deltaPhi > 1.5 * TMath::Pi()) 
              deltaPhi -= TMath::TwoPi();
            if (deltaPhi < -0.5 * TMath::Pi())
              deltaPhi += TMath::TwoPi();
              
            Int_t processID = fMcEvent->Stack()->Particle(particle->GetLabel())->GetUniqueID();
              
            ((TH2F*) fListOfHistos->FindObject("processIDs"))->Fill(deltaPhi, processID);
          }
        }
        
        delete tracksRecoMatchedSecondaries;
      }
  
      delete tracksRecoMatchedPrim;
      delete tracksRecoMatchedAll;
      delete tracks;
    }
  }
  
  delete tracksMC;
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AnalyseDataMode()
{

  // Run the analysis on DATA or MC to get raw distributions
 
  if ( fDebug > 3 ) AliInfo( " Processing event in Data mode ..." );

  // skip not selected events here (the AOD is not updated for those)
  if (!fInputHandler)
    return;
    
  if (!(fInputHandler->IsEventSelected() & (AliVEvent::kMB | AliVEvent::kUserDefined)))
    return;

  Double_t centrality = 0;
  
  AliCentrality *centralityObj = 0;
  if (fCentralityMethod.Length() > 0)
  {
    if (fAOD)
      centralityObj = fAOD->GetHeader()->GetCentralityP();
    else if (fESD)
      centralityObj = fESD->GetCentrality();
    
    if (centralityObj)
      centrality = centralityObj->GetCentralityPercentile(fCentralityMethod);
      //centrality = centralityObj->GetCentralityPercentileUnchecked(fCentralityMethod);
    else
      centrality = -1;
    AliInfo(Form("Centrality is %f", centrality));
  }
  
  // Support for ESD and AOD based analysis
  AliVEvent* inputEvent = fAOD;
  if (!inputEvent)
    inputEvent = fESD;

  fHistos->SetRunNumber(inputEvent->GetRunNumber());
  
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(centrality, AliUEHist::kCFStepAll);
  
  // Trigger selection ************************************************
  if (!fAnalyseUE->TriggerSelection(fInputHandler)) return;
  
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(centrality, AliUEHist::kCFStepTriggered);
  
  // Vertex selection *************************************************
  if(!fAnalyseUE->VertexSelection(inputEvent, fnTracksVertex, fZVertex)) return;
  
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(centrality, AliUEHist::kCFStepVertex);
 
  // optimization
  if (centrality < 0 && !fCompareCentralities)
    return;

  TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, -1, kTRUE);
  //Printf("Accepted %d tracks", tracks->GetEntries());
  
  const AliVVertex* vertex = inputEvent->GetPrimaryVertex();
  Double_t zVtx = vertex->GetZ();
  
  Float_t weight = 1;
  if (fFillpT)
    weight = -1;

  // fill two different centralities (syst study)
  // the zvtx axis is used to distinguish the results of both centralities: configured centrality in zvtx = 0, SPD in zvtx = 2
  if (fCompareCentralities)
    zVtx = 0;
  
  // Fill containers at STEP 6 (reconstructed)
  if (centrality >= 0)
    fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, 0, weight);

  // Two-track effect study
  if (fTwoTrackEfficiencyStudy)
  {
    Float_t bSign = (fESD->GetMagneticField() > 0) ? 1 : -1;
    fHistos->TwoTrackEfficiency(tracks, 0, bSign);
  }
  
  // fill second time with SPD centrality
  if (fCompareCentralities && centralityObj)
  {
    centrality = centralityObj->GetCentralityPercentile("CL1");
    if (centrality >= 0)
      fHistos->FillCorrelations(centrality, 2, AliUEHist::kCFStepReconstructed, tracks, 0, weight);
  }
    
  if (fFillMixed)
  {
    // event mixing
    
    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.

    AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx);
    
    if (!pool)
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, zVtx));
    
    //pool->SetDebug(1);
      
    if (pool->IsReady()) 
    {
      
      Int_t nMix = pool->GetCurrentNEvents();
      //cout << "nMix = " << nMix << endl;
    
      // Fill mixed-event histos here  
      for (Int_t jMix=0; jMix<nMix; jMix++) 
      {
	TObjArray* bgTracks = pool->GetEvent(jMix);
	
	if (fTwoTrackEfficiencyStudy)
	{
	  Float_t bSign = (fESD->GetMagneticField() > 0) ? 1 : -1;
	  fHistos->TwoTrackEfficiency(tracks, bgTracks, bSign);
	}

	fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, bgTracks, 1.0 / nMix, (jMix == 0));
      }
    }
    
    // clone and give ownership to event pool
    TObjArray* tracksClone = (TObjArray*) tracks->Clone();
    tracksClone->SetOwner(kTRUE);
    
    pool->UpdatePool(tracksClone);
    //pool->PrintInfo();
  }

  delete tracks;
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
