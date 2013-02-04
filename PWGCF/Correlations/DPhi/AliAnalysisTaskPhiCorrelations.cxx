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
#include <TH3F.h>
#include <TRandom.h>

#include "AliAnalysisTaskPhiCorrelations.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliUEHistograms.h"
#include "AliUEHist.h"

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
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
#include "AliAODMCHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"

#include "AliEventPoolManager.h"

#include "AliESDZDC.h"
#include "AliESDtrackCuts.h"


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
ClassImp( AliDPhiBasicParticle )

//____________________________________________________________________
AliAnalysisTaskPhiCorrelations:: AliAnalysisTaskPhiCorrelations(const char* name):
AliAnalysisTask(name,""),
// general configuration
fDebug(0),
fMode(0),
fReduceMemoryFootprint(kFALSE),
fFillMixed(kTRUE),
fMixingTracks(50000),
fCompareCentralities(kFALSE),
fTwoTrackEfficiencyStudy(kFALSE),
fTwoTrackEfficiencyCut(0),
fUseVtxAxis(kFALSE),
fCourseCentralityBinning(kFALSE),
fSkipTrigger(kFALSE),
fInjectedSignals(kFALSE),
// pointers to UE classes
fAnalyseUE(0x0),
fHistos(0x0),
fHistosMixed(0),
fEfficiencyCorrection(0),
fCorrectTriggers(kFALSE),
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
fOnlyOneEtaSide(0),
fPtMin(0.5),
fFilterBit(0xFF),
fSelectBit(AliVEvent::kMB|AliVEvent::kUserDefined),
fUseChargeHadrons(kFALSE),
fSelectParticleSpecies(-1),
fSelectCharge(0),
fTriggerSelectCharge(0),
fTriggerRestrictEta(-1),
fEtaOrdering(kFALSE),
fCutConversions(kFALSE),
fCutResonances(kFALSE),
fFillOnlyStep0(kFALSE),
fSkipStep6(kFALSE),
fRejectCentralityOutliers(kFALSE),
fRemoveWeakDecays(kFALSE),
fRemoveDuplicates(kFALSE),
fSkipFastCluster(kFALSE),
fWeightPerEvent(kFALSE),
fCustomBinning(),
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
  fAnalyseUE->SetEventSelection(fSelectBit);

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
  TString histType = "4R";
  if (fUseVtxAxis == 1)
    histType = "5R";
  else if (fUseVtxAxis == 2)
    histType = "6R";
  if (fCourseCentralityBinning)
    histType += "C";
  fHistos = new AliUEHistograms("AliUEHistogramsSame", histType, fCustomBinning);
  fHistosMixed = new AliUEHistograms("AliUEHistogramsMixed", histType, fCustomBinning);
  
  fHistos->SetSelectCharge(fSelectCharge);
  fHistosMixed->SetSelectCharge(fSelectCharge);
  
  fHistos->SetSelectTriggerCharge(fTriggerSelectCharge);
  fHistosMixed->SetSelectTriggerCharge(fTriggerSelectCharge);

  fHistos->SetTriggerRestrictEta(fTriggerRestrictEta);
  fHistosMixed->SetTriggerRestrictEta(fTriggerRestrictEta);
  
  fHistos->SetOnlyOneEtaSide(fOnlyOneEtaSide);
  fHistosMixed->SetOnlyOneEtaSide(fOnlyOneEtaSide);
  
  fHistos->SetEtaOrdering(fEtaOrdering);
  fHistosMixed->SetEtaOrdering(fEtaOrdering);

  fHistos->SetPairCuts(fCutConversions, fCutResonances);
  fHistosMixed->SetPairCuts(fCutConversions, fCutResonances);
  
  fHistos->SetTrackEtaCut(fTrackEtaCut);
  fHistosMixed->SetTrackEtaCut(fTrackEtaCut);
  
  fHistos->SetWeightPerEvent(fWeightPerEvent);
  fHistosMixed->SetWeightPerEvent(fWeightPerEvent);

  if (fEfficiencyCorrection)
  {
    fHistos->SetEfficiencyCorrection(fEfficiencyCorrection, fCorrectTriggers);
    fHistosMixed->SetEfficiencyCorrection((THnF*) fEfficiencyCorrection->Clone(), fCorrectTriggers);
  }
  
  // add histograms to list
  fListOfHistos->Add(fHistos);
  fListOfHistos->Add(fHistosMixed);
  
  fListOfHistos->Add(new TH2F("trackletsVsV0Cent", ";L1 clusters;v0 centrality", 100, -0.5, 9999.5, 101, 0, 101));
  fListOfHistos->Add(new TH2F("processIDs", ";#Delta#phi;process id", 100, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(), kPNoProcess + 1, -0.5, kPNoProcess + 0.5));
  fListOfHistos->Add(new TH1F("eventStat", ";;events", 4, -0.5, 3.5));
  fListOfHistos->Add(new TH2F("mixedDist", ";centrality;tracks;events", 101, 0, 101, 200, 0, fMixingTracks * 1.5));
  fListOfHistos->Add(new TH1F("pids", ";pdg;tracks", 2001, -1000.5, 1000.5));
  fListOfHistos->Add(new TH2F("referenceMultiplicity", ";centrality;tracks;events", 101, 0, 101, 200, 0, 200));
  
  PostData(0,fListOfHistos);
  
  // Add task configuration to output list 
  AddSettingsTree();

  // event mixing
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemention of AliEventPoolManager
   
  Int_t nCentralityBins  = fHistos->GetUEHist(2)->GetEventHist()->GetNBins(1);
  Double_t* centralityBins = (Double_t*) fHistos->GetUEHist(2)->GetEventHist()->GetAxis(1, 0)->GetXbins()->GetArray();
  
  const Int_t kNZvtxBins  = 10+(1+10)*4;
  // bins for further buffers are shifted by 100 cm
  Double_t vertexBins[kNZvtxBins+1] = { -10,   -8,  -6,  -4,  -2,   0,   2,   4,   6,   8,  10, 
				       90,  92,  94,  96,  98, 100, 102, 104, 106, 108, 110, 
				      190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 
				      290, 292, 294, 296, 298, 300, 302, 304, 306, 308, 310, 
				      390, 392, 394, 396, 398, 400, 402, 404, 406, 408, 410 };

  Int_t nZvtxBins  = kNZvtxBins;
  Double_t* zvtxbin = vertexBins;

  if (fMode == 0 && fHistos->GetUEHist(2)->GetEventHist()->GetNVar() > 2)
  {
    nZvtxBins = fHistos->GetUEHist(2)->GetEventHist()->GetNBins(2);
    zvtxbin = (Double_t*) fHistos->GetUEHist(2)->GetEventHist()->GetAxis(2, 0)->GetXbins()->GetArray();
  }

  fPoolMgr = new AliEventPoolManager(poolsize, fMixingTracks, nCentralityBins, centralityBins, nZvtxBins, zvtxbin);
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::Exec(Option_t */*option*/)
{
  // receive ESD pointer if we are not running AOD analysis
  if (!fAOD)
  {
    AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    if (handler && handler->InheritsFrom("AliESDInputHandler"))
      fESD = ((AliESDInputHandler*)handler)->GetEvent();
  }

  if (fMode)
  {
    // correction mode
    
    if (fMcHandler)
      fMcEvent = fMcHandler->MCEvent();

    if (fAOD)
    {
      // array of MC particles
      fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
	AliFatal("No array of MC particles found !!!");
    }
    
    AnalyseCorrectionMode();
  }
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
  settingsTree->Branch("fOnlyOneEtaSide", &fOnlyOneEtaSide,"OnlyOneEtaSide/I");
  settingsTree->Branch("fPtMin", &fPtMin, "PtMin/D");
  settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fSelectBit", &fSelectBit,"EventSelectionBit/I");
  settingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  settingsTree->Branch("fSelectParticleSpecies", &fSelectParticleSpecies,"ParticleSpecies/I");
  settingsTree->Branch("fSelectCharge", &fSelectCharge,"SelectCharge/I");
  settingsTree->Branch("fTriggerSelectCharge", &fTriggerSelectCharge,"TriggerSelectCharge/I");
  settingsTree->Branch("fTriggerRestrictEta", &fTriggerRestrictEta,"TriggerRestrictEta/D");
  settingsTree->Branch("fEtaOrdering", &fEtaOrdering,"EtaOrdering/O");
  settingsTree->Branch("fCutConversions", &fCutConversions,"CutConversions/O");
  settingsTree->Branch("fCutResonances", &fCutResonances,"CutResonances/O");
  settingsTree->Branch("fFillpT", &fFillpT,"FillpT/O");
  settingsTree->Branch("fMixingTracks", &fMixingTracks,"MixingTracks/I");
  settingsTree->Branch("fSkipTrigger", &fSkipTrigger,"SkipTrigger/O");
  settingsTree->Branch("fInjectedSignals", &fInjectedSignals,"SkipTrigger/O");
  settingsTree->Branch("fRejectCentralityOutliers", &fRejectCentralityOutliers,"RejectCentralityOutliers/O");
  settingsTree->Branch("fRemoveWeakDecays", &fRemoveWeakDecays,"RemoveWeakDecays/O");
  settingsTree->Branch("fRemoveDuplicates", &fRemoveDuplicates,"RemoveDuplicates/O");
  settingsTree->Branch("fSkipFastCluster", &fSkipFastCluster,"SkipFastCluster/O");
  settingsTree->Branch("fWeightPerEvent", &fWeightPerEvent,"WeightPerEvent/O");
  //fCustomBinning
  settingsTree->Branch("fCorrectTriggers", &fCorrectTriggers,"CorrectTriggers/O");
  
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

  Float_t bSign = 0;
  
  if (inputEvent)
  {
    fHistos->SetRunNumber(inputEvent->GetRunNumber());
    bSign = (inputEvent->GetMagneticField() > 0) ? 1 : -1;
  }
    
  TObject* mc = fArrayMC;
  if (!mc)
    mc = fMcEvent;
  
  // count all events
  fHistos->FillEvent(centrality, -1);
  
  if (centrality < 0)
    return;

  // Only consider MC events within the vtx-z region used also as cut on the reconstructed vertex
  TObject* vertexSupplier = fMcEvent;
  if (fAOD) // AOD
    vertexSupplier = fAOD->FindListObject(AliAODMCHeader::StdBranchName());
    
  if (!fAnalyseUE->VertexSelection(vertexSupplier, 0, fZVertex)) 
    return;
  
  Float_t zVtx = 0;
  if (fAOD)
    zVtx = ((AliAODMCHeader*) vertexSupplier)->GetVtxZ();
  else
    zVtx = fMcEvent->GetPrimaryVertex()->GetZ();
  Float_t weight = 1;
  if (fFillpT)
    weight = -1;
  
  // For productions with injected signals, figure out above which label to skip particles/tracks
  Int_t skipParticlesAbove = 0;
  if (fInjectedSignals)
  {
    AliGenEventHeader* eventHeader = 0;
    Int_t headers = 0;

    if (fMcEvent)
    {
      // ESD
      AliHeader* header = (AliHeader*) fMcEvent->Header();
      if (!header)
	AliFatal("fInjectedSignals set but no MC header found");
	
      AliGenCocktailEventHeader* cocktailHeader = dynamic_cast<AliGenCocktailEventHeader*> (header->GenEventHeader());
      if (!cocktailHeader)
      {
	header->Dump();
	AliFatal("fInjectedSignals set but no MC cocktail header found");
      }

      headers = cocktailHeader->GetHeaders()->GetEntries();
      eventHeader = dynamic_cast<AliGenEventHeader*> (cocktailHeader->GetHeaders()->First());
    }
    else
    {
      // AOD
      AliAODMCHeader* header = (AliAODMCHeader*) fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      if (!header)
	AliFatal("fInjectedSignals set but no MC header found");
      
      headers = header->GetNCocktailHeaders();
      eventHeader = header->GetCocktailHeader(0);
    }
    
    if (!eventHeader)
    {
      // We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
      // (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
      AliError("First event header not found. Skipping this event.");
      fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
      return;
    }
    
    skipParticlesAbove = eventHeader->NProduced();
    AliInfo(Form("Injected signals in this event (%d headers). Keeping events of %s. Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove));
  }
  
  // Get MC primaries
  // triggers
  TObjArray* tmpList = fAnalyseUE->GetAcceptedParticles(mc, 0, kTRUE, -1, kTRUE);
  CleanUp(tmpList, mc, skipParticlesAbove);
  TObjArray* tracksMC = CloneAndReduceTrackList(tmpList);
  delete tmpList;
  
  // associated
  TObjArray* tracksCorrelateMC = tracksMC;
  if (fSelectParticleSpecies != -1)
  {
    // TODO for MC this uses to PDG of the mother of the particle
    tracksCorrelateMC = fAnalyseUE->GetAcceptedParticles(mc, 0, kTRUE, fSelectParticleSpecies, kTRUE);
    CleanUp(tracksCorrelateMC, mc, skipParticlesAbove);
  }
  
  /*
  if (fAOD)
  {
    for (Int_t i=0; i<fArrayMC->GetEntriesFast(); i++)
      ((TH1F*) fListOfHistos->FindObject("pids"))->Fill(((AliAODMCParticle*) fArrayMC->At(i))->PdgCode());
  }
  else
  {
    for (Int_t i=0; i<fMcEvent->GetNumberOfTracks(); i++)
      ((TH1F*) fListOfHistos->FindObject("pids"))->Fill(fMcEvent->GetTrack(i)->PdgCode());
  }
  */
  
  if (fFillOnlyStep0)
    zVtx = 0;
  
  // (MC-true all particles)
  // STEP 0
  fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepAll, tracksMC, tracksCorrelateMC, weight);
  
  // mixed event
  if (fFillMixed)
  {
    AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx);
    //pool->PrintInfo();
    if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5) 
      for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) 
	fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepAll, tracksMC, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix == 0));
    pool->UpdatePool(CloneAndReduceTrackList(tracksCorrelateMC));
  }
  
//   Printf("trigger: %d", ((AliInputEventHandler*)fInputHandler)->IsEventSelected());
  
  // Trigger selection ************************************************
  if (!fFillOnlyStep0 && (fSkipTrigger || fAnalyseUE->TriggerSelection(fInputHandler)))
  {  
    // (MC-true all particles)
    // STEP 1
    if (!fReduceMemoryFootprint)
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTriggered, tracksMC, tracksCorrelateMC, weight);
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
	
	CleanUp(primMCParticles, mc, skipParticlesAbove);
	CleanUp(primRecoTracksMatched, mc, skipParticlesAbove);
	CleanUp(allRecoTracksMatched, mc, skipParticlesAbove);
      
        fHistos->FillTrackingEfficiency(primMCParticles, primRecoTracksMatched, allRecoTracksMatched, 0, particleSpecies, centrality, zVtx);
        
// 	Printf("%d --> %d %d %d", particleSpecies, primMCParticles->GetEntries(), primRecoTracksMatched->GetEntries(), allRecoTracksMatched->GetEntries());

	delete primMCParticles;
        delete primRecoTracksMatched;
        delete allRecoTracksMatched;
      }
      
      TObjArray* fakeParticles = fAnalyseUE->GetFakeParticles(inputEvent, mc, kFALSE, -1, kTRUE);
      CleanUp((TObjArray*) fakeParticles->At(0), mc, skipParticlesAbove);
      CleanUp((TObjArray*) fakeParticles->At(1), mc, skipParticlesAbove);

      fHistos->FillTrackingEfficiency(0, 0, 0, (TObjArray*) fakeParticles->At(2), -1, centrality, zVtx);
      fHistos->FillFakePt(fakeParticles, centrality);
//       Printf(">>>>> %d %d %d fakes", ((TObjArray*) fakeParticles->At(0))->GetEntriesFast(), ((TObjArray*) fakeParticles->At(1))->GetEntriesFast(), ((TObjArray*) fakeParticles->At(2))->GetEntriesFast());
      delete fakeParticles;
    
      // (MC-true all particles)
      // STEP 2
      if (!fReduceMemoryFootprint)
	fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepVertex, tracksMC, tracksCorrelateMC, weight);
      else
	fHistos->FillEvent(centrality, AliUEHist::kCFStepVertex);
      
      // Get MC primaries that match reconstructed track
      // triggers
      TObjArray* tracksRecoMatchedPrim = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kTRUE, -1, kTRUE);
      CleanUp(tracksRecoMatchedPrim, mc, skipParticlesAbove);
      // associated
      TObjArray* tracksCorrelateRecoMatchedPrim = tracksRecoMatchedPrim;
      if (fSelectParticleSpecies != -1)
      {
	tracksCorrelateRecoMatchedPrim = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kTRUE, fSelectParticleSpecies, kTRUE);
	CleanUp(tracksCorrelateRecoMatchedPrim, mc, skipParticlesAbove);
      }

      // (RECO-matched (quantities from MC particle) primary particles)
      // STEP 4
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTrackedOnlyPrim, tracksRecoMatchedPrim, tracksCorrelateRecoMatchedPrim, weight);

      // mixed event
      if (fFillMixed)
      {
	AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx + 200);
	//pool->PrintInfo();
	if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5) 
	  for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) 
	    fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTrackedOnlyPrim, tracksRecoMatchedPrim, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix == 0));
	pool->UpdatePool(CloneAndReduceTrackList(tracksCorrelateRecoMatchedPrim));
      }
      
      // Get MC primaries + secondaries that match reconstructed track
      // triggers
      TObjArray* tracksRecoMatchedAll = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kFALSE, -1, kTRUE);
      CleanUp(tracksRecoMatchedAll, mc, skipParticlesAbove);
      // associated
      TObjArray* tracksCorrelateRecoMatchedAll = tracksRecoMatchedAll;
      if (fSelectParticleSpecies != -1)
      {
	tracksCorrelateRecoMatchedAll = fAnalyseUE->GetAcceptedParticles(inputEvent, mc, kFALSE, fSelectParticleSpecies, kTRUE);
	CleanUp(tracksCorrelateRecoMatchedAll, mc, skipParticlesAbove);
      }
      
      // (RECO-matched (quantities from MC particle) all particles)
      // STEP 5
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTracked, tracksRecoMatchedAll, tracksCorrelateRecoMatchedAll, weight);
      
      // mixed event
      if (fFillMixed)
      {
	AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zVtx + 300);
	//pool->PrintInfo();
	if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5) 
	  for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) 
	    fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepTracked, tracksRecoMatchedAll, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix == 0));
	pool->UpdatePool(CloneAndReduceTrackList(tracksCorrelateRecoMatchedAll));
      }
      
      // Get RECO tracks
      // triggers
      TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, -1, kTRUE);
      CleanUp(tracks, mc, skipParticlesAbove);
      // associated
      TObjArray* tracksCorrelate = tracks;
      if (fSelectParticleSpecies != -1)
      {
	tracksCorrelate = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, fSelectParticleSpecies, kTRUE);
	CleanUp(tracksCorrelate, mc, skipParticlesAbove);
      }
     
      // (RECO all tracks)
      // STEP 6
      if (!fSkipStep6)
	fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, tracksCorrelate, weight);
      
      // two track cut, STEP 8
      if (fTwoTrackEfficiencyCut > 0)
	fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracks, tracksCorrelate, weight, kTRUE, kTRUE, bSign, fTwoTrackEfficiencyCut);

      // apply correction efficiency, STEP 10
      if (fEfficiencyCorrection)
      {
	// with or without two track efficiency depending on if fTwoTrackEfficiencyCut is set
	Bool_t twoTrackCut = (fTwoTrackEfficiencyCut > 0);
	
	fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepCorrected, tracks, tracksCorrelate, weight, kTRUE, twoTrackCut, bSign, fTwoTrackEfficiencyCut, kTRUE);
      }

      // mixed event
      if (fFillMixed)
      {
	AliEventPool* pool2 = fPoolMgr->GetEventPool(centrality, zVtx + 100);
	//pool2->PrintInfo();
	if (pool2->IsReady() || pool2->NTracksInPool() > fMixingTracks / 10 || pool2->GetCurrentNEvents() >= 5) 
	{
	  for (Int_t jMix=0; jMix<pool2->GetCurrentNEvents(); jMix++)
	  {
	    // STEP 6
	    if (!fSkipStep6)
	      fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, pool2->GetEvent(jMix), 1.0 / pool2->GetCurrentNEvents(), (jMix == 0));
	    
	    // two track cut, STEP 8
	    if (fTwoTrackEfficiencyCut > 0)
	      fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracks, pool2->GetEvent(jMix), 1.0 / pool2->GetCurrentNEvents(), (jMix == 0), kTRUE, bSign, fTwoTrackEfficiencyCut);
	    
	    // apply correction efficiency, STEP 10
	    if (fEfficiencyCorrection)
	    {
	      // with or without two track efficiency depending on if fTwoTrackEfficiencyCut is set
	      Bool_t twoTrackCut = (fTwoTrackEfficiencyCut > 0);
	      
	      fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepCorrected, tracks, pool2->GetEvent(jMix), 1.0 / pool2->GetCurrentNEvents(), (jMix == 0), twoTrackCut, bSign, fTwoTrackEfficiencyCut, kTRUE);
	    }
	  }
	}
	pool2->UpdatePool(CloneAndReduceTrackList(tracksCorrelate));
      }
      
      if (0 && !fReduceMemoryFootprint)
      {
        // make list of secondaries (matched with MC)
        TObjArray* tracksRecoMatchedSecondaries = new TObjArray;
        for (Int_t i=0; i<tracksRecoMatchedAll->GetEntriesFast(); i++)
          if (((AliAODMCParticle*)tracksRecoMatchedAll->At(i))->IsPhysicalPrimary() == kFALSE)
            tracksRecoMatchedSecondaries->Add(tracksRecoMatchedAll->At(i));
      
        // Study: Use only secondaries as trigger particles and plot the correlation vs. all particles; store in step 9
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy2, tracksRecoMatchedSecondaries, tracksRecoMatchedAll, weight);
        
        // Study: Use only primaries as trigger particles and plot the correlation vs. secondaries; store in step 8
        fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracksRecoMatchedPrim, tracksRecoMatchedSecondaries, weight);
      
        // plot delta phi vs process id of secondaries
        // trigger particles: primaries in 4 < pT < 10
        // associated particles: secondaries in 1 < pT < 10
        
        for (Int_t i=0; i<tracksRecoMatchedPrim->GetEntriesFast(); i++)
        {
          AliVParticle* triggerParticle = (AliVParticle*) tracksRecoMatchedPrim->At(i);
          
          if (triggerParticle->Pt() < 4 || triggerParticle->Pt() > 10)
            continue;
          
          for (Int_t j=0; j<tracksRecoMatchedSecondaries->GetEntriesFast(); j++)
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
  
      if (tracksCorrelateRecoMatchedPrim != tracksRecoMatchedPrim)
	delete tracksCorrelateRecoMatchedPrim;
      delete tracksRecoMatchedPrim;

      if (tracksCorrelateRecoMatchedAll != tracksRecoMatchedAll)
	delete tracksCorrelateRecoMatchedAll;
      delete tracksRecoMatchedAll;
      
      if (tracksCorrelate != tracks)
	delete tracksCorrelate;
      delete tracks;
    }
  }
  
  if (tracksMC != tracksCorrelateMC)
    delete tracksCorrelateMC;
  delete tracksMC;
}

//____________________________________________________________________
void  AliAnalysisTaskPhiCorrelations::AnalyseDataMode()
{

  // Run the analysis on DATA or MC to get raw distributions
 
  if ( fDebug > 3 ) AliInfo( " Processing event in Data mode ..." );

  ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(0);

  if (!fInputHandler)
    return;
    
  // skip not selected events here (the AOD is not updated for those)
  if (!fSkipTrigger && !(fInputHandler->IsEventSelected() & fSelectBit))
    return;
  
  // skip fast cluster events here if requested
  if (fSkipFastCluster && (fInputHandler->IsEventSelected() & AliVEvent::kFastOnly))
    return;

  // Support for ESD and AOD based analysis
  AliVEvent* inputEvent = fAOD;
  if (!inputEvent)
    inputEvent = fESD;

  Double_t centrality = 0;
  
  AliCentrality *centralityObj = 0;
  if (fCentralityMethod.Length() > 0)
  {
    if (fCentralityMethod == "ZNA_MANUAL")
    {
      Bool_t zna = kFALSE;
      for(Int_t j = 0; j < 4; ++j) {
	if (fESD->GetZDCData()->GetZDCTDCData(12,j) != 0) {
	  zna = kTRUE;
	}
      }

//       Printf("%d %f", zna, fZNAtower[0]);
      if (zna)
      {
	// code from Chiara O (23.10.12)
	const Double_t *fZNAtower = fESD->GetZDCData()->GetZN2TowerEnergy();
	Float_t znacut[4] = {681., 563., 413., 191.};
	
	if(fZNAtower[0]>znacut[0]) centrality = 1;
	else if(fZNAtower[0]>znacut[1]) centrality = 21;
	else if(fZNAtower[0]>znacut[2]) centrality = 41;
	else if(fZNAtower[0]>znacut[3]) centrality = 61;
	else centrality = 81;
      }
      else
	centrality = -1;
    }
    else if (fCentralityMethod == "TRACKS_MANUAL")
    {
      // for pp
      TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, -1, kTRUE);
      centrality = tracks->GetEntriesFast();
      if (centrality > 40)
	centrality = 41;
//       Printf("%d %f", tracks->GetEntriesFast(), centrality);
      delete tracks;
    }
    else
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

      if (fAOD)
      {
	// remove outliers
	if (centrality == 0)
	{
	  if (fAOD->GetVZEROData())
	  {
	    Float_t multV0 = 0;
	    for (Int_t i=0; i<64; i++)
	      multV0 += fAOD->GetVZEROData()->GetMultiplicity(i);
	    if (multV0 < 19500)
	    {
	      centrality = -1;
	      AliInfo("Rejecting event due to too small V0 multiplicity");
	    }
	  }
	}
      }
    }
    
    AliInfo(Form("Centrality is %f", centrality));
  }
  
  Float_t bSign = (inputEvent->GetMagneticField() > 0) ? 1 : -1;

  fHistos->SetRunNumber(inputEvent->GetRunNumber());
  
  // Fill the "event-counting-container", it is needed to get the number of events remaining after each event-selection cut
  fHistos->FillEvent(centrality, AliUEHist::kCFStepAll);
  
  // Trigger selection ************************************************
  if (!fSkipTrigger && !fAnalyseUE->TriggerSelection(fInputHandler)) return;
  
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
  
  // check for outlier in centrality vs number of tracks (rough constants extracted from correlation histgram)
  Bool_t reject = kFALSE;
  if (fRejectCentralityOutliers)
  {
    if (centrality > 40 && centrality <= 50 && tracks->GetEntriesFast() > 1160)
      reject = kTRUE;
    if (centrality > 50 && centrality <= 60 && tracks->GetEntriesFast() > 650)
      reject = kTRUE;
    if (centrality > 60 && centrality <= 70 && tracks->GetEntriesFast() > 370)
      reject = kTRUE;
    if (centrality > 70 && centrality <= 80 && tracks->GetEntriesFast() > 220)
      reject = kTRUE;
    if (centrality > 80 && centrality <= 90 && tracks->GetEntriesFast() > 130)
      reject = kTRUE;
    if (centrality > 90 && tracks->GetEntriesFast() > 75)
      reject = kTRUE;
  }
  
  if (reject)
  {
    AliInfo(Form("Rejecting event due to centrality vs tracks correlation: %f %d", centrality, tracks->GetEntriesFast()));
    fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
    delete tracks;
    return;
  }
  
  // correlate particles with...
  TObjArray* tracksCorrelate = 0;
  if (fSelectParticleSpecies != -1)
    tracksCorrelate = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, fSelectParticleSpecies, kTRUE);
  
  // reference multiplicity
  Int_t referenceMultiplicity = -1;
  if (fESD)
    referenceMultiplicity = AliESDtrackCuts::GetReferenceMultiplicity(fESD);
  else if (fAOD)
    referenceMultiplicity = tracks->GetEntriesFast(); // TODO to be replaced by the estimator once available in the AOD

  ((TH2F*) fListOfHistos->FindObject("referenceMultiplicity"))->Fill(centrality, referenceMultiplicity);
  
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
  {
    if (!fSkipStep6)
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracks, tracksCorrelate, weight, kTRUE, kFALSE, 0, 0.02, kTRUE);
    
    ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(1);
    
    if (fTwoTrackEfficiencyCut > 0)
      fHistos->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracks, tracksCorrelate, weight, kTRUE, kTRUE, bSign, fTwoTrackEfficiencyCut, kTRUE);
  }

  // fill second time with SPD centrality
  if (fCompareCentralities && centralityObj)
  {
    centrality = centralityObj->GetCentralityPercentile("CL1");
    if (centrality >= 0 && !fSkipStep6)
      fHistos->FillCorrelations(centrality, 2, AliUEHist::kCFStepReconstructed, tracks, tracksCorrelate, weight, kFALSE, kFALSE, 0, 0.02, kTRUE);
  }
  
  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = CloneAndReduceTrackList(tracks);
  delete tracks;
  
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
     
    if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5) 
    {
      
      Int_t nMix = pool->GetCurrentNEvents();
//       cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
      
      ((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
      ((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
      if (pool->IsReady())
	((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(3);
    
      // Fill mixed-event histos here  
      for (Int_t jMix=0; jMix<nMix; jMix++) 
      {
	TObjArray* bgTracks = pool->GetEvent(jMix);
	
	if (!fSkipStep6)
	  fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepReconstructed, tracksClone, bgTracks, 1.0 / nMix, (jMix == 0), kFALSE, 0, 0.02, kTRUE);

	if (fTwoTrackEfficiencyCut > 0)
	  fHistosMixed->FillCorrelations(centrality, zVtx, AliUEHist::kCFStepBiasStudy, tracksClone, bgTracks, 1.0 / nMix, (jMix == 0), kTRUE, bSign, fTwoTrackEfficiencyCut, kTRUE);
      }
    }
    
    // ownership is with the pool now
    if (tracksCorrelate)
    {
      pool->UpdatePool(CloneAndReduceTrackList(tracksCorrelate));
      delete tracksClone;
    }
    else
      pool->UpdatePool(tracksClone);
    //pool->PrintInfo();
  }
  else
  {
    delete tracksClone;
    if (tracksCorrelate)
      delete tracksCorrelate;
  }
}

TObjArray* AliAnalysisTaskPhiCorrelations::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliDPhiBasicParticle which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* particle = (AliVParticle*) tracks->At(i);
    tracksClone->Add(new AliDPhiBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge()));
  }
  
  return tracksClone;
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

//____________________________________________________________________
void AliAnalysisTaskPhiCorrelations::RemoveDuplicates(TObjArray* tracks)
{
  // remove particles with the same label
  
  Int_t before = tracks->GetEntriesFast();

  for (Int_t i=0; i<before; ++i) 
  {
    AliVParticle* part = (AliVParticle*) tracks->At(i);
    
    for (Int_t j=i+1; j<before; ++j) 
    {
      AliVParticle* part2 = (AliVParticle*) tracks->At(j);
      
      if (part->GetLabel() == part2->GetLabel())
      {
	Printf("Removing %d with label %d (duplicated in %d)", i, part->GetLabel(), j); part->Dump(); part2->Dump();
	TObject* object = tracks->RemoveAt(i);
	if (tracks->IsOwner())
	  delete object;
	break;
      }
    }
  }
 
  tracks->Compress();
  
  if (before > tracks->GetEntriesFast())
    AliInfo(Form("Reduced from %d to %d", before, tracks->GetEntriesFast())); 
}

void AliAnalysisTaskPhiCorrelations::CleanUp(TObjArray* tracks, TObject* mcObj, Int_t maxLabel)
{
  // calls RemoveInjectedSignals, RemoveWeakDecays and RemoveDuplicates
  
  if (fInjectedSignals)
    fAnalyseUE->RemoveInjectedSignals(tracks, mcObj, maxLabel);
  if (fRemoveWeakDecays)
    fAnalyseUE->RemoveWeakDecays(tracks, mcObj);
  if (fRemoveDuplicates)
    RemoveDuplicates(tracks);
}
