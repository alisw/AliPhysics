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
//
// The analysis task:
// Filling an AliCFContainer with the quantities pt, eta and phi
// for tracks which survivied the particle cuts (MC resp. ESD tracks)
// Track selection is done using the AliHFE package
// 
// Author:
//  Raphaelle Bailhache <R.Bailhache@gsi.de>
//  Markus Fasel <M.Fasel@gsi.de>
//  Matus Kalisky <matus.kalisky@cern.ch>
//  MinJung Kweon <minjung@physi.uni-heidelberg.de>
//
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>

#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliVVertex.h"

#include "AliHFEpid.h"
#include "AliHFEcollection.h"
#include "AliHFEcuts.h"
#include "AliHFEmcQA.h"
#include "AliHFEpairs.h"
#include "AliHFEpostAnalysis.h"
#include "AliHFEsecVtxs.h"
#include "AliHFEsecVtx.h"
#include "AliHFEelecbackground.h"
#include "AliHFEtools.h"
#include "AliAnalysisTaskHFE.h"

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE():
  AliAnalysisTaskSE("PID efficiency Analysis")
  , fQAlevel(0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fPlugins(0)
  , fCFM(NULL)
  , fCorrelation(NULL)
  , fPIDperformance(NULL)
  , fSignalToBackgroundMC(NULL)
  , fPID(NULL)
  , fCuts(NULL)
  , fSecVtx(NULL)
  , fElecBackGround(NULL)
  , fMCQA(NULL)
  , fNEvents(NULL)
  , fNElectronTracksEvent(NULL)
  , fQA(NULL)
  , fOutput(NULL)
  , fHistMCQA(NULL)
  , fHistSECVTX(NULL)
  , fHistELECBACKGROUND(NULL)
//  , fQAcoll(0x0)
{
  //
  // Dummy constructor
  //
  DefineOutput(1, TH1I::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
//  DefineOutput(4, TList::Class());

  // Initialize cuts
  fPID = new AliHFEpid;
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const char * name):
  AliAnalysisTaskSE(name)
  , fQAlevel(0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fPlugins(0)
  , fCFM(NULL)
  , fCorrelation(NULL)
  , fPIDperformance(NULL)
  , fSignalToBackgroundMC(NULL)
  , fPID(NULL)
  , fCuts(NULL)
  , fSecVtx(NULL)
  , fElecBackGround(NULL)
  , fMCQA(NULL)
  , fNEvents(NULL)
  , fNElectronTracksEvent(NULL)
  , fQA(NULL)
  , fOutput(NULL)
  , fHistMCQA(NULL)
  , fHistSECVTX(NULL)
  , fHistELECBACKGROUND(NULL)
//  , fQAcoll(0x0)
{
  //
  // Default constructor
  // 
  DefineOutput(1, TH1I::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
//  DefineOutput(4, TList::Class());

  // Initialize cuts
  fPID = new AliHFEpid;
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref):
  AliAnalysisTaskSE(ref)
  , fQAlevel(ref.fQAlevel)
  , fPIDdetectors(ref.fPIDdetectors)
  , fPIDstrategy(ref.fPIDstrategy)
  , fPlugins(ref.fPlugins)
  , fCFM(ref.fCFM)
  , fCorrelation(ref.fCorrelation)
  , fPIDperformance(ref.fPIDperformance)
  , fSignalToBackgroundMC(ref.fSignalToBackgroundMC)
  , fPID(ref.fPID)
  , fCuts(ref.fCuts)
  , fSecVtx(ref.fSecVtx)
  , fElecBackGround(ref.fElecBackGround)
  , fMCQA(ref.fMCQA)
  , fNEvents(ref.fNEvents)
  , fNElectronTracksEvent(ref.fNElectronTracksEvent)
  , fQA(ref.fQA)
  , fOutput(ref.fOutput)
  , fHistMCQA(ref.fHistMCQA)
  , fHistSECVTX(ref.fHistSECVTX)
  , fHistELECBACKGROUND(ref.fHistELECBACKGROUND)
//  , fQAcoll(ref.fQAcoll)
{
  //
  // Copy Constructor
  //
}

//____________________________________________________________
AliAnalysisTaskHFE &AliAnalysisTaskHFE::operator=(const AliAnalysisTaskHFE &ref){
  //
  // Assignment operator
  //
  if(this == &ref) return *this;
  AliAnalysisTask::operator=(ref);
  fQAlevel = ref.fQAlevel;
  fPIDdetectors = ref.fPIDdetectors;
  fPIDstrategy = ref.fPIDstrategy;
  fPlugins = ref.fPlugins;
  fCFM = ref.fCFM;
  fCorrelation = ref.fCorrelation;
  fPIDperformance = ref.fPIDperformance;
  fSignalToBackgroundMC = ref.fSignalToBackgroundMC;
  fPID = ref.fPID;
  fCuts = ref.fCuts;
  fSecVtx = ref.fSecVtx;
  fElecBackGround = ref.fElecBackGround;
  fMCQA = ref.fMCQA;
  fNEvents = ref.fNEvents;
  fNElectronTracksEvent = ref.fNElectronTracksEvent;
  fQA = ref.fQA;
  fOutput = ref.fOutput;
  fHistMCQA = ref.fHistMCQA;
  fHistSECVTX = ref.fHistSECVTX;
  fHistELECBACKGROUND = ref.fHistELECBACKGROUND;
  
//  fQAcoll = ref.fQAcoll;
  return *this;
}

//____________________________________________________________
AliAnalysisTaskHFE::~AliAnalysisTaskHFE(){
  //
  // Destructor
  //
  if(fPID) delete fPID;
  if(fQA){
    fQA->Clear();
    delete fQA;
  }
  if(fOutput){ 
    fOutput->Clear();
    delete fOutput;
  }
  if(fHistMCQA){
    fHistMCQA->Clear();
    delete fHistMCQA;
  }
  if(fHistSECVTX){
    fHistSECVTX->Clear();
    delete fHistSECVTX;
  }
  if(fHistELECBACKGROUND){
    fHistELECBACKGROUND->Clear();
    delete fHistELECBACKGROUND;
  }
  if(fSecVtx) delete fSecVtx;
  if(fElecBackGround) delete fElecBackGround;
  if(fMCQA) delete fMCQA;
  if(fNEvents) delete fNEvents;
  if(fCorrelation){
    fCorrelation->Clear();
    delete fCorrelation;
  }
  if(fPIDperformance) delete fPIDperformance;
  if(fSignalToBackgroundMC) delete fSignalToBackgroundMC;
//  if(fQAcoll) delete fQAcoll;
}

//____________________________________________________________
void AliAnalysisTaskHFE::UserCreateOutputObjects(){
  //
  // Creating output container and output objects
  // Here we also Initialize the correction framework container and 
  // the objects for
  // - PID
  // - MC QA
  // - SecVtx
  // QA histograms are created if requested
  // Called once per worker
  //
  AliDebug(3, "Creating Output Objects");
  // Automatic determination of the analysis mode
  AliVEventHandler *inputHandler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
    SetAODAnalysis();
  } else {
    SetESDAnalysis();
    if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())
      SetHasMCData();
  }
  printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");
  printf("MC Data available %s\n", HasMCData() ? "Yes" : "No");

  // example how to use the AliHFEcollection
  //fQAcoll = new AliHFEcollection("fQAcoll", "QA");
  //fQAcoll->CreateTH1F("fNevents", "Number of Events in the Analysis", 2, 0, 2);
  //fQAcoll->CreateProfile("fNtrdclusters", "Number of TRD clusters as function of momentum; p[GeV/c]", 20, 0, 20);

  fNEvents = new TH1I("nEvents", "Number of Events in the Analysis", 2, 0, 2); // Number of Events neccessary for the analysis and not a QA histogram
  fNElectronTracksEvent = new TH1I("nElectronTracksEvent", "Number of Electron Candidates", 100, 0, 100);
  // First Step: TRD alone
  if(!fQA) fQA = new TList;
  fQA->AddAt(new TProfile("conr", "Electron PID contamination", 20, 0, 20), 0);
  fQA->AddAt(new TH1F("alpha_rec", "Alpha from reconstructed tracks with TRD hits", 36, -TMath::Pi(), TMath::Pi()), 1);
  fQA->AddAt(new TH1F("alpha_sim", "Alpha from simulated electron tracks", 36, -TMath::Pi(), TMath::Pi()), 2);
  fQA->AddAt(new TH1F("nElectron", "Number of electrons", 100, 0, 100), 3);
  fQA->AddAt(new TProfile("pidquality", "TRD PID quality as function of momentum", 20, 0, 20), 4);
  fQA->AddAt(new TProfile("ntrdclusters", "Number of TRD clusters as function of momentum", 20, 0, 20), 5);
  fQA->AddAt(new TH1F("chi2TRD","#chi2 per TRD cluster", 20, 0, 20), 6);
  fQA->AddAt(new TH1I("mccharge", "MC Charge", 200, -100, 100), 7);

  if(!fOutput) fOutput = new TList;
  // Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  MakeParticleContainer();
  MakeEventContainer();
  // Temporary fix: Initialize particle cuts with NULL
  for(Int_t istep = 0; istep < fCFM->GetParticleContainer()->GetNStep(); istep++)
    fCFM->SetParticleCutsList(istep, NULL);
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  if(IsAODanalysis()) fCuts->SetAOD();
  fCuts->Initialize(fCFM);
  if(fCuts->IsInDebugMode()) fQA->Add(fCuts->GetQAhistograms());
 
  // add output objects to the List
  fOutput->AddAt(fCFM->GetParticleContainer(), 0);
  fOutput->AddAt(fCFM->GetEventContainer(), 1);
  fOutput->AddAt(fCorrelation, 2);
  fOutput->AddAt(fPIDperformance, 3);
  fOutput->AddAt(fSignalToBackgroundMC, 4);
  fOutput->AddAt(fNElectronTracksEvent, 5);

  // Initialize PID
  if(IsQAOn(kPIDqa)){
    AliInfo("PID QA switched on");
    //fPID->SetDebugLevel(2);
    fPID->SetQAOn();
    fQA->Add(fPID->GetQAhistograms());
  }
  fPID->SetHasMCData(HasMCData());
  if(!fPIDdetectors.Length() && ! fPIDstrategy) AddPIDdetector("TPC");
  if(fPIDstrategy)
    fPID->InitializePID(Form("Strategy%d", fPIDstrategy));
  else
    fPID->InitializePID(fPIDdetectors.Data());     // Only restrictions to TPC allowed 

  // mcQA----------------------------------
  if (HasMCData() && IsQAOn(kMCqa)) {
    AliInfo("MC QA on");
    if(!fMCQA) fMCQA = new AliHFEmcQA;
    if(!fHistMCQA) fHistMCQA = new TList();
    fHistMCQA->SetName("MCqa");
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,0,"mcqa_");               // create histograms for charm
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,0,"mcqa_");              // create histograms for beauty
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,1,"mcqa_barrel_");        // create histograms for charm 
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,1,"mcqa_barrel_");       // create histograms for beauty
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,2,"mcqa_unitY_");         // create histograms for charm 
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,2,"mcqa_unitY_");        // create histograms for beauty
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,3,"mcqa_reccut_");        // create histograms for charm 
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,3,"mcqa_reccut_");       // create histograms for beauty
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,4,"mcqa_recpidcut_");     // create histograms for charm 
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,4,"mcqa_recpidcut_");    // create histograms for beauty
    TIter next(gDirectory->GetList());
    TObject *obj;
    int counter = 0;
    TString objname;
    while ((obj = next.Next())) {
      objname = obj->GetName();
      TObjArray *toks = objname.Tokenize("_");
      if (toks->GetEntriesFast()){
        TObjString *fpart = (TObjString *)(toks->UncheckedAt(0));
        if ((fpart->String()).CompareTo("mcqa") == 0) fHistMCQA->AddAt(obj, counter++);
      }
    }
    fQA->Add(fHistMCQA);
  } 

  // secvtx----------------------------------
  if (GetPlugin(kSecVtx)) {
    AliInfo("Secondary Vertex Analysis on");
    fSecVtx = new AliHFEsecVtx;
    fSecVtx->SetHasMCData(HasMCData());

    if(!fHistSECVTX) fHistSECVTX = new TList();
    fSecVtx->CreateHistograms(fHistSECVTX);
    fOutput->Add(fHistSECVTX);
  }
  
  // background----------------------------------
  if (GetPlugin(kIsElecBackGround)) {
    AliInfo("Electron BackGround Analysis on");
    if(!fElecBackGround){
      AliWarning("ElecBackGround not available. Default elecbackground will be used");
      fElecBackGround = new AliHFEelecbackground;
    }
    fElecBackGround->SetHasMCData(HasMCData());

    if(!fHistELECBACKGROUND) fHistELECBACKGROUND = new TList();
    fElecBackGround->CreateHistograms(fHistELECBACKGROUND);
    fOutput->Add(fHistELECBACKGROUND);
  }  
}

//____________________________________________________________
void AliAnalysisTaskHFE::UserExec(Option_t *){
  //
  // Run the analysis
  // 
  AliDebug(3, "Starting Single Event Analysis");
  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return;
  }
  if(HasMCData()){
    AliDebug(4, Form("MC Event: %p", fMCEvent));
    if(!fMCEvent){
      AliError("No MC Event, but MC Data required");
      return;
    }
  }
  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }

  if(HasMCData()) ProcessMC();  // Run the MC loop + MC QA in case MC Data are available

  if(IsAODanalysis()) ProcessAOD();
  else ProcessESD();
  // Done!!!
  PostData(1, fNEvents);
  PostData(2, fOutput);
  PostData(3, fQA);
//  PostData(4, fQAcoll->GetList());
}

//____________________________________________________________
void AliAnalysisTaskHFE::Terminate(Option_t *){
  //
  // Terminate not implemented at the moment
  //
  if(GetPlugin(kPostProcess)){
    fOutput = dynamic_cast<TList *>(GetOutputData(2));
    if(!fOutput){
      AliError("Results not available");
      return;
    }
    AliHFEpostAnalysis postanalysis;
    postanalysis.SetResults(fOutput);
    if(HasMCData())postanalysis.DrawMCSignal2Background();
    postanalysis.DrawEfficiency();
    postanalysis.DrawPIDperformance();
    postanalysis.DrawCutEfficiency();

    if (GetPlugin(kIsElecBackGround)) {
      AliHFEelecbackground elecBackGround;
      TList *oe = 0x0;
      if(!(oe = (TList*)dynamic_cast<TList *>(fOutput->FindObject("HFEelecbackground")))){
	return;
      }
      elecBackGround.Load(oe);
      elecBackGround.Plot();
      elecBackGround.PostProcess();      
    }
  }
}

//____________________________________________________________
void AliAnalysisTaskHFE::ProcessMC(){
  //
  // Runs the MC Loop (filling the container for the MC Cut Steps with the observables pt, eta and phi)
  // In case MC QA is on also MC QA loop is done
  //
  AliDebug(3, "Processing MC Information");
  Double_t nContrib = 0;
  const AliVVertex *pVertex = fMCEvent->GetPrimaryVertex();
  if(pVertex) nContrib = pVertex->GetNContributors();
  if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepGenerated, fMCEvent)) return;
  fCFM->GetEventContainer()->Fill(&nContrib,AliHFEcuts::kEventStepGenerated);
  Int_t nElectrons = 0;
  if(IsESDanalysis()){
    if (HasMCData() && IsQAOn(kMCqa)) {
      AliDebug(2, "Running MC QA");

      if(fMCEvent->Stack()){
        fMCQA->SetStack(fMCEvent->Stack());
        fMCQA->SetGenEventHeader(fMCEvent->GenEventHeader());
        fMCQA->Init();

        Int_t nMCTracks = fMCEvent->Stack()->GetNtrack();

        // loop over all tracks for decayed electrons
        for (Int_t igen = 0; igen < nMCTracks; igen++){
          TParticle* mcpart = fMCEvent->Stack()->Particle(igen);
          fMCQA->GetQuarkKine(mcpart, igen, AliHFEmcQA::kCharm);
          fMCQA->GetQuarkKine(mcpart, igen, AliHFEmcQA::kBeauty);
          fMCQA->GetHadronKine(mcpart, AliHFEmcQA::kCharm);
          fMCQA->GetHadronKine(mcpart, AliHFEmcQA::kBeauty);
          fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 0); // no accept cut
          fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 0); // no accept cut
          if (TMath::Abs(mcpart->Eta()) < 0.9) {
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 1); // accept |eta|<0.9
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 1); // accept |eta|<0.9
          }
          if (TMath::Abs(AliHFEtools::GetRapidity(mcpart)) < 0.5) {
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 2); // accept |y|<0.5
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 2); // accept |y|<0.5
          }
        }
        fMCQA->EndOfEventAna(AliHFEmcQA::kCharm);
        fMCQA->EndOfEventAna(AliHFEmcQA::kBeauty);
      }

    } // end of MC QA loop
    // -----------------------------------------------------------------
    fCFM->SetMCEventInfo(fMCEvent);
    // fCFM->CheckEventCuts(AliCFManager::kEvtRecCuts, fESD);
  } else {
    fCFM->SetMCEventInfo(fInputEvent);
  }
  // Run MC loop
  AliVParticle *mctrack = NULL;
  for(Int_t imc = 0; imc <fMCEvent->GetNumberOfTracks(); imc++){
    if(!(mctrack = fMCEvent->GetTrack(imc))) continue;
    if(ProcessMCtrack(mctrack)) nElectrons++;
  }

  // fCFM->CheckEventCuts(AliCFManager::kEvtRecCuts, fESD);
  (dynamic_cast<TH1F *>(fQA->At(3)))->Fill(nElectrons);
}

//____________________________________________________________
void AliAnalysisTaskHFE::ProcessESD(){
  //
  // Run Analysis of reconstructed event in ESD Mode
  // Loop over Tracks, filter according cut steps defined in AliHFEcuts
  //
  AliDebug(3, "Processing ESD Event");
  Double_t nContrib = fInputEvent->GetPrimaryVertex()->GetNContributors();
  if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fInputEvent)) return;
  fCFM->GetEventContainer()->Fill(&nContrib, AliHFEcuts::kEventStepReconstructed);
  AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!fESD){
    AliError("ESD Event required for ESD Analysis")
    return;
  }
  if (GetPlugin(kIsElecBackGround)) { 
    fElecBackGround->SetEvent(fESD);
  }
  if (GetPlugin(kSecVtx)) {
    fSecVtx->SetEvent(fESD);
    fSecVtx->GetPrimaryCondition();
  }

  if(HasMCData()){
    if (GetPlugin(kSecVtx)) { 
      if(fMCEvent->Stack()) fSecVtx->SetStack(fMCEvent->Stack());
    }
    if (GetPlugin(kIsElecBackGround)) { 
      fElecBackGround->SetMCEvent(fMCEvent);
    }
  }


  Double_t container[8];
  memset(container, 0, sizeof(Double_t) * 8);
  // container for the output THnSparse
  Double_t dataE[5]; // [pT, eta, Phi, type, 'C' or 'B']
  Int_t nElectronCandidates = 0;
  AliESDtrack *track = NULL, *htrack = NULL;
  AliMCParticle *mctrack = NULL;
  TParticle* mctrack4QA = NULL;
  Int_t pid = 0;
  // For double counted tracks
  LabelContainer cont(fESD->GetNumberOfTracks());
  Bool_t alreadyseen = kFALSE;

  Bool_t signal = kTRUE;

  fCFM->SetRecEventInfo(fESD);
  // Electron background analysis 
  if (GetPlugin(kIsElecBackGround)) {
    
    AliDebug(2, "Running BackGround Analysis");
    
    fElecBackGround->Reset();
    
  } // end of electron background analysis
  //
  // Loop ESD
  //
  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    
    track = fESD->GetTrack(itrack);
          
    container[0] = track->Pt();
    container[1] = track->Eta();
    container[2] = track->Phi();
    container[3] = track->Charge();

    dataE[0] = track->Pt();
    dataE[1] = track->Eta();
    dataE[2] = track->Phi();
    dataE[3] = track->Charge();
    dataE[4] = -1;
    dataE[5] = -1;

    signal = kTRUE;

    // Fill step without any cut
          
    if(HasMCData()){
      // Check if it is electrons near the vertex
      if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(track->GetLabel()))))) continue;
      mctrack4QA = mctrack->Particle();//fMCEvent->Stack()->Particle(TMath::Abs(track->GetLabel()));

      container[4] = mctrack->Pt();
      container[5] = mctrack->Eta();
      container[6] = mctrack->Phi();
      container[7] = mctrack->Charge()/3.;
    
      if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE;
    }
    if(signal) {
      alreadyseen = cont.Find(TMath::Abs(track->GetLabel()));
      cont.Append(TMath::Abs(track->GetLabel()));
      
      fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepRecNoCut);
      fCFM->GetParticleContainer()->Fill(&container[0], AliHFEcuts::kStepRecNoCut + 2*AliHFEcuts::kNcutStepsESDtrack);
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepRecNoCut + AliHFEcuts::kNcutStepsESDtrack);
      }
    }

    // RecKine: ITSTPC cuts  
    if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track, container, signal, alreadyseen)) continue;
    
    // Check TRD criterions (outside the correction framework)
    if(track->GetTRDncls()){
      (dynamic_cast<TH1F *>(fQA->At(6)))->Fill(track->GetTRDchi2()/track->GetTRDncls());
      (dynamic_cast<TH1F *>(fQA->At(1)))->Fill(track->GetAlpha());    // Check the acceptance without tight cuts
      (dynamic_cast<TProfile *>(fQA->At(4)))->Fill(container[0], track->GetTRDpidQuality());
      (dynamic_cast<TProfile *>(fQA->At(5)))->Fill(container[0], track->GetTRDncls());
      //fQAcoll->Fill("fNtrdclusters", container[0], track->GetTRDncls());
    }

    
    // RecPrim
    if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track, container, signal, alreadyseen)) continue;

    // HFEcuts: ITS layers cuts
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track, container, signal, alreadyseen)) continue;

    // HFEcuts: Nb of tracklets TRD0
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTRD, track, container, signal, alreadyseen)) continue;
    if(signal) {
      // dimensions 3&4&5 : pt,eta,phi (MC)
      ((THnSparseF *)fCorrelation->At(0))->Fill(container);
    }

    if(HasMCData() && IsQAOn(kMCqa)) {
      // mc qa for after the reconstruction cuts  
      AliDebug(2, "Running MC QA");
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 3);  // charm
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 3); // beauty 
    } 

    // track accepted, do PID
    AliHFEpidObject hfetrack;
    hfetrack.fAnalysisType = AliHFEpidObject::kESDanalysis;
    hfetrack.fRecTrack = track;
    if(HasMCData()) hfetrack.fMCtrack = mctrack;
    if(!fPID->IsSelected(&hfetrack)) continue;
    nElectronCandidates++;

    if (HasMCData() && IsQAOn(kMCqa)) {
      // mc qa for after the reconstruction and pid cuts  
      AliDebug(2, "Running MC QA");
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 4);  // charm
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 4); // beauty 
    } 

    // Fill Containers
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepPID + 2*AliHFEcuts::kNcutStepsESDtrack);
      fCFM->GetParticleContainer()->Fill(&container[4], AliHFEcuts::kStepPID);
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[4], (AliHFEcuts::kStepPID + (AliHFEcuts::kNcutStepsESDtrack)));
      }
      // dimensions 3&4&5 : pt,eta,phi (MC)
      ((THnSparseF *)fCorrelation->At(1))->Fill(container);
    }

    if(GetPlugin(kSecVtx) && fMCEvent->Stack()) {
      AliDebug(2, "Running Secondary Vertex Analysis");
      if(track->Pt()>1.0){
        fSecVtx->InitHFEpairs();
        fSecVtx->InitHFEsecvtxs();
        for(Int_t jtrack = 0; jtrack < fESD->GetNumberOfTracks(); jtrack++){
          htrack = fESD->GetTrack(jtrack);
          if ( itrack == jtrack ) continue; // since it is for tagging single electron, don't need additional condition 
          if (htrack->Pt()<1.0) continue;
          if (!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC, htrack)) continue;
          if (!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim, htrack)) continue;
          fSecVtx->PairAnalysis(track, htrack, jtrack); // e-h pairing
        }
        /*for(int ip=0; ip<fSecVtx->HFEpairs()->GetEntriesFast(); ip++){
          if(HasMCData()){
            AliHFEpairs *pair = (AliHFEpairs*) (fSecVtx->HFEpairs()->UncheckedAt(ip));
            if(!(pair->GetPairCode()>1. && pair->GetPairCode()<4.))  // apply various cuts
              fSecVtx->HFEpairs()->RemoveAt(ip);
          }
        }*/
        fSecVtx->HFEpairs()->Compress();
        fSecVtx->RunSECVTX(track); // secondary vertexing with e,h1,h2,.. tracks
        for(int ip=0; ip<fSecVtx->HFEsecvtxs()->GetEntriesFast(); ip++){
          AliHFEsecVtxs *secvtx=0x0;
          secvtx = (AliHFEsecVtxs*) (fSecVtx->HFEsecvtxs()->UncheckedAt(ip));
          // here you apply cuts, then if it doesn't pass the cut, remove it from the fSecVtx->HFEsecvtxs() 
        }
        fSecVtx->DeleteHFEpairs();
        fSecVtx->DeleteHFEsecvtxs();
      }
    }

    if(HasMCData()){
      // Track selected: distinguish between true and fake
      AliDebug(1, Form("Candidate Selected, filling THnSparse, PID: %d\n", mctrack->Particle()->GetPdgCode()));
      if((pid = TMath::Abs(mctrack->Particle()->GetPdgCode())) == 11){
        Int_t type = IsSignalElectron(track);
        AliDebug(1, Form("Type: %d\n", type));
        if(type){
	        dataE[5] = type; // beauty[1] or charm[2]
	        dataE[4] = 2;  // signal electron
        }
        else{
	        dataE[4] = 1; // not a signal electron
	        dataE[5] = 0;
        }
      } 
      else {
        // Fill THnSparse with the information for Fake Electrons
        dataE[4] = 0;
        dataE[5] = 0;
      }
      // fill the performance THnSparse, if the mc origin could be defined
      if(dataE[4] > -1){
        AliDebug(1, Form("Entries: [%.3f|%.3f|%.3f|%f|%f|%f]\n", dataE[0],dataE[1],dataE[2],dataE[3],dataE[4],dataE[5]));
        fPIDperformance->Fill(dataE);
      }
    }
    // Electron background analysis 
    if (GetPlugin(kIsElecBackGround)) {
      
      AliDebug(2, "Running BackGround Analysis");
      
      for(Int_t jtrack = 0; jtrack < fESD->GetNumberOfTracks(); jtrack++){
	      htrack = fESD->GetTrack(jtrack);
      	if ( itrack == jtrack ) continue;  
      	fElecBackGround->PairAnalysis(track, htrack); 
      }
    } // end of electron background analysis

  }
  fNEvents->Fill(1);
  //fQAcoll->Fill("fNevents", 1);
  fNElectronTracksEvent->Fill(nElectronCandidates);
}

//____________________________________________________________
void AliAnalysisTaskHFE::ProcessAOD(){
  //
  // Run Analysis in AOD Mode
  // Function is still in development
  //
  AliDebug(3, "Processing AOD Event");
  Double_t nContrib = fInputEvent->GetPrimaryVertex()->GetNContributors();
  if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fInputEvent)) return;
  fCFM->GetEventContainer()->Fill(&nContrib,AliHFEcuts::kEventStepReconstructed);
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
  if(!fAOD){
    AliError("AOD Event required for AOD Analysis")
    return;
  }
 
  AliAODTrack *track = NULL;
  AliAODMCParticle *mctrack = NULL;
  Double_t container[8]; memset(container, 0, sizeof(Double_t) * 8);
  Double_t dataE[6]; // [pT, eta, Phi, Charge, type, 'C' or 'B']
  Int_t nElectronCandidates = 0;
  Int_t pid;
  for(Int_t itrack = 0; itrack < fAOD->GetNumberOfTracks(); itrack++){
    track = fAOD->GetTrack(itrack);
    if(!track) continue;
    if(track->GetFlags() != 1<<4) continue;  // Only process AOD tracks where the HFE is set

    container[0] = track->Pt();
    container[1] = track->Eta();
    container[2] = track->Phi();
    container[3] = track->Charge();

    dataE[0] = track->Pt();
    dataE[1] = track->Eta();
    dataE[2] = track->Phi();
    dataE[3] = track->Charge();
    dataE[4] = -1;
    dataE[5] = -1;
    
    if(HasMCData()){
      Int_t label = TMath::Abs(track->GetLabel());
      if(label){
        mctrack = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(label));
        container[4] = mctrack->Pt();
        container[5] = mctrack->Eta();
        container[6] = mctrack->Phi();
        container[7] = mctrack->Charge();
      }
    }
    // track accepted, do PID
    AliHFEpidObject hfetrack;
    hfetrack.fAnalysisType = AliHFEpidObject::kAODanalysis;
    hfetrack.fRecTrack = track;
    if(HasMCData()) hfetrack.fMCtrack = mctrack;
    //if(!fPID->IsSelected(&hfetrack)) continue;    // we will do PID here as soon as possible
    // Particle identified - Fill CF Container
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepPID + 2*AliHFEcuts::kNcutStepsESDtrack);
    nElectronCandidates++;    
    if(HasMCData()){
      // Track selected: distinguish between true and fake
      AliDebug(1, Form("Candidate Selected, filling THnSparse, PID: %d\n", mctrack->GetPdgCode()));
      if((pid = TMath::Abs(mctrack->GetPdgCode())) == 11){
        Int_t type = IsSignalElectron(track);
        AliDebug(1, Form("Type: %d\n", type));
        if(type){
	        dataE[5] = type; // beauty[1] or charm[2]
	        dataE[4] = 2;  // signal electron
        }
        else{
	        dataE[4] = 1; // not a signal electron
	        dataE[5] = 0;
        }
      } 
      else {
        // Fill THnSparse with the information for Fake Electrons
        dataE[4] = 0;
        dataE[5] = 0;
      }
      // fill the performance THnSparse, if the mc origin could be defined
      if(dataE[4] > -1){
        AliDebug(1, Form("Entries: [%.3f|%.3f|%.3f|%f|%f|%f]\n", dataE[0],dataE[1],dataE[2],dataE[3],dataE[4],dataE[5]));
        fPIDperformance->Fill(dataE);
      }
    }
  }
  fNEvents->Fill(1);
  fNElectronTracksEvent->Fill(nElectronCandidates);
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFE::ProcessMCtrack(AliVParticle *track){
  //
  // Filter the Monte Carlo Track
  // Additionally Fill a THnSparse for Signal To Background Studies
  // Works for AOD and MC analysis Type
  //
  Double_t container[4], signalContainer[6];
  Double_t vertex[3]; // Production vertex cut to mask gammas which are NOT supposed to have hits in the first ITS layer(s)
  if(IsESDanalysis()){
    AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(track);
    container[0] = mctrack->Pt();
    container[1] = mctrack->Eta();
    container[2] = mctrack->Phi();
    container[3] = mctrack->Charge()/3;

    signalContainer[0] = mctrack->Pt();
    signalContainer[1] = mctrack->Eta();
    signalContainer[2] = mctrack->Phi();
    signalContainer[3] = mctrack->Charge()/3;

    vertex[0] = mctrack->Particle()->Vx();
    vertex[1] = mctrack->Particle()->Vy();
  } else {
    AliAODMCParticle *aodmctrack = dynamic_cast<AliAODMCParticle *>(track);
    container[0] = aodmctrack->Pt();
    container[1] = aodmctrack->Eta();
    container[2] = aodmctrack->Phi();
    container[3] = aodmctrack->Charge()/3;

    signalContainer[0] = aodmctrack->Pt();
    signalContainer[1] = aodmctrack->Eta();
    signalContainer[2] = aodmctrack->Phi();
    signalContainer[3] = aodmctrack->Charge()/3;

    aodmctrack->XvYvZv(vertex);
  }
  if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, track)) return kFALSE;
  TH1 *test = dynamic_cast<TH1I*>(fQA->FindObject("mccharge"));
  test->Fill(signalContainer[3]);
 fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCGenerated);
  if((signalContainer[4] = static_cast<Double_t >(IsSignalElectron(track))) > 1e-3) fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCsignal);
  signalContainer[5] = 0;
  // apply cut on the sqrt of the production vertex
  Double_t radVertex = TMath::Sqrt(vertex[0]*vertex[0] + vertex[1] * vertex[1]);
  if(radVertex < 3.5){
    // Within first ITS layer(2) -> Background we cannot reject by ITS cut, let it pass
    signalContainer[5] = 1;
  } else if (radVertex < 7.5){
    signalContainer[5] = 2;
  }
  fSignalToBackgroundMC->Fill(signalContainer);
  (dynamic_cast<TH1F *>(fQA->At(2)))->Fill(container[2] - TMath::Pi());
  //if(IsESDanalysis()){
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCInAcceptance, track)) return kFALSE;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCInAcceptance);
  //}
  return kTRUE;
}

//____________________________________________________________
void AliAnalysisTaskHFE::MakeEventContainer(){
  //
  // Create the event container for the correction framework and link it
  //
  const Int_t kNvar = 1;  // number of variables on the grid: number of tracks per event
  const Double_t kNTrackBound[2] = {-0.5, 200.5};
  const Int_t kNBins = 201;

  AliCFContainer *evCont = new AliCFContainer("eventContainer", "Container for events", AliHFEcuts::kNcutStepsEvent, kNvar, &kNBins);

  Double_t *trackBins = AliHFEtools::MakeLinearBinning(kNBins, kNTrackBound[0], kNTrackBound[1]);
  evCont->SetBinLimits(0,trackBins);
  delete[] trackBins;

  fCFM->SetEventContainer(evCont);
}

//____________________________________________________________
void AliAnalysisTaskHFE::MakeParticleContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
  //
  const Int_t kNvar   = 4 ; //number of variables on the grid:pt,eta, phi, charge
  const Double_t kPtbound[2] = {0.1, 10.};
  const Double_t kEtabound[2] = {-0.8, 0.8};
  const Double_t kPhibound[2] = {0., 2. * TMath::Pi()};

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 40; // bins in pt
  iBin[1] =  8; // bins in eta 
  iBin[2] = 18; // bins in phi
  iBin[3] =  2; // bins in charge

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  binEdges[0] = AliHFEtools::MakeLogarithmicBinning(iBin[0], kPtbound[0], kPtbound[1]);
  binEdges[1] = AliHFEtools::MakeLinearBinning(iBin[1], kEtabound[0], kEtabound[1]);
  binEdges[2] = AliHFEtools::MakeLinearBinning(iBin[2], kPhibound[0], kPhibound[1]);
  binEdges[3] = AliHFEtools::MakeLinearBinning(iBin[3], -1.1, 1.1); // Numeric precision

  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("trackContainer", "Container for tracks", (AliHFEcuts::kNcutStepsTrack + 2*AliHFEcuts::kNcutStepsESDtrack), kNvar, iBin);

  //setting the bin limits
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    container -> SetBinLimits(ivar, binEdges[ivar]);
  fCFM->SetParticleContainer(container);

  //create correlation matrix for unfolding
  Int_t thnDim[2*kNvar];
  for (int k=0; k<kNvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = iBin[k];
    thnDim[k+kNvar] = iBin[k];
  }

  if(!fCorrelation) fCorrelation = new TList();
  fCorrelation->SetName("correlation");

  THnSparseF *correlation0 = new THnSparseF("correlationstepbeforePID","THnSparse with correlations",2*kNvar,thnDim);
  THnSparseF *correlation1 = new THnSparseF("correlationstepafterPID","THnSparse with correlations",2*kNvar,thnDim);
  for (int k=0; k<kNvar; k++) {
    correlation0->SetBinEdges(k,binEdges[k]);
    correlation0->SetBinEdges(k+kNvar,binEdges[k]);
    correlation1->SetBinEdges(k,binEdges[k]);
    correlation1->SetBinEdges(k+kNvar,binEdges[k]);
  }
  correlation0->Sumw2();
  correlation1->Sumw2();
  
  fCorrelation->AddAt(correlation0,0);
  fCorrelation->AddAt(correlation1,1);
  
  // Add a histogram for Fake electrons
  const Int_t nDim=6;
  Int_t nBin[nDim] = {40, 8, 18, 2, 3, 3};
  Double_t* binEdges2[nDim];

  //values for bin lower bounds
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    binEdges2[ivar] = binEdges[ivar];
  binEdges2[4] = AliHFEtools::MakeLinearBinning(nBin[4], 0, nBin[4]);
  binEdges2[5] = AliHFEtools::MakeLinearBinning(nBin[5], 0, nBin[5]);

  fPIDperformance = new THnSparseF("PIDperformance", "PID performance; pT [GeV/c]; theta [rad]; phi [rad]; charge; type (0 - not el, 1 - other el, 2 - HF el; flavor (0 - no, 1 - charm, 2 - bottom)", nDim, nBin);
  fPIDperformance->Sumw2();
  fSignalToBackgroundMC = new THnSparseF("SignalToBackgroundMC", "PID performance; pT [GeV/c]; theta [rad]; phi [rad]; charge; flavor (0 - no, 1 - charm, 2 - bottom); ITS Cluster (0 - no, 1 - first (and maybe second), 2 - second)", nDim, nBin);
  fSignalToBackgroundMC->Sumw2();
  for(Int_t idim = 0; idim < nDim; idim++){
    fPIDperformance->SetBinEdges(idim, binEdges2[idim]);
    fSignalToBackgroundMC->SetBinEdges(idim, binEdges2[idim]); 
  }
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    delete binEdges[ivar];
  for(Int_t ivar = kNvar; ivar < nDim; ivar++)
    delete binEdges2[ivar];
}

//____________________________________________________________
void AliAnalysisTaskHFE::AddPIDdetector(TString detector){
  //
  // Adding PID detector to the task
  //
  if(!fPIDdetectors.Length()) 
    fPIDdetectors = detector;
  else
    fPIDdetectors += ":" + detector;
}

//____________________________________________________________
void AliAnalysisTaskHFE::PrintStatus() const {
  //
  // Print Analysis status
  //
  printf("\n\tAnalysis Settings\n\t========================================\n\n");
  printf("\tSecondary Vertex finding: %s\n", GetPlugin(kSecVtx) ? "YES" : "NO");
  printf("\tPrimary Vertex resolution: %s\n", GetPlugin(kPriVtx) ? "YES" : "NO");
  printf("\n");
  printf("\tParticle Identification Detectors:\n");
  TObjArray *detectors = fPIDdetectors.Tokenize(":");
  for(Int_t idet = 0; idet < detectors->GetEntries(); idet++)
    printf("\t\t%s\n", (dynamic_cast<TObjString *>(detectors->At(idet)))->String().Data());
  printf("\n");
  printf("\tQA: \n");
  printf("\t\tPID: %s\n", IsQAOn(kPIDqa) ? "YES" :  "NO");
  printf("\t\tCUTS: %s\n", (fCuts != NULL && fCuts->IsInDebugMode()) ? "YES" : "NO");
  printf("\t\tMC: %s\n", IsQAOn(kMCqa) ? "YES" : "NO");
  printf("\n");
}

//____________________________________________________________
AliAnalysisTaskHFE::LabelContainer::LabelContainer(Int_t capacity):
  fContainer(NULL),
  fBegin(NULL),
  fEnd(NULL),
  fLast(NULL),
  fCurrent(NULL)
{
  //
  // Default constructor
  //
  fContainer = new Int_t[capacity];
  fBegin = &fContainer[0];
  fEnd = &fContainer[capacity - 1];
  fLast = fCurrent = fBegin;
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFE::LabelContainer::Append(Int_t label){
  //
  // Add Label to the container
  //
  if(fLast > fEnd) return kFALSE;
  *fLast++ = label;
  return kTRUE;
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFE::LabelContainer::Find(Int_t label) const {
  //
  // Find track in the list of labels
  //
  for(Int_t *entry = fBegin; entry <= fLast; entry++) 
    if(*entry == label) return kTRUE;
  return kFALSE;
}

//____________________________________________________________
Int_t AliAnalysisTaskHFE::LabelContainer::Next(){
  //
  // Mimic iterator
  //
  if(fCurrent > fLast) return -1; 
  return *fCurrent++;
}

//____________________________________________________________
Int_t AliAnalysisTaskHFE::IsSignalElectron(AliVParticle *fTrack) const{
  //
  // Checks whether the identified electron track is coming from heavy flavour
  // returns 0 in case of no signal, 1 in case of charm and 2 in case of Bottom
  //
  enum{
    kNoSignal = 0,
    kCharm = 1,
    kBeauty = 2
  };
  TString objname = fTrack->IsA()->GetName();
  Int_t pid = 0;
  if(IsESDanalysis()){
    // ESD Analysis
    AliMCParticle *mctrack = NULL;
    if(!objname.CompareTo("AliESDtrack")){
      AliDebug(2, "Checking signal for ESD track");
      AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(fTrack);
      mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(esdtrack->GetLabel())));
    }
    else if(!objname.CompareTo("AliMCParticle")){
      AliDebug(2, "Checking signal for MC track");
      mctrack = dynamic_cast<AliMCParticle *>(fTrack);
    }
    else{
      AliError("Input object not supported");
      return kNoSignal;
    }
    if(!mctrack) return kNoSignal;
    TParticle *ecand = mctrack->Particle(); 
    if(TMath::Abs(ecand->GetPdgCode()) != 11) return kNoSignal; // electron candidate not true electron
    Int_t motherLabel = TMath::Abs(ecand->GetFirstMother());
    AliDebug(3, Form("mother label: %d\n", motherLabel));
    if(!motherLabel) return kNoSignal; // mother track unknown
    AliMCParticle *motherTrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(motherLabel));
    if(!motherTrack) return kNoSignal;
    TParticle *mparticle = motherTrack->Particle();
    pid = TMath::Abs(mparticle->GetPdgCode());
  } else {
    // AOD Analysis - Different Data handling
    AliAODMCParticle *aodmc = NULL;
    if(!objname.CompareTo("AliAODTrack")){
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(fTrack);
      Int_t aodlabel = TMath::Abs(aodtrack->GetLabel());
      if(aodlabel >= fMCEvent->GetNumberOfTracks()) return kNoSignal;
      aodmc = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(aodlabel));
    } else if(!objname.CompareTo("AliAODMCParticle")){
      aodmc = dynamic_cast<AliAODMCParticle *>(fTrack);
    } else{
      AliError("Input object not supported");
      return kNoSignal;
    }
    if(!aodmc) return kNoSignal;
    Int_t motherLabel = TMath::Abs(aodmc->GetMother());
    AliDebug(3, Form("mother label: %d\n", motherLabel));
    if(!motherLabel || motherLabel >= fMCEvent->GetNumberOfTracks()) return kNoSignal;
    AliAODMCParticle *aodmother = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(motherLabel));
    pid = aodmother->GetPdgCode();
  }
  // From here the two analysis modes go together
  AliDebug(3, Form("PDG code: %d\n", pid));

  // identify signal according to Pdg Code 
  if((pid % 1000) / 100 == 4) return kCharm;    // charmed meson, 3rd position in pdg code == 4
  if(pid / 1000 == 4) return kCharm;            // charmed baryon, 4th position in pdg code == 4
  if((pid % 1000) / 100 == 5) return kBeauty;   // beauty meson, 3rd position in pdg code == 5
  if(pid / 1000 == 5) return kBeauty;           // beauty baryon, 4th position in pdg code == 5   
  return kNoSignal;
}

//__________________________________________
void AliAnalysisTaskHFE::SwitchOnPlugin(Int_t plug){
  //
  // Switch on Plugin
  // Available:
  //  - Primary vertex studies
  //  - Secondary vertex Studies
  //  - Post Processing
  //
  switch(plug){
    case kPriVtx: SETBIT(fPlugins, plug); break;
    case kSecVtx: SETBIT(fPlugins, plug); break;
    case kIsElecBackGround: SETBIT(fPlugins, plug); break;
    case kPostProcess: SETBIT(fPlugins, plug); break;
    default: AliError("Unknown Plugin");
  };
}

//__________________________________________
Bool_t AliAnalysisTaskHFE::ProcessCutStep(Int_t cutStep, AliVParticle *track, Double_t *container, Bool_t signal, Bool_t alreadyseen){
  //
  // Check single track cuts for a given cut step
  // Fill the particle container
  //
  if(!fCFM->CheckParticleCuts(cutStep, track)) return kFALSE;
  if(signal) {
    fCFM->GetParticleContainer()->Fill(container, cutStep + 2*AliHFEcuts::kNcutStepsESDtrack);
    fCFM->GetParticleContainer()->Fill(&container[4], cutStep);
    if(alreadyseen) {
      fCFM->GetParticleContainer()->Fill(&container[4], cutStep + AliHFEcuts::kNcutStepsESDtrack);
    }
  }
  return kTRUE;
}

//__________________________________________
void AliAnalysisTaskHFE::SetTPCBetheBlochParameters(Double_t *pars){
  //
  // Set Bethe-Bloch Parameters for TPC PID
  //
  fPID->SetTPCBetheBlochParameters(pars);
}

