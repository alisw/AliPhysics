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
#include <TH3D.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TProfile.h>
#include <TString.h>
#include <TF1.h>
#include <TTree.h>

#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTriggerAnalysis.h"
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

ClassImp(AliAnalysisTaskHFE)

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE():
  AliAnalysisTaskSE("PID efficiency Analysis")
  , fQAlevel(0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fPlugins(0)
  , fWeighting(kFALSE)
  , fWeightFactors(NULL)
  , fWeightFactorsFunction(NULL)
  , fBackGroundFactorsFunction(NULL)
  , fCFM(NULL)
  , fV0CF(NULL)
  , fTriggerAnalysis(NULL)
  , fHadronicBackground(NULL)
  , fCorrelation(NULL)
  , fPIDperformance(NULL)
  , fSignalToBackgroundMC(NULL)
  , fPID(NULL)
  , fPIDtagged(NULL)
  , fPIDpreselect(NULL)
  , fCuts(NULL)
  , fCutsTagged(NULL)
  , fCutspreselect(NULL)
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
  , fQACollection(NULL)
{
  //
  // Dummy constructor
  //
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const char * name):
  AliAnalysisTaskSE(name)
  , fQAlevel(0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fPlugins(0)
  , fWeighting(kFALSE)
  , fWeightFactors(NULL)
  , fWeightFactorsFunction(NULL)
  , fBackGroundFactorsFunction(NULL)
  , fCFM(NULL)
  , fV0CF(NULL)
  , fTriggerAnalysis(NULL)
  , fHadronicBackground(NULL)
  , fCorrelation(NULL)
  , fPIDperformance(NULL)
  , fSignalToBackgroundMC(NULL)
  , fPID(NULL)
  , fPIDtagged(NULL)
  , fPIDpreselect(NULL)
  , fCuts(NULL)
  , fCutsTagged(NULL)
  , fCutspreselect(NULL)
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
  , fQACollection(0x0)
{
  //
  // Default constructor
  // 
  DefineOutput(1, TH1I::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());

  // Initialize cuts
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref):
  AliAnalysisTaskSE(ref)
  , fQAlevel(0)
  , fPIDdetectors()
  , fPIDstrategy(0)
  , fPlugins(0)
  , fWeighting(kFALSE)
  , fWeightFactors(NULL)
  , fWeightFactorsFunction(NULL)
  , fBackGroundFactorsFunction(NULL)
  , fCFM(NULL)
  , fV0CF(NULL)
  , fTriggerAnalysis(NULL)
  , fHadronicBackground(NULL)
  , fCorrelation(NULL)
  , fPIDperformance(NULL)
  , fSignalToBackgroundMC(NULL)
  , fPID(NULL)
  , fPIDtagged(NULL)
  , fPIDpreselect(NULL)
  , fCuts(NULL)
  , fCutsTagged(NULL)
  , fCutspreselect(NULL)
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
  , fQACollection(NULL)
{
  //
  // Copy Constructor
  //
  ref.Copy(*this);
}

//____________________________________________________________
AliAnalysisTaskHFE &AliAnalysisTaskHFE::operator=(const AliAnalysisTaskHFE &ref){
  //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}

//____________________________________________________________
void AliAnalysisTaskHFE::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliAnalysisTaskHFE &target = dynamic_cast<AliAnalysisTaskHFE &>(o);
  target.fQAlevel = fQAlevel;
  target.fPIDdetectors = fPIDdetectors;
  target.fPIDstrategy = fPIDstrategy;
  target.fPlugins = fPlugins;
  target.fWeighting = fWeighting;
  target.fWeightFactors = fWeightFactors;
  target.fWeightFactorsFunction = fWeightFactorsFunction;
  target.fBackGroundFactorsFunction = fBackGroundFactorsFunction;
  target.fCFM = fCFM;
  target.fV0CF = fV0CF;
  target.fTriggerAnalysis = fTriggerAnalysis;
  target.fHadronicBackground = fHadronicBackground;
  target.fCorrelation = fCorrelation;
  target.fPIDperformance = fPIDperformance;
  target.fSignalToBackgroundMC = fSignalToBackgroundMC;
  target.fPID = fPID;
  target.fPIDtagged = fPIDtagged;
  target.fPIDpreselect = fPIDpreselect;
  target.fCuts = fCuts;
  target.fCutspreselect = fCutspreselect;
  target.fSecVtx = fSecVtx;
  target.fElecBackGround = fElecBackGround;
  target.fMCQA = fMCQA;
  target.fNEvents = fNEvents;
  target.fNElectronTracksEvent = fNElectronTracksEvent;
  target.fQA = fQA;
  target.fOutput = fOutput;
  target.fHistMCQA = fHistMCQA;
  target.fHistSECVTX = fHistSECVTX;
  target.fHistELECBACKGROUND = fHistELECBACKGROUND;
  target.fQACollection = fQACollection;
}

//____________________________________________________________
AliAnalysisTaskHFE::~AliAnalysisTaskHFE(){
  //
  // Destructor
  //
  return;
  if(fPID) delete fPID;
  if(fPIDtagged) delete fPIDtagged;
  if(fQA){
    fQA->Clear();
    delete fQA;
  }
  if(fOutput){ 
    fOutput->Clear();
    delete fOutput;
  }
  if(fWeightFactors) delete fWeightFactors;
  if(fWeightFactorsFunction) delete fWeightFactorsFunction;
  if(fBackGroundFactorsFunction) delete fBackGroundFactorsFunction;
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
  fPID = new AliHFEpid("standardPID"); fPIDtagged = new AliHFEpid("taggedPID");
  AliDebug(3, "Creating Output Objects");
  // Automatic determination of the analysis mode
  AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
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

  // Enable Trigger Analysis
  fTriggerAnalysis = new AliTriggerAnalysis;
  fTriggerAnalysis->EnableHistograms();
  fTriggerAnalysis->SetAnalyzeMC(HasMCData());

  // Make QA histograms
  fNEvents = new TH1I("nEvents", "Number of Events in the Analysis", 2, 0, 2); // Number of Events neccessary for the analysis and not a QA histogram
  fNElectronTracksEvent = new TH1I("nElectronTracksEvent", "Number of Electron Candidates", 100, 0, 100);
  // First Step: TRD alone
  if(!fQA) fQA = new TList;

  fQACollection = new AliHFEcollection("TaskQA", "QA histos from the Electron Task");
  fQACollection->CreateProfile("conr", "Electron PID contamination", 20, 0, 20);
  fQACollection->CreateTH1F("alpha_rec", "Alpha from reconstructed tracks with TRD hits", 36, -TMath::Pi(), TMath::Pi());
  fQACollection->CreateTH1F("alpha_sim", "Alpha from simulated electron tracks", 36, -TMath::Pi(), TMath::Pi());
  fQACollection->CreateTH1F("nElectron", "Number of electrons", 100, 0, 100);
  fQACollection->CreateProfile("pidquality", "TRD PID quality as function of momentum", 20, 0, 20);
  fQACollection->CreateProfile("ntrdclusters", "Number of TRD clusters as function of momentum", 20, 0, 20);
  fQACollection->CreateTH1F("chi2TRD","#chi2 per TRD cluster", 20, 0, 20);
  fQACollection->CreateTH1F("mccharge", "MC Charge", 200, -100, 100);
  fQACollection->CreateTH2F("radius", "Production Vertex", 100, 0.0, 5.0, 100, 0.0, 5.0);
  fQACollection->CreateTH1F("secvtxept", "pT of tagged e", 500, 0, 50); // mj: will move to another place soon
  fQACollection->CreateTH2F("secvtxeTPCsig", "TPC signal for tagged e",125, 0, 25, 200, 0, 200 ); // mj: will move to another place soon 
  fQA->Add(fQACollection->GetList());

  if(!fOutput) fOutput = new TList;
  // Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  fV0CF = new AliCFManager;
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
  // Make clone for V0 tagging step
  fCutsTagged = new AliHFEcuts(*fCuts);
  fCutsTagged->SetName("hfeV0Cuts");
  fCutsTagged->SetTitle("Cuts for tagged Particles");
  fCuts->Initialize(fCFM);
  fCutsTagged->Initialize(fV0CF);
  if(fCuts->IsQAOn()) fQA->Add(fCuts->GetQAhistograms());
  if(fCutsTagged->IsQAOn()) fQA->Add(fCutsTagged->GetQAhistograms());
 
  // add output objects to the List
  fOutput->AddAt(fCFM->GetParticleContainer(), 0);
  fOutput->AddAt(fCFM->GetEventContainer(), 1);
  fOutput->AddAt(fCorrelation, 2);
  fOutput->AddAt(fPIDperformance, 3);
  fOutput->AddAt(fSignalToBackgroundMC, 4);
  fOutput->AddAt(fNElectronTracksEvent, 5);
  fOutput->AddAt(fHadronicBackground, 6);
  fOutput->AddAt(fV0CF, 7);
  // Initialize PID
  if(IsQAOn(kPIDqa)){
    AliInfo("PID QA switched on");
    //fPID->SetDebugLevel(2);
    fPID->SetQAOn();
    fQA->Add(fPID->GetQAhistograms());
  }
  fPID->SetHasMCData(HasMCData());
  if(!fPIDdetectors.Length() && ! fPIDstrategy) AddPIDdetector("TPC");
  if(fPIDstrategy){
    fPID->InitializePID(Form("Strategy%d", fPIDstrategy));
    fPIDtagged->InitializePID(Form("Strategy%d", fPIDstrategy));
  }
  else{
    fPID->InitializePID(fPIDdetectors.Data());     // Only restrictions to TPC allowed 
    fPIDtagged->InitializePID(fPIDdetectors.Data());     // Only restrictions to TPC allowed 
  }

  // mcQA----------------------------------
  if (HasMCData() && IsQAOn(kMCqa)) {
    AliInfo("MC QA on");
    if(!fMCQA) fMCQA = new AliHFEmcQA;
    if(!fHistMCQA) fHistMCQA = new TList();
    fMCQA->CreatDefaultHistograms(fHistMCQA);
    fQA->Add(fHistMCQA);
  } 

  // secvtx----------------------------------
  if (GetPlugin(kSecVtx)) {
    AliInfo("Secondary Vertex Analysis on");
    if(!fSecVtx) fSecVtx = new AliHFEsecVtx;
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

  if(IsESDanalysis() && HasMCData()){
    // Protect against missing MC trees
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH){ 
      AliError("No MC Event Handler available");
      return;
    }
    if(!mcH->InitOk()) return;
    if(!mcH->TreeK()) return;
    if(!mcH->TreeTR()) return;
  }

  // Protect agains missing 
  if(HasMCData()) ProcessMC();  // Run the MC loop + MC QA in case MC Data are available

  if(IsAODanalysis()) ProcessAOD();
  else{
    AliESDInputHandler *inH = dynamic_cast<AliESDInputHandler *>(fInputHandler);
    if(!inH){
      AliError("No ESD Input handler available");
      return;
    }
    AliESDpid *workingPID = inH->GetESDpid();
    if(workingPID){
      AliDebug(1, "Using ESD PID from the input handler");
      fPID->SetESDpid(workingPID);
      fPIDtagged->SetESDpid(workingPID);
      if(fPIDpreselect) fPIDpreselect->SetESDpid(workingPID);
    } else { 
      AliDebug(1, "Using default ESD PID");
      fPID->SetESDpid(AliHFEtools::GetDefaultPID(HasMCData()));
      fPIDtagged->SetESDpid(AliHFEtools::GetDefaultPID(HasMCData()));
      if(fPIDpreselect) fPIDpreselect->SetESDpid(AliHFEtools::GetDefaultPID(HasMCData())); 
    }
    ProcessESD();
  }
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
//_______________________________________________________________
Bool_t AliAnalysisTaskHFE::IsEventInBinZero() {
  //
  //
  //

  //printf("test in IsEventInBinZero\n");
  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return kFALSE;
  }

  // check vertex
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  if(!vertex) return kTRUE;
  //if(vertex) return kTRUE;

  // check tracks
  if(fInputEvent->GetNumberOfTracks()<=0) return kTRUE;
  //if(fInputEvent->GetNumberOfTracks()>0) return kTRUE;
  
  
  return kFALSE;
  
}
//____________________________________________________________
void AliAnalysisTaskHFE::ProcessMC(){
  //
  // Runs the MC Loop (filling the container for the MC Cut Steps with the observables pt, eta and phi)
  // In case MC QA is on also MC QA loop is done
  //
  AliDebug(3, "Processing MC Information");
  Double_t eventContainer [2];
  eventContainer[0] = fMCEvent->GetPrimaryVertex()->GetZ();
  if(fCFM->CheckEventCuts(AliHFEcuts::kEventStepGenerated, fMCEvent)) 
    fCFM->GetEventContainer()->Fill(eventContainer,AliHFEcuts::kEventStepGenerated);
  Int_t nElectrons = 0;
  if(IsESDanalysis()){
    if (HasMCData() && IsQAOn(kMCqa)) {
      AliDebug(2, "Running MC QA");

      if(fMCEvent->Stack()){
	      fMCQA->SetMCEvent(fMCEvent);
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
          fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kOthers, AliHFEmcQA::kElectronPDG, 0); // no accept cut
          if (TMath::Abs(mcpart->Eta()) < 0.9) {
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 1); // accept |eta|<0.9
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 1); // accept |eta|<0.9
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kOthers, AliHFEmcQA::kElectronPDG, 1); // accept |eta|<0.9
          }
          if (TMath::Abs(AliHFEtools::GetRapidity(mcpart)) < 0.5) {
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 2); // accept |y|<0.5
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 2); // accept |y|<0.5
            fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kOthers, AliHFEmcQA::kElectronPDG, 2); // accept |y|<0.5
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
  AliDebug(3, Form("Number of Tracks: %d", fMCEvent->GetNumberOfTracks()));
  for(Int_t imc = 0; imc <fMCEvent->GetNumberOfTracks(); imc++){
    if(!(mctrack = fMCEvent->GetTrack(imc))) continue;
    AliDebug(4, "Next Track");
    if(ProcessMCtrack(mctrack)) nElectrons++;
  }

  // fCFM->CheckEventCuts(AliCFManager::kEvtRecCuts, fESD);
  fQACollection->Fill("nElectron", nElectrons);
}

//____________________________________________________________
void AliAnalysisTaskHFE::ProcessESD(){
  //
  // Run Analysis of reconstructed event in ESD Mode
  // Loop over Tracks, filter according cut steps defined in AliHFEcuts
  //
  AliDebug(3, "Processing ESD Event");
  AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!fESD){
    AliError("ESD Event required for ESD Analysis")
    return;
  }

  // Do event Normalization
  Double_t eventContainer[2];
  eventContainer[0] = fInputEvent->GetPrimaryVertex()->GetZ();
  eventContainer[1] = 0.;
  if(fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0AND))
    eventContainer[1] = 1.;
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoCut);
  if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fInputEvent)) return;
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepReconstructed);

  if (GetPlugin(kIsElecBackGround)) { 
    fElecBackGround->SetEvent(fESD);
  }
  if (GetPlugin(kSecVtx)) {
    fSecVtx->SetEvent(fESD);
    fSecVtx->GetPrimaryCondition();
  }

  if(HasMCData()){
    if (GetPlugin(kSecVtx)) { 
      fSecVtx->SetMCEvent(fMCEvent);
    }
    if (GetPlugin(kIsElecBackGround)) { 
      fElecBackGround->SetMCEvent(fMCEvent);
    }
  }

  Double_t nContrib = fInputEvent->GetPrimaryVertex()->GetNContributors();
  Double_t container[10];
  memset(container, 0, sizeof(Double_t) * 10);
  // container for the output THnSparse
  Double_t dataE[6]; // [pT, eta, Phi, type, 'C' or 'B']
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
  AliDebug(3, Form("Number of Tracks: %d", fESD->GetNumberOfTracks()));
  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    track = fESD->GetTrack(itrack);

    // fill counts of v0-identified particles
    Int_t v0pid = -1;
    if(track->TestBit(BIT(14))) v0pid = AliPID::kElectron;
    else if(track->TestBit(BIT(15))) v0pid = AliPID::kPion;
    else if(track->TestBit(BIT(16))) v0pid = AliPID::kProton;
    if(v0pid > -1)
      FilterTaggedTrack(track, v0pid);
 
    AliDebug(3, Form("Doing track %d, %p", itrack, track));
     
    //////////////////////////////////////
    // preselect
    /////////////////////////////////////
    if(fPIDpreselect && fCutspreselect) {
      if(!PreSelectTrack(track)) continue;
    }
     
    container[0] = track->Pt();
    container[1] = track->Eta();
    container[2] = track->Phi();
    container[3] = track->Charge();
    container[4] = 0;

    dataE[0] = track->Pt();
    dataE[1] = track->Eta();
    dataE[2] = track->Phi();
    dataE[3] = track->Charge();
    dataE[4] = -1;
    dataE[5] = -1;

    signal = kTRUE;
    Double_t weight = 1.0;
    
    // Fill step without any cut
          
    if(HasMCData()){
      container[4] = container[9] = kOther;
      // Check if it is electrons near the vertex
      if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(track->GetLabel()))))) continue;
      mctrack4QA = mctrack->Particle();

      container[5] = mctrack->Pt();
      container[6] = mctrack->Eta();
      container[7] = mctrack->Phi();
      container[8] = mctrack->Charge()/3.;

      if(fWeighting) weight = FindWeight(container[5],container[6],container[7]);    

      if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE;
      else AliDebug(3, "Signal Electron");
      
      Int_t signalTrack = 0;
      if((signalTrack = IsSignalElectron(track))){
        AliDebug(3, Form("Signal: Index = %d\n", signalTrack));
        switch(signalTrack){
          case 1: container[4] = container[9] = kSignalCharm; break;
          case 2: container[4] = container[9] = kSignalBeauty; break;
          default: container[4] = container[9] = kOther; break;
        };
      } else if(IsGammaElectron(track)) container[4] = container[9] = kGammaConv;
      AliDebug(3, Form("Signal Decision(%f/%f)", container[4], container[9]));
    } 
    AliDebug(3, Form("Weight? %f", weight));
    if(signal) {
      alreadyseen = cont.Find(TMath::Abs(track->GetLabel()));
      cont.Append(TMath::Abs(track->GetLabel()));
      
      fCFM->GetParticleContainer()->Fill(&container[5], AliHFEcuts::kStepRecNoCut,weight);
      fCFM->GetParticleContainer()->Fill(&container[0], AliHFEcuts::kStepRecNoCut + 2*AliHFEcuts::kNcutStepsESDtrack,weight);
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[5], AliHFEcuts::kStepRecNoCut + AliHFEcuts::kNcutStepsESDtrack,weight);
      }
    }

    // RecKine: ITSTPC cuts  
    if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track, container, signal, alreadyseen, weight)) continue;
    
    // Check TRD criterions (outside the correction framework)
    if(track->GetTRDncls()){
      fQACollection->Fill("chi2TRD", track->GetTRDchi2()/track->GetTRDncls());
      fQACollection->Fill("alpha_rec", track->GetAlpha());
      fQACollection->Fill("pidquality", container[0], track->GetTRDpidQuality());
      fQACollection->Fill("ntrdclusters", container[0], track->GetTRDncls());
    }

    
    // RecPrim
    if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track, container, signal, alreadyseen,weight)) continue;

    // HFEcuts: ITS layers cuts
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track, container, signal, alreadyseen,weight)) continue;

    // HFEcuts: Nb of tracklets TRD0
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTRD, track, container, signal, alreadyseen,weight)) continue;
    if(signal) {
      // dimensions 3&4&5 : pt,eta,phi (MC)
      ((THnSparseF *)fCorrelation->At(0))->Fill(container);
    }

    if(HasMCData() && IsQAOn(kMCqa)) {
      // mc qa for after the reconstruction cuts  
      AliDebug(2, "Running MC QA");
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 3);  // charm
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 3); // beauty 
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kOthers,  AliHFEmcQA::kElectronPDG, 3); // beauty 
    } 

    if(HasMCData()){
      FillProductionVertex(track);
    }

    // track accepted, do PID
    AliHFEpidObject hfetrack;
    hfetrack.fAnalysisType = AliHFEpidObject::kESDanalysis;
    hfetrack.fRecTrack = track;
    if(HasMCData()) hfetrack.fMCtrack = mctrack;
    if(!fPID->IsSelected(&hfetrack)) continue;
    nElectronCandidates++;

    // Fill Histogram for Hadronic Background
    if(HasMCData()){
      if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) != 11))
        fHadronicBackground->Fill(container, 0);
    }

    if (HasMCData() && IsQAOn(kMCqa)) {
      // mc qa for after the reconstruction and pid cuts  
      AliDebug(2, "Running MC QA");
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 4);  // charm
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 4); // beauty 
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kOthers,  AliHFEmcQA::kElectronPDG, 4); // beauty 
    } 

    // Fill Containers
    if(signal) {
      // Apply weight for background contamination
      if(fBackGroundFactorsFunction) {
	      Double_t weightBackGround = fBackGroundFactorsFunction->Eval(TMath::Abs(track->P()));
	      if(weightBackGround < 0.0) weightBackGround = 0.0;
        else if(weightBackGround > 1.0) weightBackGround = 0.0;
        fHadronicBackground->Fill(container, 1, weight * weightBackGround);
      }
      //      
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepPID + 2*AliHFEcuts::kNcutStepsESDtrack, weight);
      fCFM->GetParticleContainer()->Fill(&container[5], AliHFEcuts::kStepPID, weight);
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[5], (AliHFEcuts::kStepPID + (AliHFEcuts::kNcutStepsESDtrack)),weight);
      }
      // dimensions 3&4&5 : pt,eta,phi (MC)
      ((THnSparseF *)fCorrelation->At(1))->Fill(container);
    }

    if(GetPlugin(kSecVtx)) {
      AliDebug(2, "Running Secondary Vertex Analysis");
      if(track->Pt()>2.0 && nContrib > 1){ 
        fSecVtx->InitHFEpairs();
        fSecVtx->InitHFEsecvtxs();
        for(Int_t jtrack = 0; jtrack < fESD->GetNumberOfTracks(); jtrack++){
          htrack = fESD->GetTrack(jtrack);
          if ( itrack == jtrack ) continue; // since it is for tagging single electron, don't need additional condition 
          if (htrack->Pt()<2.0) continue;
          if (!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC, htrack)) continue;
          if (!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim, htrack)) continue;
          fSecVtx->PairAnalysis(track, htrack, jtrack); // e-h pairing
        }
        for(int ip=0; ip<fSecVtx->HFEpairs()->GetEntriesFast(); ip++){
          //if(HasMCData()){
            AliHFEpairs *pair = (AliHFEpairs*) (fSecVtx->HFEpairs()->UncheckedAt(ip));
            //if(!(pair->GetPairCode()>1. && pair->GetPairCode()<4.))  // apply various cuts
            // apply various cuts
            if(pair->GetKFChi2()>5.) // only apply vertex chi2 cut for the moment
            //if((pair->GetKFChi2()>5.) || !(pair->GetSignedLxy()>0. && pair->GetSignedLxy()<2.)) 
              fSecVtx->HFEpairs()->RemoveAt(ip);
          //}
        }
        fSecVtx->HFEpairs()->Compress();
        if(fSecVtx->HFEpairs()->GetEntriesFast()) fSecVtx->RunSECVTX(track); // secondary vertexing with e,h1,h2,.. tracks
        for(int ip=0; ip<fSecVtx->HFEsecvtxs()->GetEntriesFast(); ip++){
          AliHFEsecVtxs *secvtx=0x0;
          secvtx = (AliHFEsecVtxs*) (fSecVtx->HFEsecvtxs()->UncheckedAt(ip));
          if(!(secvtx->GetInvmass()>2.0 && secvtx->GetInvmass()<5.2) || !(secvtx->GetSignedLxy2()>0.08 && secvtx->GetSignedLxy2()<1.5) || !(secvtx->GetKFIP2()>-0.1 && secvtx->GetKFIP2()<0.1))
            fSecVtx->HFEsecvtxs()->RemoveAt(ip);
          // here you apply cuts, then if it doesn't pass the cut, remove it from the fSecVtx->HFEsecvtxs() 
        }
        if(fSecVtx->HFEsecvtxs()->GetEntriesFast()) {
          fQACollection->Fill("secvtxpt", track->Pt());
          fQACollection->Fill("secvtxTPCsig", track->P(),track->GetTPCsignal());
          if (HasMCData() && IsQAOn(kMCqa)) {
            // mc qa for after the reconstruction and pid cuts  
            AliDebug(2, "Running MC QA");
            fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 5);  // charm
            fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 5); // beauty 
            fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kOthers,  AliHFEmcQA::kElectronPDG, 5); // beauty 
          }
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
  Double_t eventContainer[2];
  eventContainer[0] = fInputEvent->GetPrimaryVertex()->GetZ();
  eventContainer[1] = 1.; // No Information available in AOD analysis, assume all events have V0AND
  fCFM->GetEventContainer()->Fill(eventContainer,AliHFEcuts::kEventStepRecNoCut);
  if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fInputEvent)) return;
  fCFM->GetEventContainer()->Fill(eventContainer,AliHFEcuts::kEventStepReconstructed);
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
  if(!fAOD){
    AliError("AOD Event required for AOD Analysis")
    return;
  }
 
  AliAODTrack *track = NULL;
  AliAODMCParticle *mctrack = NULL;
  Double_t container[10]; memset(container, 0, sizeof(Double_t) * 10);
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
      Int_t signalTrack = 0;
      if((signalTrack = IsSignalElectron(track))){
        switch(signalTrack){
          case 1: container[4] = container[9] = kSignalCharm; break;
          case 2: container[4] = container[9] = kSignalBeauty; break;
        };
      } else if(IsGammaElectron(track)) 
        container[4] = container[9] = kGammaConv;
      else container[4] = container[9] = kOther;

      Int_t label = TMath::Abs(track->GetLabel());
      if(label){
        mctrack = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(label));
        container[5] = mctrack->Pt();
        container[6] = mctrack->Eta();
        container[7] = mctrack->Phi();
        container[8] = mctrack->Charge();
      }
    }
    // track accepted, do PID
    AliHFEpidObject hfetrack;
    hfetrack.fAnalysisType = AliHFEpidObject::kAODanalysis;
    hfetrack.fRecTrack = track;
    if(HasMCData()) hfetrack.fMCtrack = mctrack;
    //if(!fPID->IsSelected(&hfetrack)) continue;    // we will do PID here as soon as possible
    // Particle identified - Fill CF Container
    // Apply weight for background contamination
    Double_t weightBackGround = 1.0;
    if(fBackGroundFactorsFunction) {
      weightBackGround = weightBackGround - fBackGroundFactorsFunction->Eval(TMath::Abs(track->P()));
      if(weightBackGround < 0.0) weightBackGround = 1.0;
    }
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepPID + 2*AliHFEcuts::kNcutStepsESDtrack, weightBackGround);
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
  Double_t container[5], signalContainer[6];
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
  Int_t signal = 0;
  if((signal = IsSignalElectron(track))){
    switch(signal){
      case 1: container[4] = kSignalCharm; break;
      case 2: container[4] = kSignalBeauty; break;
    };
  }else if(IsGammaElectron(track)) container[4] = kGammaConv;
  else container[4] = kOther;

  // weight
  Double_t weight = 1.0;
  if(fWeighting) weight = FindWeight(container[0],container[1],container[2]);

 if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, track)) return kFALSE;
 fQACollection->Fill("mccharge", signalContainer[3]);
 fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCGenerated,weight);
  if((signalContainer[4] = static_cast<Double_t >(IsSignalElectron(track))) > 1e-3) fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCsignal,weight);
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
  fQACollection->Fill("alpha_sim", container[2] - TMath::Pi());
  //if(IsESDanalysis()){
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCInAcceptance, track)) return kFALSE;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCInAcceptance,weight);
  //}
  return kTRUE;
}

//____________________________________________________________
void AliAnalysisTaskHFE::FilterTaggedTrack(AliESDtrack *track, Int_t species){
  //
  // Filter tracks tagged by V0 PID class
  //
  Int_t offset = AliHFEcuts::kStepRecKineITSTPC;
  Double_t container[5] ={track->Pt(), track->Eta(), track->Phi(), track->Charge(), species};
  fV0CF->GetParticleContainer()->Fill(container, 0); // Fill Container without filtering
  Bool_t survived = kTRUE;
  for(Int_t icut = AliHFEcuts::kStepRecKineITSTPC; icut < AliHFEcuts::kStepPID; icut++){
    AliDebug(2, Form("Checking cut %d for species %s", icut, AliPID::ParticleName(species)));
    if(!fV0CF->CheckParticleCuts(icut, track)){
      survived = kFALSE;
      break;
    }
    AliDebug(2, Form("Cut passed, filling container %d", icut - offset + 1));
    fV0CF->GetParticleContainer()->Fill(container, icut - offset + 1);
  }
  if(survived){
    // Apply PID
    AliHFEpidObject hfetrack;
    hfetrack.fAnalysisType = AliHFEpidObject::kESDanalysis;
    hfetrack.fRecTrack = track;
    if(fPIDtagged->IsSelected(&hfetrack)) fV0CF->GetParticleContainer()->Fill(container, AliHFEcuts::kStepPID - offset + 1); 
  } 
}
//____________________________________________________________
Bool_t AliAnalysisTaskHFE::PreSelectTrack(AliESDtrack *track) const {
  //
  // Preselect tracks
  //
  

  Bool_t survived = kTRUE;
  
  if(!fCutspreselect->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC, track)) {
    survived = kFALSE;
    //printf("Did not pass AliHFEcuts::kStepRecKineITSTPC\n");
  }
  //else printf("Pass AliHFEcuts::kStepRecKineITSTPC\n");
  if(!fCutspreselect->CheckParticleCuts(AliHFEcuts::kStepRecPrim, track)) {
    survived = kFALSE;
    //printf("Did not pass AliHFEcuts::kStepRecPrim\n");
  }
  //else printf("Pass AliHFEcuts::kStepRecPrim\n");
  if(!fCutspreselect->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS, track)) {
    survived = kFALSE;
    //printf("Did not pass AliHFEcuts::kStepHFEcutsITS\n");
  }
  //else printf("Pass AliHFEcuts::kStepHFEcutsITS\n");
  if(!fCutspreselect->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD, track)) {
    survived = kFALSE;
    //printf("Did not pass AliHFEcuts::kStepHFEcutsTRD\n");
  }
  //else printf("Pass AliHFEcuts::kStepHFEcutsTRD\n");
  
  if(survived){
    // Apply PID
    AliHFEpidObject hfetrack;
    hfetrack.fAnalysisType = AliHFEpidObject::kESDanalysis;
    hfetrack.fRecTrack = track;
    if(!fPIDpreselect->IsSelected(&hfetrack)) {
      //printf("Did not pass AliHFEcuts::kPID\n");
      survived = kFALSE;
    }
    //else printf("Pass AliHFEcuts::kPID\n");
  }

  return survived; 
      
}
//____________________________________________________________
void AliAnalysisTaskHFE::MakeEventContainer(){
  //
  // Create the event container for the correction framework and link it
  //
  const Int_t kNvar = 2;  // number of variables on the grid: 
  Int_t nBins[kNvar] = {120, 2};
  Double_t binMin[kNvar] = {-30. , 0.};
  Double_t binMax[kNvar] = {30., 2.};

  AliCFContainer *evCont = new AliCFContainer("eventContainer", "Container for events", AliHFEcuts::kNcutStepsEvent, kNvar, nBins);

  Double_t *vertexBins = AliHFEtools::MakeLinearBinning(nBins[0], binMin[0], binMax[0]);
  Double_t *v0andBins = AliHFEtools::MakeLinearBinning(nBins[1], binMin[1], binMax[1]);
  evCont->SetBinLimits(0, vertexBins);
  evCont->SetBinLimits(1, v0andBins);
  delete[] vertexBins; delete[] v0andBins;

  fCFM->SetEventContainer(evCont);
}

//____________________________________________________________
void AliAnalysisTaskHFE::MakeParticleContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
  //
  const Int_t kNvar   = 5;
  //number of variables on the grid:pt,eta, phi, charge
  const Double_t kPtbound[2] = {0.1, 20.};
  const Double_t kEtabound[2] = {-0.8, 0.8};
  const Double_t kPhibound[2] = {0., 2. * TMath::Pi()};

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 44; // bins in pt
  iBin[1] =  8; // bins in eta 
  iBin[2] = 18; // bins in phi
  iBin[3] =  2; // bins in charge
  iBin[4] =  4; // creation process of the electron

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  binEdges[0] = AliHFEtools::MakeLogarithmicBinning(iBin[0], kPtbound[0], kPtbound[1]);
  binEdges[1] = AliHFEtools::MakeLinearBinning(iBin[1], kEtabound[0], kEtabound[1]);
  binEdges[2] = AliHFEtools::MakeLinearBinning(iBin[2], kPhibound[0], kPhibound[1]);
  binEdges[3] = AliHFEtools::MakeLinearBinning(iBin[3], -1.1, 1.1); // Numeric precision
  binEdges[4] = AliHFEtools::MakeLinearBinning(iBin[4], 0, iBin[4]); // Numeric precision
  //for(Int_t ib = 0; ib <= iBin[4]; ib++) printf("%f\t", binEdges[4][ib]);
  //printf("\n");

  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("trackContainer", "Container for tracks", (AliHFEcuts::kNcutStepsTrack + 2*AliHFEcuts::kNcutStepsESDtrack), kNvar, iBin);
  fHadronicBackground = new AliCFContainer("hadronicBackground", "Container for hadronic Background", 2, kNvar, iBin);

  //setting the bin limits
  for(Int_t ivar = 0; ivar < kNvar; ivar++){
    container -> SetBinLimits(ivar, binEdges[ivar]);
    fHadronicBackground -> SetBinLimits(ivar, binEdges[ivar]);
  }
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

  // create correction framework container for V0-tagged particles, new bin limits in 4th bin
  iBin[4] = 5;
  delete binEdges[4];
  binEdges[4] = AliHFEtools::MakeLinearBinning(iBin[4], 0, iBin[4]); // Numeric precision
  AliCFContainer *tagged = new AliCFContainer("taggedTrackContainer", "Correction Framework Container for tagged tracks", AliHFEcuts::kNcutStepsESDtrack, kNvar, iBin);
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    tagged->SetBinLimits(ivar, binEdges[ivar]);
  fV0CF->SetParticleContainer(tagged);

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
  printf("\t\tCUTS: %s\n", (fCuts != NULL && fCuts->IsQAOn()) ? "YES" : "NO");
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
Int_t AliAnalysisTaskHFE::IsSignalElectron(const AliVParticle * const track) const{
  //
  // Checks whether the identified electron track is coming from heavy flavour
  // returns 0 in case of no signal, 1 in case of charm and 2 in case of Bottom
  //
  enum{
    kNoSignal = 0,
    kCharm = 1,
    kBeauty = 2
  };

  if(!fMCEvent) return kNoSignal;
  const AliVParticle *motherParticle = NULL, *mctrack = NULL;
  TString objectType = track->IsA()->GetName();
  Int_t label = 0;
  if(objectType.CompareTo("AliESDtrack") == 0 || objectType.CompareTo("AliAODTrack") == 0){
    // Reconstructed track
    if((label = TMath::Abs(track->GetLabel())) && label < fMCEvent->GetNumberOfTracks())
      mctrack = fMCEvent->GetTrack(label);
  } else {
    // MCParticle
    mctrack = track;
  }

  if(!mctrack) return kNoSignal;
  
  Int_t pid = 0;
  Int_t daughterPDG = 0, motherLabel = 0;
  if(TString(mctrack->IsA()->GetName()).CompareTo("AliMCParticle") == 0){
    // case MC Particle
    daughterPDG = TMath::Abs((dynamic_cast<const AliMCParticle *>(mctrack))->Particle()->GetPdgCode());
    motherLabel = (dynamic_cast<const AliMCParticle *>(mctrack))->Particle()->GetFirstMother();
    if(motherLabel >= 0 && motherLabel < fMCEvent->GetNumberOfTracks())
      motherParticle = fMCEvent->GetTrack(motherLabel);
    if(motherParticle)
      pid = TMath::Abs((dynamic_cast<const AliMCParticle *>(motherParticle))->Particle()->GetPdgCode());
  } else {
    // case AODMCParticle
    daughterPDG = TMath::Abs((dynamic_cast<const AliAODMCParticle *>(mctrack))->GetPdgCode());
    motherLabel = (dynamic_cast<const AliAODMCParticle *>(mctrack))->GetMother();
    if(motherLabel >= 0 && motherLabel < fMCEvent->GetNumberOfTracks())
      motherParticle = fMCEvent->GetTrack(motherLabel);
    if(motherParticle)
      pid = TMath::Abs((dynamic_cast<const AliAODMCParticle *>(motherParticle))->GetPdgCode());
  }
  AliDebug(5, Form("Daughter PDG code: %d", daughterPDG));

  if(!pid) return kNoSignal;

  // From here the two analysis modes go together
  AliDebug(5, Form("Mother PDG code: %d", pid));

  // identify signal according to Pdg Code - barions higher ranked than mesons 
  if(pid / 1000 == 4) return kCharm;            // charmed baryon, 4th position in pdg code == 4
  if(pid / 1000 == 5) return kBeauty;           // beauty baryon, 4th position in pdg code == 5   
  if((pid % 1000) / 100 == 4) return kCharm;    // charmed meson, 3rd position in pdg code == 4
  if((pid % 1000) / 100 == 5) return kBeauty;   // beauty meson, 3rd position in pdg code == 5
  return kNoSignal;
}

//__________________________________________
Bool_t AliAnalysisTaskHFE::IsGammaElectron(const AliVParticle * const track) const {
  //
  // Check for MC if the electron is coming from Gamma
  //
  if(!fMCEvent) return kFALSE;
  const AliVParticle *motherParticle = NULL, *mctrack = NULL;
  TString objectType = track->IsA()->GetName();
  if(objectType.CompareTo("AliESDtrack") == 0 || objectType.CompareTo("AliAODTrack") == 0){
    // Reconstructed track
    if(track->GetLabel())
      mctrack = fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
  } else {
    // MCParticle
    mctrack = track;
  }

  if(!mctrack) return kFALSE;
  
  Int_t motherPDG = 0;
  if(TString(mctrack->IsA()->GetName()).CompareTo("AliMCParticle") == 0){
    // case MC Particle
    motherParticle = fMCEvent->GetTrack((dynamic_cast<const AliMCParticle *>(mctrack)->Particle()->GetFirstMother()));
    if(motherParticle)
      motherPDG = TMath::Abs((dynamic_cast<const AliMCParticle *>(motherParticle))->Particle()->GetPdgCode());
  } else {
    // case AODMCParticle
    motherParticle = fMCEvent->GetTrack((dynamic_cast<const AliAODMCParticle *>(mctrack))->GetMother());
    if(motherParticle)
      motherPDG = TMath::Abs((dynamic_cast<const AliAODMCParticle *>(motherParticle))->GetPdgCode());
  }
  if(motherPDG!=22) return kFALSE;
  else return kTRUE;
}
//____________________________________________________________
Bool_t AliAnalysisTaskHFE::FillProductionVertex(const AliVParticle * const track) const{
  //
  // Find the production vertex of the associated MC track
  //
  if(!fMCEvent) return kFALSE;
  const AliVParticle *mctrack = NULL;
  TString objectType = track->IsA()->GetName();
  if(objectType.CompareTo("AliESDtrack") == 0 || objectType.CompareTo("AliAODTrack") == 0){
    // Reconstructed track
    mctrack = fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
  } else {
    // MCParticle
    mctrack = track;
  }

  if(!mctrack) return kFALSE;

  Double_t xv = 0.0;
  Double_t yv = 0.0;
 
  if(TString(mctrack->IsA()->GetName()).CompareTo("AliMCParticle") == 0){
    // case MCParticle
    xv =  (dynamic_cast<const AliMCParticle *>(mctrack)->Xv());
    yv =  (dynamic_cast<const AliMCParticle *>(mctrack)->Yv());
       
  } else {
    // case AODMCParticle
    xv =  (dynamic_cast<const AliAODMCParticle *>(mctrack)->Xv());
    yv =  (dynamic_cast<const AliAODMCParticle *>(mctrack)->Yv());
  }

  //printf("xv %f, yv %f\n",xv,yv);
  fQACollection->Fill("radius", TMath::Abs(xv),TMath::Abs(yv));

  return kTRUE;

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
//_______________________________________________
void AliAnalysisTaskHFE::SetWeightFactors(TH3D * const weightFactors){
  //
  // Set the histos with the weights for the efficiency maps
  //
  fWeighting = kTRUE;
  fWeightFactors = weightFactors;
}
//_______________________________________________
void AliAnalysisTaskHFE::SetWeightFactorsFunction(TF1 * const weightFactorsFunction){
  //
  // Set the histos with the weights for the efficiency maps
  //
  fWeighting = kTRUE;
  fWeightFactorsFunction = weightFactorsFunction;
  //printf("SetWeightFactors\n");
}
//_______________________________________________
Double_t AliAnalysisTaskHFE::FindWeight(Double_t pt, Double_t eta, Double_t phi) const {
  //
  // Find the weight corresponding to pt eta and phi in the TH3D
  //
  Double_t weight = 1.0;
  if(fWeightFactors) {
    
    TAxis *ptaxis = fWeightFactors->GetXaxis();
    TAxis *etaaxis = fWeightFactors->GetYaxis();
    TAxis *phiaxis = fWeightFactors->GetZaxis();
    
    Int_t ptbin = ptaxis->FindBin(pt);
    Int_t etabin = etaaxis->FindBin(eta);
    Int_t phibin = phiaxis->FindBin(phi);


    weight = fWeightFactors->GetBinContent(ptbin,etabin,phibin);
  }
  else if(fWeightFactorsFunction) {
    
    weight = fWeightFactorsFunction->Eval(pt,eta,phi);
    //printf("pt %f and weight %f\n",pt,weight);

  }

  //printf("pt %f, eta %f, phi %f, weight %f\n",pt,eta,phi,weight);
  
  return weight;  

}
//__________________________________________
Bool_t AliAnalysisTaskHFE::ProcessCutStep(Int_t cutStep, AliVParticle *track, Double_t *container, Bool_t signal, Bool_t alreadyseen,Double_t weight){
  //
  // Check single track cuts for a given cut step
  // Fill the particle container
  //
  if(!fCFM->CheckParticleCuts(cutStep, track)) return kFALSE;
  if(signal) {
    fCFM->GetParticleContainer()->Fill(container, cutStep + 2*AliHFEcuts::kNcutStepsESDtrack,weight);
    fCFM->GetParticleContainer()->Fill(&container[5], cutStep,weight);
    if(alreadyseen) {
      fCFM->GetParticleContainer()->Fill(&container[5], cutStep + AliHFEcuts::kNcutStepsESDtrack,weight);
    }
  }
  return kTRUE;
}

