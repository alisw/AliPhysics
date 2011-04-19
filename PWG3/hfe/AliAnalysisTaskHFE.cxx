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
#include <TH3D.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TProfile.h>
//#include <TString.h>
#include <TF1.h>
#include <TTree.h>

#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTriggerAnalysis.h"
#include "AliVVertex.h"

#include "AliHFEcollection.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEelecbackground.h"
#include "AliHFEmcQA.h"
#include "AliHFEpairs.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEpostAnalysis.h"
#include "AliHFEsecVtxs.h"
#include "AliHFEsecVtx.h"
#include "AliHFEsignalCuts.h"
#include "AliHFEtaggedTrackAnalysis.h"
#include "AliHFEtools.h"
#include "AliHFEvarManager.h"
#include "AliAnalysisTaskHFE.h"

ClassImp(AliAnalysisTaskHFE)

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE():
  AliAnalysisTaskSE("PID efficiency Analysis")
  , fQAlevel(0)
  , fPlugins(0)
  , fFillSignalOnly(kTRUE)
  , fBackGroundFactorApply(kFALSE)
  , fRemovePileUp(kFALSE)
  , fIdentifiedAsPileUp(kFALSE)
  , fIdentifiedAsOutInz(kFALSE)
  , fPassTheEventCut(kFALSE)
  , fHasSpecialTriggerSelection(kFALSE)
  , fSpecialTrigger("NONE")
  , fCentralityF(99.0)
  , fContributors(0.5)
  , fWeightBackGround(0.)
  , fVz(0.0)
  , fContainer(NULL)
  , fVarManager(NULL)
  , fSignalCuts(NULL)
  , fCFM(NULL)
  , fTriggerAnalysis(NULL)
  , fPID(NULL)
  , fPIDqa(NULL)
  , fPIDpreselect(NULL)
  , fCuts(NULL)
  , fTaggedTrackCuts(NULL)
  , fCleanTaggedTrack(kFALSE)
  , fVariablesTRDTaggedTrack(kFALSE)
  , fCutspreselect(NULL)
  , fSecVtx(NULL)
  , fElecBackGround(NULL)
  , fMCQA(NULL)
  , fTaggedTrackAnalysis(NULL)
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
  , fPlugins(0)
  , fFillSignalOnly(kTRUE)
  , fBackGroundFactorApply(kFALSE)
  , fRemovePileUp(kFALSE)
  , fIdentifiedAsPileUp(kFALSE)
  , fIdentifiedAsOutInz(kFALSE)
  , fPassTheEventCut(kFALSE)  
  , fHasSpecialTriggerSelection(kFALSE)
  , fSpecialTrigger("NONE")
  , fCentralityF(99.0)
  , fContributors(0.5)
  , fWeightBackGround(0.)
  , fVz(0.0)
  , fContainer(NULL)
  , fVarManager(NULL)
  , fSignalCuts(NULL)
  , fCFM(NULL)
  , fTriggerAnalysis(NULL)
  , fPID(NULL)
  , fPIDqa(NULL)
  , fPIDpreselect(NULL)
  , fCuts(NULL)
  , fTaggedTrackCuts(NULL)
  , fCleanTaggedTrack(kFALSE)
  , fVariablesTRDTaggedTrack(kFALSE)
  , fCutspreselect(NULL)
  , fSecVtx(NULL)
  , fElecBackGround(NULL)
  , fMCQA(NULL)
  , fTaggedTrackAnalysis(NULL)
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
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  fPID = new AliHFEpid("hfePid");
  fPIDqa = new AliHFEpidQAmanager;
  fVarManager = new AliHFEvarManager("hfeVarManager");
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref):
  AliAnalysisTaskSE(ref)
  , fQAlevel(0)
  , fPlugins(0)
  , fFillSignalOnly(ref.fFillSignalOnly)
  , fBackGroundFactorApply(ref.fBackGroundFactorApply)
  , fRemovePileUp(ref.fRemovePileUp)
  , fIdentifiedAsPileUp(ref.fIdentifiedAsPileUp)
  , fIdentifiedAsOutInz(ref.fIdentifiedAsOutInz)
  , fPassTheEventCut(ref.fPassTheEventCut)
  , fHasSpecialTriggerSelection(ref.fHasSpecialTriggerSelection)
  , fSpecialTrigger(ref.fSpecialTrigger)
  , fCentralityF(ref.fCentralityF)
  , fContributors(ref.fContributors)
  , fWeightBackGround(ref.fWeightBackGround)
  , fVz(ref.fVz)
  , fContainer(NULL)
  , fVarManager(NULL)
  , fSignalCuts(NULL)
  , fCFM(NULL)
  , fTriggerAnalysis(NULL)
  , fPID(NULL)
  , fPIDqa(NULL)
  , fPIDpreselect(NULL)
  , fCuts(NULL)
  , fTaggedTrackCuts(NULL)
  , fCleanTaggedTrack(ref.fCleanTaggedTrack)
  , fVariablesTRDTaggedTrack(ref.fVariablesTRDTaggedTrack)
  , fCutspreselect(NULL)
  , fSecVtx(NULL)
  , fElecBackGround(NULL)
  , fMCQA(NULL)
  , fTaggedTrackAnalysis(NULL)
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
  target.fPlugins = fPlugins;
  target.fFillSignalOnly = fFillSignalOnly;
  target.fBackGroundFactorApply = fBackGroundFactorApply;
  target.fRemovePileUp = fRemovePileUp;
  target.fIdentifiedAsPileUp = fIdentifiedAsPileUp;
  target.fIdentifiedAsOutInz = fIdentifiedAsOutInz;
  target.fPassTheEventCut = fPassTheEventCut;
  target.fHasSpecialTriggerSelection = fHasSpecialTriggerSelection;
  target.fSpecialTrigger = fSpecialTrigger;
  target.fCentralityF = fCentralityF;
  target.fContributors = fContributors;
  target.fWeightBackGround = fWeightBackGround;
  target.fVz = fVz;
  target.fContainer = fContainer;
  target.fVarManager = fVarManager;
  target.fSignalCuts = fSignalCuts;
  target.fCFM = fCFM;
  target.fTriggerAnalysis = fTriggerAnalysis;
  target.fPID = fPID;
  target.fPIDqa = fPIDqa;
  target.fPIDpreselect = fPIDpreselect;
  target.fCuts = fCuts;
  target.fTaggedTrackCuts = fTaggedTrackCuts;
  target.fCleanTaggedTrack = fCleanTaggedTrack;
  target.fVariablesTRDTaggedTrack = fVariablesTRDTaggedTrack;
  target.fCutspreselect = fCutspreselect;
  target.fSecVtx = fSecVtx;
  target.fElecBackGround = fElecBackGround;
  target.fMCQA = fMCQA;
  target.fTaggedTrackAnalysis = fTaggedTrackAnalysis;
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
  if(fPID) delete fPID;
  if(fVarManager) delete fVarManager;
  if(fPIDqa) delete fPIDqa;
  if(fSignalCuts) delete fSignalCuts;
  if(fCFM) delete fCFM;
  if(fSecVtx) delete fSecVtx;
  if(fMCQA) delete fMCQA;
  if(fElecBackGround) delete fElecBackGround;
  if(fTriggerAnalysis) delete fTriggerAnalysis;
  if(fPIDpreselect) delete fPIDpreselect;
  if(fQA) delete fQA;
  if(fOutput) delete fOutput;
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

  // Enable Trigger Analysis
  fTriggerAnalysis = new AliTriggerAnalysis;
  fTriggerAnalysis->EnableHistograms();
  fTriggerAnalysis->SetAnalyzeMC(HasMCData());


  // Make lists for Output
  if(!fQA) fQA = new TList;
  fQA->SetOwner();
  if(!fOutput) fOutput = new TList;
  fOutput->SetOwner();

  // First Part: Make QA histograms
  fQACollection = new AliHFEcollection("TaskQA", "QA histos from the Electron Task");
  fQACollection->CreateTH1F("nElectronTracksEvent", "Number of Electron Candidates", 100, 0, 100);
  fQACollection->CreateProfile("conr", "Electron PID contamination", 20, 0, 20);
  fQACollection->CreateTH1F("alpha_rec", "Alpha from reconstructed tracks with TRD hits", 36, -TMath::Pi(), TMath::Pi());
  fQACollection->CreateTH1F("alpha_sim", "Alpha from simulated electron tracks", 36, -TMath::Pi(), TMath::Pi());
  fQACollection->CreateTH1F("nElectron", "Number of electrons", 100, 0, 100);
  fQACollection->CreateProfile("pidquality", "TRD PID quality as function of momentum", 20, 0, 20);
  fQACollection->CreateProfile("ntrdclusters", "Number of TRD clusters as function of momentum", 20, 0, 20);
  fQACollection->CreateTH1F("chi2TRD","#chi2 per TRD cluster", 20, 0, 20);
  fQACollection->CreateTH1F("mccharge", "MC Charge", 200, -100, 100);
  fQACollection->CreateTH2F("radius", "Production Vertex", 100, 0.0, 5.0, 100, 0.0, 5.0);
  // Temporary histograms for TPC number of clusters for all signal tracks (MC true electrons) and for selected tracks (Markus Fasel)
  fQACollection->CreateTH2F("TPCclusters2_1_Signal", "TPCclusterInfo for findable clusters for 2 neighbors for signal tracks", 30, 0.1, 10., 162, 0., 161.);
  fQACollection->CreateTH2F("TPCclusters2_0_Signal", "TPCclusterInfo for the ratio for 2 neighbors for signal tracks", 30, 0.1, 10., 100, 0., 1.);
  fQACollection->CreateTH2F("TPCclusters2_1_Selected", "TPCclusterInfo for findable clusters for 2 neighbors for selected tracks", 30, 0.1, 10., 162, 0., 161.);
  fQACollection->CreateTH2F("TPCclusters2_0_Selected", "TPCclusterInfo for the ratio for 2 neighbors for selected tracks", 30, 0.1, 10., 110, 0., 1.1);
  fQACollection->CreateTH2F("TPCncls_Signal", "TPC Number of clusters for signal tracks", 30, 0.1, 10., 162, 0., 161.);
  fQACollection->CreateTH2F("TPCclr_Signal", "TPC cluster ratio for signal tracks", 30, 0.1, 10., 110, 0., 1.1);
  fQACollection->BinLogAxis("TPCclusters2_1_Signal", 0); 
  fQACollection->BinLogAxis("TPCclusters2_0_Signal", 0);
  fQACollection->BinLogAxis("TPCclusters2_1_Selected", 0); 
  fQACollection->BinLogAxis("TPCclusters2_0_Selected", 0);
  fQACollection->BinLogAxis("TPCncls_Signal", 0); 
  fQACollection->BinLogAxis("TPCclr_Signal", 0);

  InitPIDperformanceQA();
  InitContaminationQA();
  fQA->Add(fQACollection->GetList());

  // Initialize PID
  fPID->SetHasMCData(HasMCData());
  if(!fPID->GetNumberOfPIDdetectors()) fPID->AddDetector("TPC", 0);
  fPID->InitializePID();
  if(IsQAOn(kPIDqa)){
    AliInfo("PID QA switched on");
    fPIDqa->Initialize(fPID);
    fQA->Add(fPIDqa->MakeList("HFEpidQA"));
  }

  // Initialize correction Framework and Cuts
  const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack + AliHFEcuts::kNcutStepsSecvtxTrack;
  fCFM = new AliCFManager;
  fCFM->SetNStepParticle(kNcutSteps);
  MakeParticleContainer();
  MakeEventContainer();
  // Temporary fix: Initialize particle cuts with NULL
  for(Int_t istep = 0; istep < kNcutSteps; istep++)
    fCFM->SetParticleCutsList(istep, NULL);
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  if(IsAODanalysis()) fCuts->SetAOD();
  // Make clone for V0 tagging step
  fCuts->Initialize(fCFM);
  if(fCuts->IsQAOn()) fQA->Add(fCuts->GetQAhistograms());
  fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
  fVarManager->SetSignalCuts(fSignalCuts);
 
  // add output objects to the List
  fOutput->AddAt(fContainer, 0);
  fOutput->AddAt(fCFM->GetEventContainer(), 1);
  
  // mcQA----------------------------------
  if (HasMCData() && IsQAOn(kMCqa)) {
    AliInfo("MC QA on");
    if(!fMCQA) fMCQA = new AliHFEmcQA;
    if(!fHistMCQA) fHistMCQA = new TList();
    fHistMCQA->SetOwner();
    fMCQA->CreatDefaultHistograms(fHistMCQA);
    fQA->Add(fHistMCQA);
  } 

  // secvtx----------------------------------
  if (GetPlugin(kSecVtx)) {
    AliInfo("Secondary Vertex Analysis on");
    if(!fSecVtx) fSecVtx = new AliHFEsecVtx;
    fSecVtx->SetHasMCData(HasMCData());

    if(!fHistSECVTX) fHistSECVTX = new TList();
    fHistSECVTX->SetOwner();
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
    fHistELECBACKGROUND->SetOwner();
    fElecBackGround->CreateHistograms(fHistELECBACKGROUND);
    fOutput->Add(fHistELECBACKGROUND);
  }  

  // tagged tracks
  if(GetPlugin(kTaggedTrackAnalysis)){
    AliInfo("Analysis on V0-tagged tracks enabled");
    fTaggedTrackAnalysis = new AliHFEtaggedTrackAnalysis(Form("taggedTrackAnalysis%s", GetName()));
    fTaggedTrackAnalysis->SetCuts(fTaggedTrackCuts);
    fTaggedTrackAnalysis->SetClean(fCleanTaggedTrack);
    AliHFEvarManager *varManager = fTaggedTrackAnalysis->GetVarManager();
    TObjArray *array = fVarManager->GetVariables();
    Int_t nvars = array->GetEntriesFast();
    TString namee;
    for(Int_t v = 0; v < nvars; v++) {
      AliHFEvarManager::AliHFEvariable *variable = (AliHFEvarManager::AliHFEvariable *) array->At(v);
      if(!variable) continue;
      TString name(((AliHFEvarManager::AliHFEvariable *)variable)->GetName());
      if(!name.CompareTo("source")) namee = TString("species");
      else namee = TString(name);
      varManager->AddVariable(namee);
      //printf("For AliTaggedTrackAnalysis, had the variable %s and the one used %s\n",(const char*)variable->GetName(),(const char*) namee);
    }
    if(fPIDqa->HasHighResolutionHistos()) 
      fTaggedTrackAnalysis->GetPIDqa()->SetHighResolutionHistos();
    fTaggedTrackAnalysis->SetPID(fPID);
    fTaggedTrackAnalysis->SetVariablesTRD(fVariablesTRDTaggedTrack);
    fTaggedTrackAnalysis->InitContainer();
    fOutput->Add(fTaggedTrackAnalysis->GetContainer());
    fQA->Add(fTaggedTrackAnalysis->GetPIDQA());
    fQA->Add(fTaggedTrackAnalysis->GetCutQA());
    fQA->Add(fTaggedTrackAnalysis->GetQAcollection());
  }
  PrintStatus();
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

  // need the centrality for everything (MC also)
  fCentralityF = -100.0;
  if(!ReadCentrality()) fCentralityF = -100.0;
  //printf("pass centrality\n");
  //printf("Reading fCentralityF %f\n",fCentralityF);
  
  // See if pile up and z in the range
  RejectionPileUpVertexRangeEventCut();

  // Protect agains missing 
  if(HasMCData()){
    //printf("Has MC data\n");
    fSignalCuts->SetMCEvent(fMCEvent);
    ProcessMC();  // Run the MC loop + MC QA in case MC Data are available
  }

  if(IsAODanalysis()){
    AliAODpidUtil *aodworkingpid = AliHFEtools::GetDefaultAODPID(HasMCData());
    fPID->SetAODpid(aodworkingpid); 
    ProcessAOD();
  } else {
    // Check Trigger selection
    if(fHasSpecialTriggerSelection){
      AliESDEvent *ev = dynamic_cast<AliESDEvent *>(fInputEvent);
      if(!(ev && ev->IsTriggerClassFired(fSpecialTrigger.Data()))) return;
    }
    AliESDInputHandler *inH = dynamic_cast<AliESDInputHandler *>(fInputHandler);
    if(!inH){
      AliError("No ESD Input handler available");
      return;
    }
    AliESDpid *workingPID = inH->GetESDpid();
    if(!workingPID){
      AliDebug(1, "Using default ESD PID");
      workingPID = AliHFEtools::GetDefaultPID(HasMCData());
    } else { 
      AliDebug(1, "Using ESD PID from the input handler");
    }
    fPID->SetESDpid(workingPID);
    if(fPIDpreselect) fPIDpreselect->SetESDpid(workingPID);
    
    ProcessESD();
  }
  // Done!!!
  PostData(1, fOutput);
  PostData(2, fQA);
}

//____________________________________________________________
void AliAnalysisTaskHFE::Terminate(Option_t *){
  //
  // Terminate not implemented at the moment
  //
  if(GetPlugin(kPostProcess)){
    fOutput = dynamic_cast<TList *>(GetOutputData(1));
    fQA = dynamic_cast<TList *>(GetOutputData(2));
    if(!fOutput){
      AliError("Results not available");
      return;
    }
    if(!fQA){
      AliError("QA output not available");
      return;
    }
    fContainer = dynamic_cast<AliHFEcontainer *>(fOutput->FindObject("trackContainer")); 
    if(!fContainer){
      AliError("Track container not found");
      return;
    }
    AliHFEpostAnalysis postanalysis;
    postanalysis.SetTaskResults(fContainer);
    TList *qalist = dynamic_cast<TList *>(fQA->FindObject("list_TaskQA"));
    if(!qalist){
      AliError("QA List not found");
      return;
    }
    postanalysis.SetTaskQA(qalist);
    printf("Running post analysis\n");
    //if(HasMCData())
    postanalysis.DrawMCSignal2Background();
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
  Double_t eventContainer [4];
  eventContainer[0] = fMCEvent->GetPrimaryVertex()->GetZ();
  eventContainer[2] = fCentralityF;
  eventContainer[3] = fContributors;
  fVz = eventContainer[0];
  //printf("z position is %f\n",eventContainer[0]);
  //if(fCFM->CheckEventCuts(AliHFEcuts::kEventStepGenerated, fMCEvent)) 
  fCFM->GetEventContainer()->Fill(eventContainer,AliHFEcuts::kEventStepGenerated);
  Int_t nElectrons = 0;
  if(IsESDanalysis()){
   if(!((fIdentifiedAsPileUp) || (TMath::Abs(fVz) > fCuts->GetVertexRange()) || (TMath::Abs(fCentralityF+100.0) < 0.01))){ //kStepMCGeneratedZOutNoPileUpCentralityFine
    if (HasMCData() && IsQAOn(kMCqa)) {
      AliDebug(2, "Running MC QA");

      if(fMCEvent->Stack()){
	      fMCQA->SetMCEvent(fMCEvent);
        fMCQA->SetGenEventHeader(fMCEvent->GenEventHeader());
        fMCQA->Init();
       
        fMCQA->GetMesonKine();

        // loop over all tracks for decayed electrons
        for (Int_t igen = 0; igen < fMCEvent->GetNumberOfTracks(); igen++){
          TParticle* mcpart = fMCEvent->Stack()->Particle(igen);
	  if(!mcpart) continue;
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
   }
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
    AliDebug(4, "Next MC Track");
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

  // Set magnetic field if V0 task on
  if(fTaggedTrackAnalysis) {
    fTaggedTrackAnalysis->SetMagneticField(fESD->GetMagneticField());
    fTaggedTrackAnalysis->SetCentrality(fCentralityF);
  }

  // Do event Normalization
  Double_t eventContainer[4];
  eventContainer[0] = 0.0;
  if(HasMCData()) eventContainer[0] = fVz;
  else {
    if(fESD->GetPrimaryVertexTracks()) eventContainer[0] = fESD->GetPrimaryVertexTracks()->GetZ();
  }
  eventContainer[1] = 0.;
  eventContainer[2] = fCentralityF;
  eventContainer[3] = fContributors;
  if(fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0AND))
    eventContainer[1] = 1.;

  //
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoCut);

  //
  if(fIdentifiedAsPileUp) return; 
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoPileUp);

  //
  if(TMath::Abs(fCentralityF+100.0) < 0.01) return; 
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecCentralityOk);
  //printf("In ProcessESD %f\n",fCentralityF);

  //
  if(fIdentifiedAsOutInz) return;
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepZRange);  

  //
  if(!fPassTheEventCut) return;
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepReconstructed);

 

  fContainer->NewEvent();

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
      fSecVtx->SetMCQA(fMCQA); 
    }
    if (GetPlugin(kIsElecBackGround)) { 
      fElecBackGround->SetMCEvent(fMCEvent);
    }
  }

  Double_t container[10];
  memset(container, 0, sizeof(Double_t) * 10);
  // container for the output THnSparse
  Double_t dataE[6]; // [pT, eta, Phi, type, 'C' or 'B']
  Int_t nElectronCandidates = 0;
  AliESDtrack *track = NULL, *htrack = NULL;
  AliMCParticle *mctrack = NULL;
  TParticle* mctrack4QA = NULL;
  Int_t pid = 0;

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
    AliDebug(4, "New ESD track");
    track = fESD->GetTrack(itrack);
    track->SetESDEvent(fESD);

    // fill counts of v0-identified particles
    Int_t v0pid = -1;
    if(track->TestBit(BIT(14))) v0pid = AliPID::kElectron;
    else if(track->TestBit(BIT(15))) v0pid = AliPID::kPion;
    else if(track->TestBit(BIT(16))) v0pid = AliPID::kProton;
    // here the tagged track analysis will run
    if(fTaggedTrackAnalysis && v0pid > -1){ 
      AliDebug(1, Form("Track identified as %s", AliPID::ParticleName(v0pid)));
      fTaggedTrackAnalysis->ProcessTrack(track, v0pid);
    }
 
    AliDebug(3, Form("Doing track %d, %p", itrack, track));
     
    //////////////////////////////////////
    // preselect
    /////////////////////////////////////
    if(fPIDpreselect && fCutspreselect) {
      if(!PreSelectTrack(track)) continue;
    }

    signal = kTRUE;
    
    // Fill step without any cut
          
    if(HasMCData()){
      // Check if it is electrons near the vertex
      if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(track->GetLabel()))))) continue;
      mctrack4QA = mctrack->Particle();

      if(fFillSignalOnly && !fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE; 
      else AliDebug(3, "Signal Electron");
    } 
    // Cache new Track information inside the var manager
    fVarManager->NewTrack(track, mctrack, fCentralityF, -1, signal);

    if(signal) {
      fVarManager->FillContainer(fContainer, "recTrackContReco", AliHFEcuts::kStepRecNoCut, kFALSE);
      fVarManager->FillContainer(fContainer, "recTrackContMC", AliHFEcuts::kStepRecNoCut, kTRUE);
      if((track->GetStatus() & AliESDtrack::kTPCout) 
          && (TMath::Abs(track->Eta()) < 0.8) 
          && (track->GetKinkIndex(0) == 0)){
        fQACollection->Fill("TPCclusters2_1_Signal", track->Pt(), track->GetTPCClusterInfo(2,1));
        fQACollection->Fill("TPCclusters2_0_Signal", track->Pt(), track->GetTPCNclsF() > 0 ?  track->GetTPCClusterInfo(2,1)/track->GetTPCNclsF() : 0.);
        fQACollection->Fill("TPCncls_Signal", track->Pt(), track->GetTPCNcls());
        fQACollection->Fill("TPCclr_Signal", track->Pt(), track->GetTPCNclsF() > 0 ? static_cast<Double_t>(track->GetTPCNcls())/static_cast<Double_t>(track->GetTPCNclsF()) : 0.);
      }
    }

    // RecKine: ITSTPC cuts  
    if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    
    // Check TRD criterions (outside the correction framework)
    if(track->GetTRDncls()){
      fQACollection->Fill("chi2TRD", track->GetTRDchi2()/track->GetTRDncls());
      fQACollection->Fill("alpha_rec", track->GetAlpha());
      fQACollection->Fill("pidquality", container[0], track->GetTRDpidQuality());
      fQACollection->Fill("ntrdclusters", container[0], track->GetTRDncls());
    }

    
    // RecPrim
    if(track->GetKinkIndex(0) != 0) continue; // Quick and dirty fix to reject both kink mothers and daughters
    if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;

    // HFEcuts: ITS layers cuts
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;

    if(HasMCData() && IsQAOn(kMCqa)) {
      // mc qa for after the reconstruction cuts  
      AliDebug(2, "Running MC QA");
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 3);  // charm
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 3); // beauty 
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kOthers,  AliHFEmcQA::kElectronPDG, 3); // beauty 
    }

    // HFEcuts: Nb of tracklets TRD0
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTRD, track)) continue;

    // Fill correlation maps before PID
    if(signal && fContainer->GetCorrelationMatrix("correlationstepbeforePID")) {
      //printf("Fill correlation maps before PID\n");
      fVarManager->FillCorrelationMatrix(fContainer->GetCorrelationMatrix("correlationstepbeforePID"));
    }

    if (HasMCData() && IsQAOn(kMCqa)) {
      // mc qa for after the reconstruction and pid cuts  
      AliDebug(2, "Running MC QA");
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 4);  // charm
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 4); // beauty 
      fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kOthers,  AliHFEmcQA::kElectronPDG, 4); // beauty 
    }

    if(HasMCData()){
      FillProductionVertex(track);
    }

    // track accepted, do PID
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    if(HasMCData()) hfetrack.SetMCTrack(mctrack);
    hfetrack.SetCentrality(fCentralityF);
    fPID->SetVarManager(fVarManager);
    if(!fPID->IsSelected(&hfetrack, fContainer, "recTrackCont", fPIDqa)) continue;
    nElectronCandidates++;
    fQACollection->Fill("TPCclusters2_1_Selected", track->Pt(), track->GetTPCClusterInfo(2,1));
    fQACollection->Fill("TPCclusters2_0_Selected", track->Pt(), track->GetTPCClusterInfo(2,0));

    // Fill Histogram for Hadronic Background
    if(HasMCData()){
      if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) != 11))
        fVarManager->FillContainer(fContainer, "hadronicBackground", UInt_t(0), kFALSE);
    }

    // Fill Containers
    if(signal) {
      // Apply weight for background contamination
	    if(fBackGroundFactorApply==kTRUE) {
	      if(IsPbPb()) fWeightBackGround =  fBackGroundFactorArray[(Int_t)fCentralityF]->Eval(TMath::Abs(track->P()));
	      else    fWeightBackGround =  fBackGroundFactorArray[0]->Eval(TMath::Abs(track->P())); // pp case

	      if(fWeightBackGround < 0.0) fWeightBackGround = 0.0;
	      else if(fWeightBackGround > 1.0) fWeightBackGround = 1.0;
        // weightBackGround as special weight
        fVarManager->FillContainer(fContainer, "hadronicBackground", 1, kFALSE, fWeightBackGround);
      }
      fVarManager->FillCorrelationMatrix(fContainer->GetCorrelationMatrix("correlationstepafterPID"));
    }

    if(GetPlugin(kSecVtx)) {
      AliDebug(2, "Running Secondary Vertex Analysis");
      if(fSecVtx->Process(track) && signal) {
        fVarManager->FillContainer(fContainer, "recTrackContSecvtxReco", AliHFEcuts::kStepHFEcutsSecvtx, kFALSE);
        fVarManager->FillContainer(fContainer, "recTrackContSecvtxMC", AliHFEcuts::kStepHFEcutsSecvtx, kTRUE);
      }
    }

    if(HasMCData()){
      dataE[0] = track->Pt();
      dataE[1] = track->Eta();
      dataE[2] = track->Phi();
      dataE[3] = track->Charge();
      dataE[4] = -1.;
      dataE[5] = -1.;

      // Track selected: distinguish between true and fake
      AliDebug(1, Form("Candidate Selected, filling THnSparse, PID: %d\n", mctrack->Particle()->GetPdgCode()));
      if((pid = TMath::Abs(mctrack->Particle()->GetPdgCode())) == 11){
        Int_t type = 0;
        if(fSignalCuts->IsCharmElectron(track))
          type = 1;
        else if(fSignalCuts->IsBeautyElectron(track))
          type = 2;
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
        fQACollection->Fill("PIDperformance", dataE);
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

    if (GetPlugin(kDEstep)) { 
      if(HasMCData()){
        if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) != 11)){
          fQACollection->Fill("hadronsBeforeIPcut",track->Pt());
          fQACollection->Fill("hadronsBeforeIPcutMC",mctrack->Pt());
        }
      }
      // Fill Containers for impact parameter analysis
      if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsDca + AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack,track)) continue;

      if(signal) {
        fVarManager->FillContainer(fContainer, "recTrackContDEReco", AliHFEcuts::kStepHFEcutsDca, kFALSE);
        fVarManager->FillContainer(fContainer, "recTrackContDEMC", AliHFEcuts::kStepHFEcutsDca, kTRUE);
        fVarManager->FillCorrelationMatrix(fContainer->GetCorrelationMatrix("correlationstepafterDE"));
      }
      if(HasMCData()){
        if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) != 11)){
          fQACollection->Fill("hadronsAfterIPcut",track->Pt());
          fQACollection->Fill("hadronsAfterIPcutMC",mctrack->Pt());
        }
      }
    }

  }
  fQACollection->Fill("nElectronTracksEvent", nElectronCandidates);
}

//____________________________________________________________
void AliAnalysisTaskHFE::ProcessAOD(){
  //
  // Run Analysis in AOD Mode
  // Function is still in development
  //
  AliDebug(3, "Processing AOD Event");
  Double_t eventContainer[4];
  if(HasMCData()) eventContainer[0] = fVz;
  else {
    eventContainer[0] = fInputEvent->GetPrimaryVertex()->GetZ();
  }
  eventContainer[1] = 1.; // No Information available in AOD analysis, assume all events have V0AND
  eventContainer[2] = fCentralityF; 
  eventContainer[3] = fContributors; 
  
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
  if(!fAOD){
    AliError("AOD Event required for AOD Analysis")
      return;
  }
  
  //
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoCut);

  //
  if(fIdentifiedAsPileUp) return; 
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoPileUp);

  //
  if(fIdentifiedAsOutInz) return;
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepZRange);  

  //
  if(!fPassTheEventCut) return;
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepReconstructed);

  fContainer->NewEvent();
 
  AliAODTrack *track = NULL;
  AliAODMCParticle *mctrack = NULL;
  Double_t dataE[6]; // [pT, eta, Phi, Charge, type, 'C' or 'B']
  Int_t nElectronCandidates = 0;
  Int_t pid;
  Bool_t signal;
  for(Int_t itrack = 0; itrack < fAOD->GetNumberOfTracks(); itrack++){
    track = fAOD->GetTrack(itrack);
    if(!track) continue;
    if(track->GetFlags() != 1<<4) continue;  // Only process AOD tracks where the HFE is set

    signal = kTRUE;
    if(HasMCData()){

      Int_t label = TMath::Abs(track->GetLabel());
      if(label)
        mctrack = dynamic_cast<AliAODMCParticle *>(fMCEvent->GetTrack(label));
        if(fFillSignalOnly && !fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE;
    }
    fVarManager->NewTrack(track, mctrack, fCentralityF, -1, kTRUE);
    // track accepted, do PID
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
    hfetrack.SetRecTrack(track);
    if(HasMCData()) hfetrack.SetMCTrack(mctrack);
    hfetrack.SetCentrality(fCentralityF);
    fPID->SetVarManager(fVarManager);
    if(!fPID->IsSelected(&hfetrack, fContainer, "recTrackCont", fPIDqa)) continue;    // we will do PID here as soon as possible
    // Apply weight for background contamination
    Double_t weightBackGround = 1.0;

    // not correct treatment for pp
    if(fBackGroundFactorApply==kTRUE) {
	    if(IsPbPb()) weightBackGround =  fBackGroundFactorArray[(Int_t)fCentralityF]->Eval(TMath::Abs(track->P()));
      else weightBackGround =  fBackGroundFactorArray[0]->Eval(TMath::Abs(track->P()));

      if(weightBackGround < 0.0) weightBackGround = 0.0;
	    else if(weightBackGround > 1.0) weightBackGround = 1.0;
	    fVarManager->FillContainer(fContainer, "hadronicBackground", 1, kFALSE, weightBackGround);
    }



    nElectronCandidates++;    
    if(HasMCData()){
      dataE[0] = track->Pt();
      dataE[1] = track->Eta();
      dataE[2] = track->Phi();
      dataE[3] = track->Charge();
      dataE[4] = -1;
      dataE[5] = -1;
      // Track selected: distinguish between true and fake
      // COVERITY: missing test if mctrack != 0
      AliDebug(1, Form("Candidate Selected, filling THnSparse, PID: %d\n", mctrack->GetPdgCode()));
      if((pid = TMath::Abs(mctrack->GetPdgCode())) == 11){
      
        Int_t type = 0;
        if(fSignalCuts->IsCharmElectron(track))
          type = 1;
        else if(fSignalCuts->IsBeautyElectron(track))
          type = 2;
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
        fQACollection->Fill("PIDperformance", dataE);
      }
    }
  }
  fQACollection->Fill("nElectronTracksEvent", nElectronCandidates);
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFE::ProcessMCtrack(AliVParticle *track){
  //
  // Filter the Monte Carlo Track
  // Additionally Fill a THnSparse for Signal To Background Studies
  // Works for AOD and MC analysis Type
  //
  fVarManager->NewTrack(track, NULL, fCentralityF, -1, kTRUE);
  Double_t signalContainer[6];

  signalContainer[0] = track->Pt();
  signalContainer[1] = track->Eta();
  signalContainer[2] = track->Phi();
  signalContainer[3] = track->Charge()/3;

 Double_t vertex[3]; // Production vertex cut to mask gammas which are NOT supposed to have hits in the first ITS layer(s)
  if(IsESDanalysis()){
    AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(track);
    if(mctrack){
      vertex[0] = mctrack->Particle()->Vx();
      vertex[1] = mctrack->Particle()->Vy();
    }
  } else {
    AliAODMCParticle *aodmctrack = dynamic_cast<AliAODMCParticle *>(track);
    if(aodmctrack) aodmctrack->XvYvZv(vertex);
  }

  if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, track)) return kFALSE;
  fQACollection->Fill("mccharge", signalContainer[3]);
  fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCGenerated, kFALSE);
  signalContainer[4] = 0;
  if(fSignalCuts->IsSelected(track)){
    //fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCsignal, kFALSE);
    // Filling of the Signal/Background histogram using the 
    // definition of the codes for charm and beauty as below in
    // th crearion of the histogram
    if(fSignalCuts->IsCharmElectron(track))
      signalContainer[4] = 1;
    else 
      signalContainer[4] = 2;
  } else {
    signalContainer[4] = 0; // (and other background)
  }
  signalContainer[5] = 0;
  // apply cut on the sqrt of the production vertex
  Double_t radVertex = TMath::Sqrt(vertex[0]*vertex[0] + vertex[1] * vertex[1]);
  if(radVertex < 3.5){
    // Within first ITS layer(2) -> Background we cannot reject by ITS cut, let it pass
    signalContainer[5] = 1;
  } else if (radVertex < 7.5){
    signalContainer[5] = 2;
  }
  fQACollection->Fill("SignalToBackgroundMC", signalContainer);
  fQACollection->Fill("alpha_sim", track->Phi() - TMath::Pi());

  // Step GeneratedZOutNoPileUp
  if((fIdentifiedAsPileUp) || (TMath::Abs(fVz) > fCuts->GetVertexRange()) || (TMath::Abs(fCentralityF+100.0) < 0.01)) return kFALSE;
  fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCGeneratedZOutNoPileUpCentralityFine, kFALSE);
  //printf("In ProcessMCtrack %f\n",fCentralityF);

  // Step Generated Event Cut
  if(!fPassTheEventCut) return kFALSE;
  fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCGeneratedEventCut, kFALSE);

  if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCInAcceptance, track)) return kFALSE;
  fVarManager->FillContainer(fContainer, "MCTrackCont", AliHFEcuts::kStepMCInAcceptance, kFALSE);
  return kTRUE;
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
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
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
  // 1st bin: Vertex z-position
  // 2nd bin: V0AND decision (normalization to sigma_inel)
  // 3rd bin: Centrality class (for pp defined as number of contributors in vertex.)
  //
  
  if(IsPbPb()) {

    //printf("This is PbPb!!!!!!!!!!!\n");

    const Int_t kNvar = 4;  // number of variables on the grid: 
    Int_t nBins[kNvar] = {120, 2, 11, 2};
    Double_t binMin[kNvar] = {-30. , 0., 0.0, 0.};
    Double_t binMax[kNvar] = {30., 2., 11.0, 2.};
    
    AliCFContainer *evCont = new AliCFContainer("eventContainer", "Container for events", AliHFEcuts::kNcutStepsEvent, kNvar, nBins);
    
    Double_t *vertexBins = AliHFEtools::MakeLinearBinning(nBins[0], binMin[0], binMax[0]);
    Double_t *v0andBins = AliHFEtools::MakeLinearBinning(nBins[1], binMin[1], binMax[1]);
    Double_t *centralityBins = AliHFEtools::MakeLinearBinning(nBins[2], binMin[2], binMax[2]);
    Double_t *contributorsBins = AliHFEtools::MakeLinearBinning(nBins[3], binMin[3], binMax[3]);
    evCont->SetBinLimits(0, vertexBins);
    evCont->SetBinLimits(1, v0andBins);
    evCont->SetBinLimits(2, centralityBins);
    evCont->SetBinLimits(3, contributorsBins);
    delete[] vertexBins; delete[] v0andBins; delete[] centralityBins; delete[] contributorsBins;
    
    fCFM->SetEventContainer(evCont);
  }
  else {

    //printf("This is pp!!!!!!!!!!!\n");

    const Int_t kNvar = 3;  // number of variables on the grid: 
    Int_t nBins[kNvar] = {120, 2, 11};
    Double_t binMin[kNvar] = {-30. , 0., 0.0};
    Double_t binMax[kNvar] = {30., 2., 11.0};
    
    AliCFContainer *evCont = new AliCFContainer("eventContainer", "Container for events", AliHFEcuts::kNcutStepsEvent, kNvar, nBins);
    
    Double_t *vertexBins = AliHFEtools::MakeLinearBinning(nBins[0], binMin[0], binMax[0]);
    Double_t *v0andBins = AliHFEtools::MakeLinearBinning(nBins[1], binMin[1], binMax[1]);
    Double_t *centralityBins = AliHFEtools::MakeLinearBinning(nBins[2], binMin[2], binMax[2]);
    evCont->SetBinLimits(0, vertexBins);
    evCont->SetBinLimits(1, v0andBins);
    evCont->SetBinLimits(2, centralityBins);
    delete[] vertexBins; delete[] v0andBins; delete[] centralityBins;
    
    fCFM->SetEventContainer(evCont);
  }
  
}

//____________________________________________________________
void AliAnalysisTaskHFE::MakeParticleContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
  //
  
  if(!fContainer) fContainer = new AliHFEcontainer("trackContainer");
  fVarManager->DefineVariables(fContainer);

  // Create Correction Framework containers
  fContainer->CreateContainer("MCTrackCont", "Track Container filled with MC information", AliHFEcuts::kNcutStepsMCTrack);
  fContainer->CreateContainer("recTrackContReco", "Track Container filled with MC information", AliHFEcuts::kNcutStepsRecTrack + fPID->GetNumberOfPIDdetectors());
  fContainer->CreateContainer("recTrackContMC", "Track Container filled with MC information", AliHFEcuts::kNcutStepsRecTrack + fPID->GetNumberOfPIDdetectors());
  
  fContainer->CreateContainer("hadronicBackground", "Container for Hadronic Background", 2);
  fContainer->CreateContainer("recTrackContDEReco", "Container for displaced electron analysis with Reco information", 1);
  fContainer->CreateContainer("recTrackContDEMC", "Container for displaced electron analysis with MC information", 1);
  fContainer->CreateContainer("recTrackContSecvtxReco", "Container for secondary vertexing analysis with Reco information", 1);
  fContainer->CreateContainer("recTrackContSecvtxMC", "Container for secondary vertexing analysis with MC information", 1);
  fContainer->CreateCorrelationMatrix("correlationstepafterPID","THnSparse with correlations");
  fContainer->CreateCorrelationMatrix("correlationstepafterDE","THnSparse with correlations");
  if(!fVarManager->IsVariableDefined("centrality")) {
    //printf("Create the two other correlation maps\n");
    fContainer->CreateCorrelationMatrix("correlationstepbeforePID","THnSparse with correlations");
    fContainer->CreateCorrelationMatrix("correlationstepafterTOF","THnSparse with correlations");
  }

  // Define the step names
  for(UInt_t istep = 0; istep < AliHFEcuts::kNcutStepsMCTrack; istep++){
    fContainer->SetStepTitle("MCTrackCont", AliHFEcuts::MCCutName(istep), istep);
  }
  for(UInt_t istep = 0; istep < AliHFEcuts::kNcutStepsRecTrack; istep++){
    fContainer->SetStepTitle("recTrackContReco", AliHFEcuts::RecoCutName(istep), istep);
    fContainer->SetStepTitle("recTrackContMC", AliHFEcuts::RecoCutName(istep), istep);
  }
  for(UInt_t ipid = 0; ipid < fPID->GetNumberOfPIDdetectors(); ipid++){
    fContainer->SetStepTitle("recTrackContReco", fPID->SortedDetectorName(ipid), AliHFEcuts::kNcutStepsRecTrack + ipid);
    fContainer->SetStepTitle("recTrackContMC", fPID->SortedDetectorName(ipid), AliHFEcuts::kNcutStepsRecTrack + ipid);
  }
}

//____________________________________________________________
void AliAnalysisTaskHFE::InitPIDperformanceQA(){
  // Add a histogram for Fake electrons
  const Int_t nDim=6;
  Int_t nBin[nDim] = {40, 8, 18, 2, 3, 3};
  //number of variables on the grid:pt,eta,phi,charge,
  const Double_t kPtbound[2] = {0.1, 20.};
  const Double_t kEtabound[2] = {-0.8, 0.8};
  const Double_t kPhibound[2] = {0., 2. * TMath::Pi()}; 
  const Double_t kChargebound[2] = {-1.1, 1.1};
  const Double_t kAddInf1bound[2] = {0., 3.};
  const Double_t kAddInf2bound[2] = {0., 3.};
  Double_t minima[nDim] = {kPtbound[0], kEtabound[0], kPhibound[0], kChargebound[0], kAddInf1bound[0], kAddInf2bound[0]}; 
  Double_t maxima[nDim] = {kPtbound[1], kEtabound[1], kPhibound[1], kChargebound[1], kAddInf1bound[1], kAddInf2bound[1]}; 
  
  fQACollection->CreateTHnSparse("PIDperformance", "PID performance; pT [GeV/c]; theta [rad]; phi [rad]; charge; type (0 - not el, 1 - other el, 2 - HF el; flavor (0 - no, 1 - charm, 2 - bottom)", nDim, nBin, minima, maxima);
  fQACollection->CreateTHnSparse("SignalToBackgroundMC", "PID performance; pT [GeV/c]; theta [rad]; phi [rad]; charge; flavor (0 - no, 1 - charm, 2 - bottom); ITS Cluster (0 - no, 1 - first (and maybe second), 2 - second)", nDim, nBin, minima, maxima);

  fQACollection->BinLogAxis("PIDperformance", 0);
  fQACollection->BinLogAxis("SignalToBackgroundMC", 0);
  fQACollection->Sumw2("PIDperformance");
  fQACollection->Sumw2("SignalToBackgroundMC");
}

//____________________________________________________________
void AliAnalysisTaskHFE::InitContaminationQA(){
  // 
  // Add QA for Impact Parameter cut
  //
  const Double_t kPtbound[2] = {0.1, 20.};
  Int_t iBin[1];
  iBin[0] = 44; // bins in pt
  fQACollection->CreateTH1F("hadronsBeforeIPcut", "Hadrons before IP cut", iBin[0], kPtbound[0], kPtbound[1], 1);
  fQACollection->CreateTH1F("hadronsAfterIPcut", "Hadrons after IP cut", iBin[0], kPtbound[0], kPtbound[1], 1);
  fQACollection->CreateTH1F("hadronsBeforeIPcutMC", "Hadrons before IP cut: MC p_{t}", iBin[0], kPtbound[0], kPtbound[1], 1);
  fQACollection->CreateTH1F("hadronsAfterIPcutMC", "Hadrons after IP cut: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
}

//____________________________________________________________
void AliAnalysisTaskHFE::PrintStatus() const {
  //
  // Print Analysis status
  //
  printf("\n\tAnalysis Settings\n\t========================================\n\n");
  printf("\tSecondary Vertex finding: %s\n", GetPlugin(kSecVtx) ? "YES" : "NO");
  printf("\tPrimary Vertex resolution: %s\n", GetPlugin(kPriVtx) ? "YES" : "NO");
  printf("\tDisplaced electron analysis step: %s\n", GetPlugin(kDEstep) ? "YES" : "NO");
  printf("\tTagged Track Analysis: %s\n", GetPlugin(kTaggedTrackAnalysis) ? "YES" : "NO");
  printf("\n");
  printf("\tParticle Identification Detectors:\n");
  fPID->PrintStatus();
  printf("\n");
  printf("\tQA: \n");
  printf("\t\tPID: %s\n", IsQAOn(kPIDqa) ? "YES" :  "NO");
  printf("\t\tCUTS: %s\n", (fCuts != NULL && fCuts->IsQAOn()) ? "YES" : "NO");
  printf("\t\tMC: %s\n", IsQAOn(kMCqa) ? "YES" : "NO");
  printf("\n");
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
    const AliMCParticle *mcpart = dynamic_cast<const AliMCParticle *>(mctrack);
    if(mcpart){
      xv =  mcpart->Xv();
      yv =  mcpart->Yv();
    }
  } else {
    // case AODMCParticle
    const AliAODMCParticle *mcpart = dynamic_cast<const AliAODMCParticle *>(mctrack);
    if(mcpart){
      xv =  mcpart->Xv();
      yv =  mcpart->Yv();
    }
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
    case kDEstep: SETBIT(fPlugins, plug); break;
    case kTaggedTrackAnalysis: SETBIT(fPlugins, plug); break;
    default: AliError("Unknown Plugin");
  };
}
//__________________________________________
Bool_t AliAnalysisTaskHFE::ProcessCutStep(Int_t cutStep, AliVParticle *track){
  //
  // Check single track cuts for a given cut step
  // Fill the particle container
  //
  const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
  if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
  if(fVarManager->IsSignalTrack()) {
    fVarManager->FillContainer(fContainer, "recTrackContReco", cutStep, kFALSE);
    fVarManager->FillContainer(fContainer, "recTrackContMC", cutStep, kTRUE);
  }
  return kTRUE;
}
//___________________________________________________
Bool_t AliAnalysisTaskHFE::ReadCentrality() {
  //
  // Recover the centrality of the event from ESD or AOD
  //
 if(IsAODanalysis()){

   AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
   if(!fAOD){
     AliError("AOD Event required for AOD Analysis")
       return kFALSE;
   }

   if(IsPbPb()) {
     // Centrality
     AliCentrality *aodCentrality = fAOD->GetCentrality();
     Float_t fCentralityFtemp = aodCentrality->GetCentralityPercentile("V0M");
     
     if( fCentralityFtemp >=  0. && fCentralityFtemp <  10.) fCentralityF =  0;
     else if ( fCentralityFtemp >=  10. && fCentralityFtemp <  20.) fCentralityF =  1;
     else if ( fCentralityFtemp >= 20. && fCentralityFtemp <  30.) fCentralityF = 2;
     else if ( fCentralityFtemp >= 30. && fCentralityFtemp <  40.) fCentralityF = 3;
     else if ( fCentralityFtemp >= 40. && fCentralityFtemp <  50.) fCentralityF = 4;
     else if ( fCentralityFtemp >= 50. && fCentralityFtemp <  60.) fCentralityF = 5;
     else if ( fCentralityFtemp >= 60. && fCentralityFtemp <  90.) fCentralityF = 6;
     else if ( fCentralityFtemp >= 90. && fCentralityFtemp <=  100.) fCentralityF = 7;
     //else if ( fCentralityF_temp >= 90. && fCentralityF_temp <  95.) fCentralityF = 8;
     //else if ( fCentralityF_temp >= 95. && fCentralityF_temp <  90.) fCentralityF = 9;
     //else if ( fCentralityF_temp >= 90. && fCentralityF_temp <=100.) fCentralityF = 10;
     else return kFALSE;
 
     // contributors
     fContributors = 0.5;
     Int_t contributorstemp = 0;
     const AliAODVertex *vtxAOD = fAOD->GetPrimaryVertex();
     if(vtxAOD) contributorstemp = vtxAOD->GetNContributors();
     
     //printf("PbPb contributors_temp %d\n",contributors_temp);
     
     if( contributorstemp <=  0) fContributors =  0.5;
     else fContributors = 1.5;   
   


   }
   else {
     fCentralityF = 0;
     Int_t centralityFtemp = 0;
     const AliAODVertex *vtxAOD = fAOD->GetPrimaryVertex();
     if(vtxAOD) centralityFtemp = vtxAOD->GetNContributors();
     
     //printf("pp centralityF_temp %d\n",centralityF_temp);
     
     if( centralityFtemp <=  0) fCentralityF =  0;
     else if ( centralityFtemp >  0 && centralityFtemp <  2) fCentralityF = 1;
     else if ( centralityFtemp >=  2 && centralityFtemp <  3) fCentralityF = 2;
     else if ( centralityFtemp >= 3 && centralityFtemp <  4) fCentralityF = 3;
     else if ( centralityFtemp >= 4 && centralityFtemp <  5) fCentralityF = 4;
     else if ( centralityFtemp >= 5 && centralityFtemp <  10) fCentralityF = 5;
     else if ( centralityFtemp >= 10 && centralityFtemp <  20) fCentralityF = 6;
     else if ( centralityFtemp >= 20 && centralityFtemp <  30) fCentralityF = 7;
     else if ( centralityFtemp >= 30 && centralityFtemp <  40) fCentralityF = 8;
     else if ( centralityFtemp >= 40 && centralityFtemp <  50) fCentralityF = 9;
     else if ( centralityFtemp >= 50) fCentralityF = 10;
     
   }

   return kTRUE;

   
 } else {
   
   AliDebug(3, "Processing ESD Centrality");
   AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
   if(!fESD){
     AliError("ESD Event required for ESD Analysis")
       return kFALSE;
   }
   const char* type = fESD->GetBeamType();

   if (strstr(type,"Pb-Pb")) {
   
   // Centrality
   AliCentrality *esdCentrality = fESD->GetCentrality();
   Float_t fCentralityFtemp = esdCentrality->GetCentralityPercentile("V0M");
   //printf("PbPb fCentralityF_temp %f\n",fCentralityF_temp);

   if( fCentralityFtemp >=  0. && fCentralityFtemp <   10.) fCentralityF =  0;
   else if ( fCentralityFtemp >= 10. && fCentralityFtemp <  20.) fCentralityF =  1;
   else if ( fCentralityFtemp >= 20. && fCentralityFtemp <  30.) fCentralityF = 2;
   else if ( fCentralityFtemp >= 30. && fCentralityFtemp <  40.) fCentralityF = 3;
   else if ( fCentralityFtemp >= 40. && fCentralityFtemp <  50.) fCentralityF = 4;
   else if ( fCentralityFtemp >= 50. && fCentralityFtemp <  60.) fCentralityF = 5;
   else if ( fCentralityFtemp >= 60. && fCentralityFtemp <  90.) fCentralityF = 6;
   else if ( fCentralityFtemp >= 90. && fCentralityFtemp <=  100.) fCentralityF = 7;
   //else if ( fCentralityF_temp >= 70. && fCentralityF_temp <  80.) fCentralityF = 8;
   //else if ( fCentralityF_temp >= 80. && fCentralityF_temp <  90.) fCentralityF = 9;
   //else if ( fCentralityF_temp >= 90. && fCentralityF_temp <=100.) fCentralityF = 10;
   else return kFALSE;

   //   Float_t fCentralityF_temp10 = esdCentrality->GetCentralityClass10("V0M");
   //   printf("PbPb fCentralityF_temp %f %f %f \n",fCentralityF_temp, fCentralityF_temp10, fCentralityF);

   // contributors
   fContributors = 0.5;
   Int_t contributorstemp = 0;
   const AliESDVertex *vtxESD = fESD->GetPrimaryVertexTracks();
   if(vtxESD && vtxESD->GetStatus()) contributorstemp = vtxESD->GetNContributors();
   
   //printf("PbPb contributors_temp %d\n",contributors_temp);
   
   if( contributorstemp <=  0) fContributors =  0.5;
   else fContributors = 1.5;   
   
   return kTRUE;

   }

   
   if (strstr(type,"p-p")) {
     fCentralityF = 0;
     Int_t centralityFtemp = 0;
     const AliESDVertex *vtxESD = fESD->GetPrimaryVertexTracks();
     if(vtxESD && vtxESD->GetStatus()) centralityFtemp = vtxESD->GetNContributors();
     
     //printf("pp centralityF_temp %d\n",centralityF_temp);
     
     if( centralityFtemp <=  0) fCentralityF =  0;
     else if ( centralityFtemp >  0 && centralityFtemp <  2) fCentralityF = 1;
     else if ( centralityFtemp >=  2 && centralityFtemp <  3) fCentralityF = 2;
     else if ( centralityFtemp >= 3 && centralityFtemp <  4) fCentralityF = 3;
     else if ( centralityFtemp >= 4 && centralityFtemp <  5) fCentralityF = 4;
     else if ( centralityFtemp >= 5 && centralityFtemp <  10) fCentralityF = 5;
     else if ( centralityFtemp >= 10 && centralityFtemp <  20) fCentralityF = 6;
     else if ( centralityFtemp >= 20 && centralityFtemp <  30) fCentralityF = 7;
     else if ( centralityFtemp >= 30 && centralityFtemp <  40) fCentralityF = 8;
     else if ( centralityFtemp >= 40 && centralityFtemp <  50) fCentralityF = 9;
     else if ( centralityFtemp >= 50) fCentralityF = 10;
    
     return kTRUE; 
     
   }

   return kFALSE;

  //printf("centrality %f\n",fCentralityF);
   
 }

 //printf("centrality %f\n",fCentralityF);

}
//___________________________________________________
void AliAnalysisTaskHFE::RejectionPileUpVertexRangeEventCut() {
  //
  // Recover the centrality of the event from ESD or AOD
  //
 if(IsAODanalysis()){

   AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
   if(!fAOD){
     AliError("AOD Event required for AOD Analysis")
       return;
   }
   // PileUp
   if(fRemovePileUp && fAOD->IsPileupFromSPD()) fIdentifiedAsPileUp = kTRUE; 
   // Z vertex
   if(TMath::Abs(fAOD->GetPrimaryVertex()->GetZ()) > fCuts->GetVertexRange()) fIdentifiedAsOutInz = kTRUE;
   // Event Cut
   fPassTheEventCut = kTRUE;
   if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fAOD)) fPassTheEventCut = kFALSE; 
   
   
 } else {
   
   AliDebug(3, "Processing ESD Centrality");
   AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
   if(!fESD){
     AliError("ESD Event required for ESD Analysis")
       return;
   }
   // PileUp
   fIdentifiedAsPileUp = kFALSE;
   if(fRemovePileUp && fESD->IsPileupFromSPD()) fIdentifiedAsPileUp = kTRUE; 
   // Z vertex
   fIdentifiedAsOutInz = kFALSE;
   if(fESD->GetPrimaryVertexTracks()){
       if(TMath::Abs(fESD->GetPrimaryVertexTracks()->GetZ()) > fCuts->GetVertexRange()) fIdentifiedAsOutInz = kTRUE;
   }
   //Event Cut
   fPassTheEventCut = kTRUE;
   if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) fPassTheEventCut = kFALSE;   
  

 }

}

