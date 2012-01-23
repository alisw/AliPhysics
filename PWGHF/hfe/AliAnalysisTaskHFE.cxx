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
#include <TBits.h>
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
#include <TObjString.h>
#include <TParticle.h>
#include <TProfile.h>
#include <TString.h>
#include <TF1.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliOADBContainer.h"
#include "AliStack.h"
#include "AliTriggerAnalysis.h"
#include "AliVVertex.h"
#include "TTreeStream.h"
#include "AliESDtrackCuts.h"

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
  , fFillNoCuts(kFALSE)
  , fBackGroundFactorApply(kFALSE)
  , fRemovePileUp(kFALSE)
  , fIdentifiedAsPileUp(kFALSE)
  , fIdentifiedAsOutInz(kFALSE)
  , fPassTheEventCut(kFALSE)
  , fRejectKinkMother(kTRUE)
  , fisppMultiBin(kFALSE)
  , fisNonHFEsystematics(kFALSE)
  , fSpecialTrigger(NULL)
  , fCentralityF(-1)
  , fContributors(0.5)
  , fWeightBackGround(0.)
  , fVz(0.0)
  , fHadronBackgroundOADB(NULL)
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
  , fExtraCuts(NULL)
  , fQA(NULL)
  , fOutput(NULL)
  , fHistMCQA(NULL)
  , fHistSECVTX(NULL)
  , fHistELECBACKGROUND(NULL)
  , fQACollection(NULL)
  , fDebugLevel(0)
  , fTreeStream(NULL)
{
  //
  // Dummy constructor
  //
  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fkBackGroundFactorArray, 0, sizeof(TF1 *) * 12);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
  memset(&fisppMultiBin, kFALSE, sizeof(fisppMultiBin));
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const char * name):
  AliAnalysisTaskSE(name)
  , fQAlevel(0)
  , fPlugins(0)
  , fFillSignalOnly(kTRUE)
  , fFillNoCuts(kFALSE)
  , fBackGroundFactorApply(kFALSE)
  , fRemovePileUp(kFALSE)
  , fIdentifiedAsPileUp(kFALSE)
  , fIdentifiedAsOutInz(kFALSE)
  , fPassTheEventCut(kFALSE)  
  , fRejectKinkMother(kTRUE)
  , fisppMultiBin(kFALSE)
  , fisNonHFEsystematics(kFALSE)
  , fSpecialTrigger(NULL)
  , fCentralityF(-1)
  , fContributors(0.5)
  , fWeightBackGround(0.)
  , fVz(0.0)
  , fHadronBackgroundOADB(NULL)
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
  , fExtraCuts(NULL)
  , fQA(NULL)
  , fOutput(NULL)
  , fHistMCQA(NULL)
  , fHistSECVTX(NULL)
  , fHistELECBACKGROUND(NULL)
  , fQACollection(0x0)
  , fDebugLevel(0)
  , fTreeStream(NULL)
{
  //
  // Default constructor
  // 
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  fPID = new AliHFEpid("hfePid");
  fPIDqa = new AliHFEpidQAmanager;
  fVarManager = new AliHFEvarManager("hfeVarManager");

  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fkBackGroundFactorArray, 0, sizeof(TF1 *) * 12);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
  memset(&fisppMultiBin, kFALSE, sizeof(fisppMultiBin));
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref):
  AliAnalysisTaskSE(ref)
  , fQAlevel(0)
  , fPlugins(0)
  , fFillSignalOnly(ref.fFillSignalOnly)
  , fFillNoCuts(ref.fFillNoCuts)
  , fBackGroundFactorApply(ref.fBackGroundFactorApply)
  , fRemovePileUp(ref.fRemovePileUp)
  , fIdentifiedAsPileUp(ref.fIdentifiedAsPileUp)
  , fIdentifiedAsOutInz(ref.fIdentifiedAsOutInz)
  , fPassTheEventCut(ref.fPassTheEventCut)
  , fRejectKinkMother(ref.fRejectKinkMother)
  , fisppMultiBin(ref.fisppMultiBin)
  , fisNonHFEsystematics(ref.fisNonHFEsystematics)
  , fSpecialTrigger(ref.fSpecialTrigger)
  , fCentralityF(ref.fCentralityF)
  , fContributors(ref.fContributors)
  , fWeightBackGround(ref.fWeightBackGround)
  , fVz(ref.fVz)
  , fHadronBackgroundOADB(ref.fHadronBackgroundOADB)
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
  , fExtraCuts(NULL)
  , fQA(NULL)
  , fOutput(NULL)
  , fHistMCQA(NULL)
  , fHistSECVTX(NULL)
  , fHistELECBACKGROUND(NULL)
  , fQACollection(NULL)
  , fDebugLevel(ref.fDebugLevel)
  , fTreeStream(NULL)
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
  target.fFillNoCuts = fFillNoCuts;
  target.fBackGroundFactorApply = fBackGroundFactorApply;
  target.fRemovePileUp = fRemovePileUp;
  target.fIdentifiedAsPileUp = fIdentifiedAsPileUp;
  target.fIdentifiedAsOutInz = fIdentifiedAsOutInz;
  target.fPassTheEventCut = fPassTheEventCut;
  target.fRejectKinkMother = fRejectKinkMother;
  target.fisppMultiBin =   fisppMultiBin;
  target.fisNonHFEsystematics = fisNonHFEsystematics;
  target.fSpecialTrigger = fSpecialTrigger;
  target.fCentralityF = fCentralityF;
  target.fContributors = fContributors;
  target.fWeightBackGround = fWeightBackGround;
  target.fVz = fVz;
  target.fHadronBackgroundOADB = fHadronBackgroundOADB;
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
  target.fExtraCuts = fExtraCuts;
  target.fQA = fQA;
  target.fOutput = fOutput;
  target.fHistMCQA = fHistMCQA;
  target.fHistSECVTX = fHistSECVTX;
  target.fHistELECBACKGROUND = fHistELECBACKGROUND;
  target.fQACollection = fQACollection;
  target.fDebugLevel = fDebugLevel;
  target.fTreeStream = fTreeStream;
}

//____________________________________________________________
AliAnalysisTaskHFE::~AliAnalysisTaskHFE(){
  //
  // Destructor
  //
  if(fPID) delete fPID;
  if(fPIDpreselect) delete fPIDpreselect;
  if(fVarManager) delete fVarManager;
  if(fCFM) delete fCFM;
  if(fTriggerAnalysis) delete fTriggerAnalysis;
  if(fSignalCuts) delete fSignalCuts;
  if(fSecVtx) delete fSecVtx;
  if(fMCQA) delete fMCQA;
  if(fElecBackGround) delete fElecBackGround;
  if(fSpecialTrigger) delete fSpecialTrigger;
  // Delete output objects only if we are not running in PROOF mode because otherwise this produces a crash during merging
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(mgr && mgr->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
    if(fPIDqa) delete fPIDqa;
    if(fOutput) delete fOutput;
    if(fQA) delete fQA;
    if(fTreeStream) delete fTreeStream;
  }
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
  fQACollection->CreateTH1F("nElectron", "Number of electrons", 100, 0, 100);
  fQACollection->CreateTH2F("radius", "Production Vertex", 100, 0.0, 5.0, 100, 0.0, 5.0);
 
  InitPIDperformanceQA();
  InitContaminationQA();
  InitHistoITScluster();
  fQA->Add(fQACollection);

  // Initialize PID
  fPID->SetHasMCData(HasMCData());
  if(!fPID->GetNumberOfPIDdetectors()) fPID->AddDetector("TPC", 0);
  if(IsQAOn(kPIDqa)){
    AliInfo("PID QA switched on");
    fPIDqa->Initialize(fPID);
    fQA->Add(fPIDqa->MakeList("HFEpidQA"));
  }
  fPID->SortDetectors();

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
    if(IsPbPb()) fMCQA->SetPbPb();
    if(fisppMultiBin) fMCQA->SetPPMultiBin();
    fMCQA->CreatDefaultHistograms(fHistMCQA);
    if(!fFillSignalOnly) fMCQA->SetBackgroundWeightFactor(fElecBackgroundFactor[0][0][0],fBinLimit);
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
      Int_t nbins = variable->GetNumberOfBins();
      if(variable->HasUserDefinedBinning()){
        varManager->AddVariable(namee, nbins, variable->GetBinning());
      } else {
        varManager->AddVariable(namee, nbins, variable->GetMinimum(), variable->GetMaximum());
      }
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

  Bool_t isProof = AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() == AliAnalysisManager::kProofAnalysis;
  if(fDebugLevel && !isProof){
    AliDebug(1,"Create OutputStream");
    fTreeStream = new TTreeSRedirector(Form("HFEdebugTree%s.root", GetName()));
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
  if(!fPID->IsInitialized()){
    // Initialize PID with the given run number
    fPID->InitializePID(fInputEvent->GetRunNumber());
  }

  // Initialize hadronic background from OADB Container
  AliDebug(2, Form("Apply background factors: %s, OADB Container %p", fBackGroundFactorApply ? "Yes" : "No", fHadronBackgroundOADB));
  if(fBackGroundFactorApply && !TestBit(kBackgroundInitialized)){ 
    AliDebug(2, "Initializing Background from OADB");
    if(!InitializeHadronBackground(fInputEvent->GetRunNumber())) AliError("Failed initializing hadronic background parameterization from OADB");
    else AliDebug(2, "Successfully loaded Background from OADB");
    SetBit(kBackgroundInitialized); 
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
  fCentralityF = -1;
  if(!ReadCentrality()) fCentralityF = -1;
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

  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
    pidResponse = AliHFEtools::GetDefaultPID(HasMCData(), fInputEvent->IsA() == AliAODEvent::Class());
  }
  fPID->SetPIDResponse(pidResponse);
  if(fPIDpreselect) fPIDpreselect->SetPIDResponse(pidResponse);

  // Event loop
  if(IsAODanalysis()){
    ProcessAOD();
  } else {
    const char *specialTrigger = GetSpecialTrigger(fInputEvent->GetRunNumber());
    // Check Trigger selection
    if(specialTrigger){
      AliDebug(2, Form("Special Trigger requested: %s", specialTrigger));
      AliESDEvent *ev = dynamic_cast<AliESDEvent *>(fInputEvent);
      if(!(ev && ev->IsTriggerClassFired(specialTrigger))){
        AliDebug(2, "Event not selected"); 
        return;
      } else AliDebug(2, "Event Selected");
    } else AliDebug(2, "No Special Trigger requested");

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
   if(!((fIdentifiedAsPileUp) || (TMath::Abs(fVz) > fCuts->GetVertexRange()) || (fCentralityF < 0))){ //kStepMCGeneratedZOutNoPileUpCentralityFine
    if (HasMCData() && IsQAOn(kMCqa)) {
      AliDebug(2, "Running MC QA");

      if(fMCEvent->Stack()){
        fMCQA->SetMCEvent(fMCEvent);
        fMCQA->SetGenEventHeader(fMCEvent->GenEventHeader());
        fMCQA->SetCentrality(fCentralityF);
        if(IsPbPb()) { fMCQA->SetPbPb();}
        else
        {
            if(fisppMultiBin) fMCQA->SetPPMultiBin();
            else fMCQA->SetPP();
        }
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
        //fMCQA->EndOfEventAna(AliHFEmcQA::kCharm);
        //fMCQA->EndOfEventAna(AliHFEmcQA::kBeauty);
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
    AliError("ESD Event required for ESD Analysis");
    return;
  }

  // Set magnetic field if V0 task on
  if(fTaggedTrackAnalysis) {
    fTaggedTrackAnalysis->SetMagneticField(fESD->GetMagneticField());
    fTaggedTrackAnalysis->SetCentrality(fCentralityF);
    if(IsPbPb()) fTaggedTrackAnalysis->SetPbPb();
    else fTaggedTrackAnalysis->SetPP();
  }

  // Do event Normalization
  Double_t eventContainer[4];
  eventContainer[0] = 0.; 
  if(HasMCData()) eventContainer[0] = fVz;
  else {
    if(fESD->GetPrimaryVertexSPD()) eventContainer[0] = fESD->GetPrimaryVertexSPD()->GetZ();
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
  if(TMath::Abs(fCentralityF) < 0) return; 
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
  AliMCParticle *mctrackmother = NULL;
  Int_t pid = 0;

  Bool_t signal = kTRUE;

  fCFM->SetRecEventInfo(fESD);

  // minjung for IP QA(temporary ~ 2weeks)
  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfetmpCuts","HFE tmp Cuts");
  }   
  fExtraCuts->SetRecEventInfo(fESD);

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

      if(fFillSignalOnly && !fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE; 
      else AliDebug(3, "Signal Electron");

      // Fill K pt for Ke3 contributions
      if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode())==321)) fQACollection->Fill("Kptspectra",mctrack->Pt());
      else if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode())==130)) fQACollection->Fill("K0Lptspectra",mctrack->Pt());
    } 
    // Cache new Track information inside the var manager
    fVarManager->NewTrack(track, mctrack, fCentralityF, -1, signal);

    if(fFillNoCuts) {
      if(signal || !fFillSignalOnly){
        fVarManager->FillContainer(fContainer, "recTrackContReco", AliHFEcuts::kStepRecNoCut, kFALSE);
        fVarManager->FillContainer(fContainer, "recTrackContMC", AliHFEcuts::kStepRecNoCut, kTRUE);
      }
    }

    // RecKine: ITSTPC cuts  
    if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    
    
    // RecPrim
    if(fRejectKinkMother) { 
      if(track->GetKinkIndex(0) != 0) continue; } // Quick and dirty fix to reject both kink mothers and daughters
    if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;

    if(fTreeStream && fDebugLevel >= 2){
      // Debug streaming of PID-related quantities
      Double_t nSigmaTOF = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
      Double_t nSigmaTPC = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
      if(TMath::Abs(nSigmaTOF) < 5 && TMath::Abs(nSigmaTPC) < 5){
        // we are not interested in tracks which are more than 5 sigma away from the electron hypothesis in either TOF or TPC
        Double_t charge = track->Charge() > 0 ? 1. : -1.;
        Char_t myv0pid = v0pid;
        Double_t momentum = track->P() * charge;
        Double_t transversemomentum = track->Pt() * charge;
        Int_t run = fInputEvent->GetRunNumber();
        Double_t eta = track->Eta();
        Double_t phi = track->Phi();
        UChar_t ntracklets = track->GetTRDntrackletsPID();
        UChar_t nclustersTPCPID = track->GetTPCsignalN();
        UChar_t nclustersTPCshared = 0;
        const TBits &sharedTPC = track->GetTPCSharedMap();
        for(Int_t ibit = 0; ibit < 160; ibit++) if(sharedTPC.TestBitNumber(ibit)) nclustersTPCshared++;
        UChar_t nclustersTPC = track->GetTPCncls();
        UChar_t nclustersITS = track->GetITSclusters(NULL);
        UChar_t nclustersTRD = track->GetTRDncls();
        UChar_t hasClusterITS[6], hasTrackletTRD[6];
        UChar_t itsPixel = track->GetITSClusterMap();
        for(Int_t icl = 0; icl < 6; icl++) hasClusterITS[icl] = TESTBIT(itsPixel, icl) ? 1 : 0;
        for(Int_t itl = 0; itl < 6; itl++){
          Int_t nSliceNonZero = 0;
          for(Int_t islice = 0; islice < 8; islice++){
            if(track->GetTRDslice(itl, islice) > 0.001) nSliceNonZero++;
          }
          hasTrackletTRD[itl] = nSliceNonZero ? 1 : 0;
        }
        Double_t pidprobs[5];
        track->GetTRDpid(pidprobs);
        Double_t likeEleTRD = pidprobs[0];
        Double_t likeEleTRDn = likeEleTRD/(likeEleTRD + pidprobs[2]);
        (*fTreeStream) << "PIDdebug"
              << "signal="              << signal
              << "v0pid="               << myv0pid
              << "run="                 << run
              << "p="                   << momentum
              << "pt="                  << transversemomentum
              << "eta="                 << eta
              << "phi="                 << phi
              << "ntracklets="          << ntracklets
              << "nclustersTPC="        << nclustersTPC
              << "nclustersTPCPID="     << nclustersTPCPID
              << "nclustersTPCshared="  << nclustersTPCshared
              << "nclustersITS="        << nclustersITS
              << "nclusters="           << nclustersTRD
              << "its0="                << hasClusterITS[0]
              << "its1="                << hasClusterITS[1]
              << "its2="                << hasClusterITS[2]
              << "its3="                << hasClusterITS[3]
              << "its4="                << hasClusterITS[4]
              << "its5="                << hasClusterITS[5]
              << "trd0="                << hasTrackletTRD[0]
              << "trd1="                << hasTrackletTRD[1]
              << "trd2="                << hasTrackletTRD[2]
              << "trd3="                << hasTrackletTRD[3]
              << "trd4="                << hasTrackletTRD[4]
              << "trd5="                << hasTrackletTRD[5]
              << "TOFsigmaEl="          << nSigmaTOF
              << "TPCsigmaEl="          << nSigmaTPC
              << "TRDlikeEl="           << likeEleTRD
              << "TRDlikeEln="          << likeEleTRDn
              << "\n";
      }
    }

    // HFEcuts: ITS layers cuts
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
  
    // HFE cuts: TOF PID and mismatch flag
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTOF, track)) continue;

    // HFE cuts: TPC PID cleanup
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;

    // HFEcuts: Nb of tracklets TRD0
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTRD, track)) continue;

    // Fill correlation maps before PID
    if(signal && fContainer->GetCorrelationMatrix("correlationstepbeforePID")) {
      //printf("Fill correlation maps before PID\n");
      fVarManager->FillCorrelationMatrix(fContainer->GetCorrelationMatrix("correlationstepbeforePID"));
    }

    if(HasMCData()){
      FillProductionVertex(track);

      if(fMCQA){
        fMCQA->SetCentrality(fCentralityF);
        if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) == 11)){
         Double_t weightElecBgV0[kBgLevels] = {0.,0.,0.};
         for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
           weightElecBgV0[iLevel] = fMCQA->GetWeightFactor(mctrack, iLevel); // positive:conversion e, negative: nonHFE 
           if(!fisNonHFEsystematics)break;   
         }
         
         if(fisNonHFEsystematics){
           //Fill additional containers for electron source distinction
           Int_t elecSource = 0;
           elecSource = fMCQA->GetElecSource(mctrack->Particle());
           const Char_t *sourceName[kElecBgSpecies]={"Pion","Eta","Omega","Phi","EtaPrime","Rho"};
           const Char_t *levelName[kBgLevels]={"Best","Lower","Upper"};
           Int_t iName = 0;
           for(Int_t iSource = AliHFEmcQA::kPi0; iSource <=  AliHFEmcQA::kGammaRho0; iSource++){
             if((iSource == AliHFEmcQA::kElse)||(iSource == AliHFEmcQA::kMisID)) continue;
             if(elecSource == iSource){
               for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
                 if(weightElecBgV0[iLevel]>0){ fVarManager->FillContainer(fContainer, Form("conversionElecs%s%s",sourceName[iName], levelName[iLevel]), 3, kFALSE, weightElecBgV0[iLevel]);}
                 else if(weightElecBgV0[iLevel]<0){ fVarManager->FillContainer(fContainer, Form("mesonElecs%s%s",sourceName[iName], levelName[iLevel]), 3, kFALSE, -1*weightElecBgV0[iLevel]);}
               }
               break;
             }
             iName++;
             if(iName == kElecBgSpecies)iName = 0;
           }
         }
         //else{
           if(weightElecBgV0[0]>0) fVarManager->FillContainer(fContainer, "conversionElecs", 3, kFALSE, weightElecBgV0[0]);
           else if(weightElecBgV0[0]<0) fVarManager->FillContainer(fContainer, "mesonElecs", 3, kFALSE, -1*weightElecBgV0[0]);
           //}
        }
      }
    }
    
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    if(HasMCData()) hfetrack.SetMCTrack(mctrack);
    hfetrack.SetCentrality(fCentralityF);
    if(IsPbPb()) hfetrack.SetPbPb();
    else hfetrack.SetPP();
    fPID->SetVarManager(fVarManager);
    if(!fPID->IsSelected(&hfetrack, fContainer, "recTrackCont", fPIDqa)) continue;
    nElectronCandidates++;

    // Temporary histogram for chi2/ITS cluster
    if(IsPbPb()) {
            TBits shared = track->GetTPCSharedMap();
          Int_t sharebit=0;
            if(shared.CountBits() >= 2) sharebit=1;

            Double_t itsChi2[7] = {track->Pt(),track->Eta(), track->Phi(),
            fCentralityF,track->GetTPCsignalN(), sharebit,
            track->GetITSchi2()/static_cast<Double_t>(track->GetNcls(0))};
            fQACollection->Fill("fChi2perITScluster", itsChi2);
    }
    else{
            Double_t itsChi2[3] = {track->Pt(), fCentralityF, track->GetITSchi2()/static_cast<Double_t>(track->GetNcls(0))};
            fQACollection->Fill("fChi2perITScluster", itsChi2);
    }

    // Fill Histogram for Hadronic Background
    if(HasMCData()){
      if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) != 11))
        fVarManager->FillContainer(fContainer, "hadronicBackground", UInt_t(0), kFALSE);
      else if(mctrack){
        // Fill Ke3 contributions
        Int_t glabel=TMath::Abs(mctrack->GetMother());
        if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
          if(TMath::Abs(mctrackmother->Particle()->GetPdgCode())==321)
            fQACollection->Fill("Ke3Kecorr",mctrackmother->Pt(),mctrack->Pt());
          else if(TMath::Abs(mctrackmother->Particle()->GetPdgCode())==130)
            fQACollection->Fill("Ke3K0Lecorr",mctrackmother->Pt(),mctrack->Pt());
        }
      }
    }

    // Fill Containers
    if(signal) {
      // Apply weight for background contamination
            if(fBackGroundFactorApply) {
              if(IsPbPb() && fCentralityF >= 0) fWeightBackGround =  fkBackGroundFactorArray[fCentralityF >= 0 ? fCentralityF : 0]->Eval(TMath::Abs(track->P()));
              else    fWeightBackGround =  fkBackGroundFactorArray[0]->Eval(TMath::Abs(track->P())); // pp case

              if(fWeightBackGround < 0.0) fWeightBackGround = 0.0;
              else if(fWeightBackGround > 1.0) fWeightBackGround = 1.0;
        // weightBackGround as special weight
        fVarManager->FillContainer(fContainer, "hadronicBackground", 1, kFALSE, fWeightBackGround);
      }
      fVarManager->FillCorrelationMatrix(fContainer->GetCorrelationMatrix("correlationstepafterPID"));
    }

    Bool_t bTagged=kFALSE;
    if(GetPlugin(kSecVtx)) {
      AliDebug(2, "Running Secondary Vertex Analysis");
      if(fSecVtx->Process(track) && signal) {
        fVarManager->FillContainer(fContainer, "recTrackContSecvtxReco", AliHFEcuts::kStepHFEcutsSecvtx, kFALSE);
        fVarManager->FillContainer(fContainer, "recTrackContSecvtxMC", AliHFEcuts::kStepHFEcutsSecvtx, kTRUE);
        bTagged=kTRUE;
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


    // high dca track study [for temporary] ---------
    if (fTreeStream && fDebugLevel >= 1){
      Float_t b[2] = {0.,0.};
      Float_t bCov[3] = {0.,0.,0.};
      track->GetImpactParameters(b,bCov);
      Double_t dataD[11];
      dataD[0] = b[0]; // impact parameter xy
      dataD[1] = b[1]; // impact parameter z
      dataD[2] = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]); // impact parameter space
      dataD[3]=0; dataD[4]=0;
      if(bCov[0]>0) dataD[3] = b[0]/TMath::Sqrt(bCov[0]); // normalised impact parameter xy
      if(bCov[2]>0) dataD[4] = b[1]/TMath::Sqrt(bCov[2]); // normalised impact parameter z
      dataD[5] = AliESDtrackCuts::GetSigmaToVertex(track); // n_sigma
      dataD[6] = track->GetTPCclusters(0x0);
      dataD[7] = track->GetITSclusters(0x0);
      dataD[8] = track->Eta();
      dataD[9] = track->Pt();
      Double_t p = track->GetP();
      Double_t phi = track->Phi();
      if(HasMCData()){
        if(fSignalCuts->IsCharmElectron(track)) dataD[10] = 0;
        else if(fSignalCuts->IsBeautyElectron(track)) dataD[10] = 1;
        else if(fSignalCuts->IsGammaElectron(track)) dataD[10] = 2;
        else if(fSignalCuts->IsNonHFElectron(track)) dataD[10] = 3;
        else if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) != 11)) dataD[10] = 4;
        else dataD[10] = 5;
      } 
      Double_t vtx[3];
      fInputEvent->GetPrimaryVertex()->GetXYZ(vtx);
      Double_t nt = fInputEvent->GetPrimaryVertex()->GetNContributors();
      Double_t nSigmaTPC = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
      Double_t runn = (Double_t)fInputEvent->GetRunNumber();

      //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",dataD[0],dataD[1],dataD[2],dataD[3],dataD[4],dataD[5],dataD[6],dataD[7],dataD[8],dataD[9]);
      (*fTreeStream)<<"dcaQA"<<
        "dcaR="<<dataD[0]<<
        "dcaZ="<<dataD[1]<<
        "dca="<<dataD[2]<<
        "dcaSR="<<dataD[3]<<
        "dcaSZ="<<dataD[4]<<
        "dcaS="<<dataD[5]<<
        "nTPCclus="<<dataD[6]<<
        "nITSclus="<<dataD[7]<<
        "eta="<<dataD[8]<<
        "pt="<<dataD[9]<<
              "p="<< p <<
              "phi="<< phi <<
        "source="<<dataD[10] <<
        "vx=" << vtx[0] <<
        "vy=" << vtx[1] <<
        "vz=" << vtx[2] << 
              "nt=" << nt <<
        "TPCnSigma=" << nSigmaTPC <<
              "run=" << runn
        << "\n";
    }
    //-------------------------------------------

    if (GetPlugin(kDEstep)) { 
      Double_t weightElecBgV0[kBgLevels] = {0.,0.,0.,};
      Int_t elecSource = 0;
      // minjung for IP QA(temporary ~ 2weeks)
      Double_t hfeimpactR=0., hfeimpactnsigmaR=0.;
      fExtraCuts->GetHFEImpactParameters(track, hfeimpactR, hfeimpactnsigmaR);
      fQACollection->Fill("dataDca",track->Pt(),hfeimpactR);
      fQACollection->Fill("dataDcaSig",track->Pt(),hfeimpactnsigmaR);
      fQACollection->Fill("dataDcaSigDca",hfeimpactR,hfeimpactnsigmaR);
      if(HasMCData()){
        // minjung for IP QA(temporary ~ 2weeks)
        if(fSignalCuts->IsCharmElectron(track)){
          fQACollection->Fill("charmDca",track->Pt(),hfeimpactR);
          fQACollection->Fill("charmDcaSig",track->Pt(),hfeimpactnsigmaR);
        }
        else if(fSignalCuts->IsBeautyElectron(track)){
          fQACollection->Fill("beautyDca",track->Pt(),hfeimpactR);
          fQACollection->Fill("beautyDcaSig",track->Pt(),hfeimpactnsigmaR);
        }
        else if(fSignalCuts->IsGammaElectron(track)){
          fQACollection->Fill("conversionDca",track->Pt(),hfeimpactR);
          fQACollection->Fill("conversionDcaSig",track->Pt(),hfeimpactnsigmaR);
          fQACollection->Fill("conversionDcaSigDca",hfeimpactR,hfeimpactnsigmaR);
        }
        else if(fSignalCuts->IsNonHFElectron(track)){
          fQACollection->Fill("nonhfeDca",track->Pt(),hfeimpactR);
          fQACollection->Fill("nonhfeDcaSig",track->Pt(),hfeimpactnsigmaR);
        }

        if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) != 11)){
          fQACollection->Fill("hadronsBeforeIPcut",track->Pt());
          fQACollection->Fill("hadronsBeforeIPcutMC",mctrack->Pt());
        }
        if(fMCQA && !fFillSignalOnly) {
          
          for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
            weightElecBgV0[iLevel] = fMCQA->GetWeightFactor(mctrack, iLevel); // positive:conversion e, negative: nonHFE 
            if(!fisNonHFEsystematics)break;        
          }
          
          if(fisNonHFEsystematics){
            //Fill additional containers for electron source distinction           
            elecSource = fMCQA->GetElecSource(mctrack->Particle());
            const Char_t *sourceName[kElecBgSpecies]={"Pion","Eta","Omega","Phi","EtaPrime","Rho"};
            const Char_t *levelName[kBgLevels]={"Best","Lower","Upper"};
            Int_t iName = 0;
            for(Int_t iSource = AliHFEmcQA::kPi0; iSource <=  AliHFEmcQA::kGammaRho0; iSource++){
              if((iSource == AliHFEmcQA::kElse)||(iSource == AliHFEmcQA::kMisID)) continue;
              if(elecSource == iSource){
                for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
                  if(weightElecBgV0[iLevel]>0) fVarManager->FillContainer(fContainer, Form("conversionElecs%s%s",sourceName[iName], levelName[iLevel]), 0, kFALSE, weightElecBgV0[iLevel]);
                  else if(weightElecBgV0[iLevel]<0) fVarManager->FillContainer(fContainer, Form("mesonElecs%s%s",sourceName[iName], levelName[iLevel]), 0, kFALSE, -1*weightElecBgV0[iLevel]);
                }
                break;
              }
              iName++;
              if(iName == kElecBgSpecies)iName = 0;
            }
          }
          //else{
            if(weightElecBgV0[0]>0) fVarManager->FillContainer(fContainer, "conversionElecs", 0, kFALSE, weightElecBgV0[0]);
            else if(weightElecBgV0[0]<0) fVarManager->FillContainer(fContainer, "mesonElecs", 0, kFALSE, -1*weightElecBgV0[0]);
            //}
        }
      }
      // Fill Containers for impact parameter analysis
      if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsDca + AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack,track)) continue;
      if(HasMCData()){
        if(fMCQA && !fFillSignalOnly) {
          for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
            weightElecBgV0[iLevel] = fMCQA->GetWeightFactor(mctrack, iLevel); // positive:conversion e, negative: nonHFE 
            if(!fisNonHFEsystematics)break;        
          }       
          if(fisNonHFEsystematics){
            //Fill additional containers for electron source distinction             
            elecSource = fMCQA->GetElecSource(mctrack->Particle());
            const Char_t *sourceName[kElecBgSpecies]={"Pion","Eta","Omega","Phi","EtaPrime","Rho"};
            const Char_t *levelName[kBgLevels]={"Best","Lower","Upper"};
            Int_t iName = 0;
            for(Int_t iSource = AliHFEmcQA::kPi0; iSource <=  AliHFEmcQA::kGammaRho0; iSource++){
              if((iSource == AliHFEmcQA::kElse)||(iSource == AliHFEmcQA::kMisID)) continue;
              if(elecSource == iSource){
                for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
                  if(weightElecBgV0[iLevel]>0) fVarManager->FillContainer(fContainer, Form("conversionElecs%s%s",sourceName[iName], levelName[iLevel]), 1, kFALSE, weightElecBgV0[iLevel]);
                  else if(weightElecBgV0[iLevel]<0) fVarManager->FillContainer(fContainer, Form("mesonElecs%s%s",sourceName[iName], levelName[iLevel]), 1, kFALSE, -1*weightElecBgV0[iLevel]);
                }
                break;
              }
              iName++;
              if(iName == kElecBgSpecies)iName = 0;
            }
          }
          // else{
            if(weightElecBgV0[0]>0) fVarManager->FillContainer(fContainer, "conversionElecs", 1, kFALSE, weightElecBgV0[0]);
            else if(weightElecBgV0[0]<0) fVarManager->FillContainer(fContainer, "mesonElecs", 1, kFALSE, -1*weightElecBgV0[0]);
            //}
        }
      }
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
    AliError("AOD Event required for AOD Analysis");
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
    track = fAOD->GetTrack(itrack); mctrack = NULL;
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
            if(IsPbPb()) weightBackGround =  fkBackGroundFactorArray[fCentralityF >= 0 ? fCentralityF : 0]->Eval(TMath::Abs(track->P()));
      else weightBackGround =  fkBackGroundFactorArray[0]->Eval(TMath::Abs(track->P()));

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
      AliDebug(1, Form("Candidate Selected, filling THnSparse, PID: %d\n", mctrack ? mctrack->GetPdgCode(): -1));
      if(mctrack && ((pid = TMath::Abs(mctrack->GetPdgCode())) == 11)){
      
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

  // Step GeneratedZOutNoPileUp
  if((fIdentifiedAsPileUp) || (TMath::Abs(fVz) > fCuts->GetVertexRange()) || (fCentralityF < 0)) return kFALSE;
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
  if(!fCutspreselect->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTOF, track)) {
    survived = kFALSE;
    //printf("Did not pass AliHFEcuts::kStepHFEcutsTOF\n");
  }
  //else printf("Pass AliHFEcuts::kStepHFEcutsTOF\n");
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
  // 4th bin: Number of contributors > 0
  //
  
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

  if(HasMCData()){
    fContainer->CreateContainer("conversionElecs", "Container for weighted conversion electrons",4);
    fContainer->CreateContainer("mesonElecs", "Container for weighted electrons from meson decays",4);
   
    if(fisNonHFEsystematics){
      const Char_t *sourceName[kElecBgSpecies]={"Pion","Eta","Omega","Phi","EtaPrime","Rho"};
      const Char_t *levelName[kBgLevels]={"Best","Lower","Upper"};
      for(Int_t iSource = 0; iSource < kElecBgSpecies; iSource++){
        for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
          fContainer->CreateContainer(Form("conversionElecs%s%s",sourceName[iSource],levelName[iLevel]), Form("Container for weighted conversion electrons from %s grandm., %s level",sourceName[iSource],levelName[iLevel]),4);
          fContainer->CreateContainer(Form("mesonElecs%s%s",sourceName[iSource],levelName[iLevel]), Form("Container for weighted electrons from %s decays, %s level",sourceName[iSource],levelName[iLevel]),4);
        }
      }
    }
    //fContainer->CreateContainer("charmElecs", "Container for weighted charm electrons",2);
  }

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
  fQACollection->CreateTH1F("hadronsBeforeIPcutMC", "Hadrons before IP cut; MC p_{t}", iBin[0], kPtbound[0], kPtbound[1], 1);
  fQACollection->CreateTH1F("hadronsAfterIPcutMC", "Hadrons after IP cut; MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);

  fQACollection->CreateTH2F("Ke3Kecorr", "Ke3 decay e and K correlation; Ke3K p_{t}; Ke3e p_{t}; ",20,0.,20.,iBin[0],kPtbound[0],kPtbound[1], 1);
  fQACollection->CreateTH2F("Ke3K0Lecorr", "Ke3 decay e and K0L correlation; Ke3K0L p_{t}; Ke3e p_{t}; ",20,0.,20.,iBin[0],kPtbound[0],kPtbound[1], 1);
  fQACollection->CreateTH1F("Kptspectra", "Charged Kaons: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
  fQACollection->CreateTH1F("K0Lptspectra", "K0L: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);

  const Double_t kDCAbound[2] = {-5., 5.};
  const Double_t kDCAsigbound[2] = {-50., 50.};

  fQACollection->CreateTH2F("dataDcaSig", "data dca significance: dca sig ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAsigbound[0], kDCAsigbound[1], 0);
  fQACollection->CreateTH2F("charmDcaSig", "charm dca significance: dca sig ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAsigbound[0], kDCAsigbound[1], 0);
  fQACollection->CreateTH2F("beautyDcaSig", "beauty dca significance: dca sig ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAsigbound[0], kDCAsigbound[1], 0);
  fQACollection->CreateTH2F("conversionDcaSig", "conversion dca significance: dca sig ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAsigbound[0], kDCAsigbound[1], 0);
  fQACollection->CreateTH2F("nonhfeDcaSig", "nonhfe dca significance: dca sig ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAsigbound[0], kDCAsigbound[1], 0);

  fQACollection->CreateTH2F("dataDca", "data dca : dca ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAbound[0], kDCAbound[1], 0);
  fQACollection->CreateTH2F("charmDca", "charm dca : dca ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAbound[0], kDCAbound[1], 0);
  fQACollection->CreateTH2F("beautyDca", "beauty dca : dca ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAbound[0], kDCAbound[1], 0);
  fQACollection->CreateTH2F("conversionDca", "conversion dca : dca ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAbound[0], kDCAbound[1], 0);
  fQACollection->CreateTH2F("nonhfeDca", "nonhfe dca : dca ",iBin[0], kPtbound[0], kPtbound[1], 2000,kDCAbound[0], kDCAbound[1], 0);

  fQACollection->CreateTH2F("dataDcaSigDca", "data dca significance and dca correlation; dca; dca sig ", 2000, kDCAbound[0], kDCAbound[1], 2000, kDCAsigbound[0], kDCAsigbound[1], 0);
  fQACollection->CreateTH2F("conversionDcaSigDca", "conversion dca significance and dca correlation; dca; dca sig ", 2000, kDCAbound[0], kDCAbound[1], 2000, kDCAsigbound[0], kDCAsigbound[1], 0);
}

//____________________________________________________________
void AliAnalysisTaskHFE::InitHistoITScluster(){
  //
    // Initialize a temporary histogram to monitor the chi2/ITS cluster
    if(IsPbPb()) {
        const Int_t kNDim = 7;
        const Int_t kNBins[kNDim] = {88, 20,90,11, 160, 2, 1000};
        const Double_t kMin[kNDim] = {0.1, -1,0,  0.,0., 0,  0.};
        const Double_t kMax[kNDim] = {20., 1, 2.*TMath::Pi(), 11.,160, 2, 100.};
        fQACollection->CreateTHnSparse("fChi2perITScluster", "chi2/ITS cluster; p_{T} (GeV/c);eta;phi; centrality class;nclus;sharebit; #chi^{2}/ITS cluster", kNDim, kNBins, kMin, kMax);
        fQACollection->BinLogAxis("fChi2perITScluster", 0);
    }
    else
    {
        const Int_t kNDim = 3;
        const Int_t kNBins[kNDim] = {44, 11, 1000};
        const Double_t kMin[kNDim] = {0.1, 0., 0.};
        const Double_t kMax[kNDim] = {20., 11., 100.};
        fQACollection->CreateTHnSparse("fChi2perITScluster", "chi2/ITS cluster; p_{T} (GeV/c); centrality class; #chi^{2}/ITS cluster", kNDim, kNBins, kMin, kMax);
        fQACollection->BinLogAxis("fChi2perITScluster", 0);
    }
}

//____________________________________________________________
void AliAnalysisTaskHFE::SelectSpecialTrigger(const Char_t *trgclust, Int_t runMin, Int_t runMax){
  //
  // Select only events triggered by a special trigeer cluster
  //
  if(!fSpecialTrigger) fSpecialTrigger = new AliOADBContainer("SpecialTrigger");
  fSpecialTrigger->AppendObject(new TObjString(trgclust), runMin, runMax);
}

//____________________________________________________________
const Char_t * AliAnalysisTaskHFE::GetSpecialTrigger(Int_t run){
  //
  // Derive selected trigger string for given run
  //
  if(!fSpecialTrigger) return NULL;
  TObjString *trg = dynamic_cast<TObjString *>(fSpecialTrigger->GetObject(run));
  if(!trg) return NULL;
  return trg->String().Data();
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
  Int_t bin = -1;
  if(IsPbPb()) {
    // Centrality
    AliCentrality *centrality = fInputEvent->GetCentrality();
    Float_t fCentralityFtemp = centrality->GetCentralityPercentile("V0M");
    Float_t centralityLimits[12] = {0.,5.,10., 20., 30., 40., 50., 60.,70.,80., 90., 100.};
    for(Int_t ibin = 0; ibin < 11; ibin++){
      if(fCentralityFtemp >= centralityLimits[ibin] && fCentralityFtemp < centralityLimits[ibin+1]){
          bin = ibin;
        break;
      }
    } 
    if(bin == -1) bin = 11; // Overflow
  } else {
    // PP: Tracklet multiplicity, use common definition
    Int_t itsMultiplicity = GetITSMultiplicity(fInputEvent);
    Int_t multiplicityLimits[8] = {0, 1, 9, 17, 25, 36, 60, 500};
    for(Int_t ibin = 0; ibin < 7; ibin++){  
      if(itsMultiplicity >= multiplicityLimits[ibin] && itsMultiplicity < multiplicityLimits[ibin + 1]){
        bin = ibin;
        break;
      }
    }
    if(bin == -1) bin = 7;  // Overflow
  }
  fCentralityF = bin;
  AliDebug(2, Form("Centrality class %d\n", fCentralityF));
 
  // contributors, to be outsourced
  const AliVVertex *vtx;
  if(IsAODanalysis()){
    AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
    if(!fAOD){
      AliError("AOD Event required for AOD Analysis");
      return kFALSE;
    }
    vtx = fAOD->GetPrimaryVertex();
  } else {
    AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
    if(!fESD){
      AliError("ESD Event required for ESD Analysis");
      return kFALSE;
    }
    vtx = fESD->GetPrimaryVertexSPD();
  }
  if(!vtx){ 
    fContributors = 0.5;
    return kFALSE;
  }
  else {
    Int_t contributorstemp = vtx->GetNContributors();
    if( contributorstemp <=  0) fContributors =  0.5;
    else fContributors = 1.5;
  }
  return kTRUE;
}

//___________________________________________________
Int_t AliAnalysisTaskHFE::GetITSMultiplicity(AliVEvent *ev){
  //
  // Definition of the Multiplicity according to the JPSI group (F. Kramer)
  //
  Int_t nTracklets = 0;
  Int_t nAcc = 0;
  Double_t etaRange = 1.6;

  if (ev->IsA() == AliAODEvent::Class()) {
    AliAODTracklets *tracklets = ((AliAODEvent*)ev)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();
    for (Int_t nn = 0; nn < nTracklets; nn++) {
      Double_t theta = tracklets->GetTheta(nn);
      Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
      if (TMath::Abs(eta) < etaRange) nAcc++;
    }
  } else if (ev->IsA() == AliESDEvent::Class()) {
    nTracklets = ((AliESDEvent*)ev)->GetMultiplicity()->GetNumberOfTracklets();
    for (Int_t nn = 0; nn < nTracklets; nn++) {
       Double_t eta = ((AliESDEvent*)ev)->GetMultiplicity()->GetEta(nn);
      if (TMath::Abs(eta) < etaRange) nAcc++;
    }
  } else return -1;

  return nAcc;
}

//___________________________________________________
Bool_t AliAnalysisTaskHFE::InitializeHadronBackground(Int_t run){
  //
  // Initialize background factors array from the AliOADBContainer
  // The container is expected to provide a TObjArray with the name
  // "hadronBackground" and the size 12 for the given run number
  //
  AliDebug(1, "Deriving hadronic background parameterization from OADB container");
  TObjArray *functions = dynamic_cast<TObjArray *>(fHadronBackgroundOADB->GetObject(run, "hadronBackground"));
  if(!functions){
    AliDebug(1, "Content in the OADB Container is not a TObjArray");
    fBackGroundFactorApply = kFALSE;
    return kFALSE;
  } 
  if(functions->GetSize() < 12){
    AliDebug(1, Form("Size not matching: 12 expected, %d provided", functions->GetSize()));
    fBackGroundFactorApply = kFALSE;
    return kFALSE;
  }
  for(Int_t icent = 0; icent < 12; icent++) fkBackGroundFactorArray[icent] = dynamic_cast<const TF1 *>(functions->UncheckedAt(icent));
  return kTRUE;
}

//___________________________________________________
void AliAnalysisTaskHFE::RejectionPileUpVertexRangeEventCut() {
  //
  // Recover the centrality of the event from ESD or AOD
  //
 if(IsAODanalysis()){

   AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
   if(!fAOD){
     AliError("AOD Event required for AOD Analysis");
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
     AliError("ESD Event required for ESD Analysis");
       return;
   }
   // PileUp
   fIdentifiedAsPileUp = kFALSE;
   if(fRemovePileUp && fESD->IsPileupFromSPD()) fIdentifiedAsPileUp = kTRUE; 
   // Z vertex
   fIdentifiedAsOutInz = kFALSE;
   if (IsPbPb()) {
     //printf("PbPb\n");
     if(fESD->GetPrimaryVertex()){
       if(TMath::Abs(fESD->GetPrimaryVertex()->GetZ()) > fCuts->GetVertexRange()) fIdentifiedAsOutInz = kTRUE;
     }
   }
   else {
     //printf("pp\n");
     if(fESD->GetPrimaryVertexTracks()){
       if(TMath::Abs(fESD->GetPrimaryVertexTracks()->GetZ()) > fCuts->GetVertexRange()) fIdentifiedAsOutInz = kTRUE;
     }
   }
   //Event Cut
   fPassTheEventCut = kTRUE;
   if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) fPassTheEventCut = kFALSE;   
  

 }

}

