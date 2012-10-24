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

#include "AliESDtrackCuts.h"
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

#include "AliHFEcollection.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEelecbackground.h"
#include "AliHFENonPhotonicElectron.h"
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
#include "AliAODMCHeader.h"
#include "TClonesArray.h"

ClassImp(AliAnalysisTaskHFE)

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE():
AliAnalysisTaskSE("PID efficiency Analysis")
  , fAODMCHeader(NULL)
  , fAODArrayMCInfo(NULL)
  , fQAlevel(0)
  , fPlugins(0)
  , fFillSignalOnly(kTRUE)
  , fFillNoCuts(kFALSE)
  , fUseFilterAOD(kTRUE)
  , fApplyCutAOD(kFALSE)
  , fFilter(1<<4)
  , fBackGroundFactorApply(kFALSE)
  , fRemovePileUp(kFALSE)
  , fIdentifiedAsPileUp(kFALSE)
  , fIdentifiedAsOutInz(kFALSE)
  , fPassTheEventCut(kFALSE)
  , fRejectKinkMother(kTRUE)
  , fisppMultiBin(kFALSE)
  , fPbPbUserCentralityBinning(kFALSE)
  , fisNonHFEsystematics(kFALSE)
  , fSpecialTrigger(NULL)
  , fCentralityF(-1)
  , fCentralityPercent(-1)
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
  , fBackgroundSubtraction(NULL)
  , fQA(NULL)
  , fOutput(NULL)
  , fHistMCQA(NULL)
  , fHistSECVTX(NULL)
  , fHistELECBACKGROUND(NULL)
  , fQACollection(NULL)
  , fQAAODCollection(NULL)
{
  //
  // Dummy constructor
  //
  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fkBackGroundFactorArray, 0, sizeof(TF1 *) * 12);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
  memset(&fisppMultiBin, kFALSE, sizeof(fisppMultiBin));
  memset(fCentralityLimits, 0, sizeof(Float_t) * 12);


}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const char * name):
  AliAnalysisTaskSE(name)
  , fAODMCHeader(NULL)
  , fAODArrayMCInfo(NULL)
  , fQAlevel(0)
  , fPlugins(0)
  , fFillSignalOnly(kTRUE)
  , fFillNoCuts(kFALSE)
  , fUseFilterAOD(kTRUE)
  , fApplyCutAOD(kFALSE)
  , fFilter(1<<4)
  , fBackGroundFactorApply(kFALSE)
  , fRemovePileUp(kFALSE)
  , fIdentifiedAsPileUp(kFALSE)
  , fIdentifiedAsOutInz(kFALSE)
  , fPassTheEventCut(kFALSE)  
  , fRejectKinkMother(kTRUE)
  , fisppMultiBin(kFALSE)
  , fPbPbUserCentralityBinning(kFALSE)
  , fisNonHFEsystematics(kFALSE)
  , fSpecialTrigger(NULL)
  , fCentralityF(-1)
  , fCentralityPercent(-1)
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
  , fBackgroundSubtraction(NULL)
  , fQA(NULL)
  , fOutput(NULL)
  , fHistMCQA(NULL)
  , fHistSECVTX(NULL)
  , fHistELECBACKGROUND(NULL)
  , fQACollection(0x0)
  , fQAAODCollection(NULL)
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
  memset(fCentralityLimits, 0, sizeof(Float_t) * 12);

}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref):
  AliAnalysisTaskSE(ref)
  , fAODMCHeader(NULL)
  , fAODArrayMCInfo(NULL)
  , fQAlevel(0)
  , fPlugins(0)
  , fFillSignalOnly(ref.fFillSignalOnly)
  , fFillNoCuts(ref.fFillNoCuts)
  , fUseFilterAOD(ref.fUseFilterAOD)
  , fApplyCutAOD(ref.fApplyCutAOD)
  , fFilter(ref.fFilter) 
  , fBackGroundFactorApply(ref.fBackGroundFactorApply)
  , fRemovePileUp(ref.fRemovePileUp)
  , fIdentifiedAsPileUp(ref.fIdentifiedAsPileUp)
  , fIdentifiedAsOutInz(ref.fIdentifiedAsOutInz)
  , fPassTheEventCut(ref.fPassTheEventCut)
  , fRejectKinkMother(ref.fRejectKinkMother)
  , fisppMultiBin(ref.fisppMultiBin)
  , fPbPbUserCentralityBinning(ref.fPbPbUserCentralityBinning)
  , fisNonHFEsystematics(ref.fisNonHFEsystematics)
  , fSpecialTrigger(ref.fSpecialTrigger)
  , fCentralityF(ref.fCentralityF)
  , fCentralityPercent(ref.fCentralityPercent)
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
  , fBackgroundSubtraction(NULL)
  , fQA(NULL)
  , fOutput(NULL)
  , fHistMCQA(NULL)
  , fHistSECVTX(NULL)
  , fHistELECBACKGROUND(NULL)
  , fQACollection(NULL)
  , fQAAODCollection(NULL)
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
  target.fAODMCHeader = fAODMCHeader;
  target.fAODArrayMCInfo = fAODArrayMCInfo;
  target.fQAlevel = fQAlevel;
  target.fPlugins = fPlugins;
  target.fFillSignalOnly = fFillSignalOnly;
  target.fFillNoCuts = fFillNoCuts;
  target.fUseFilterAOD = fUseFilterAOD;
  target.fApplyCutAOD = fApplyCutAOD;
  target.fFilter = fFilter;
  target.fBackGroundFactorApply = fBackGroundFactorApply;
  target.fRemovePileUp = fRemovePileUp;
  target.fIdentifiedAsPileUp = fIdentifiedAsPileUp;
  target.fIdentifiedAsOutInz = fIdentifiedAsOutInz;
  target.fPassTheEventCut = fPassTheEventCut;
  target.fRejectKinkMother = fRejectKinkMother;
  target.fisppMultiBin =   fisppMultiBin;
  target.fPbPbUserCentralityBinning = fPbPbUserCentralityBinning;
  target.fisNonHFEsystematics = fisNonHFEsystematics;
  target.fSpecialTrigger = fSpecialTrigger;
  target.fCentralityF = fCentralityF;
  target.fCentralityPercent = fCentralityPercent;
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
  target.fBackgroundSubtraction = fBackgroundSubtraction;
  target.fQA = fQA;
  target.fOutput = fOutput;
  target.fHistMCQA = fHistMCQA;
  target.fHistSECVTX = fHistSECVTX;
  target.fHistELECBACKGROUND = fHistELECBACKGROUND;
  target.fQACollection = fQACollection;
  target.fQAAODCollection = fQAAODCollection;
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
  if(fBackgroundSubtraction) delete fBackgroundSubtraction;
  if(fSpecialTrigger) delete fSpecialTrigger;
  // Delete output objects only if we are not running in PROOF mode because otherwise this produces a crash during merging
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(mgr && mgr->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
    if(fPIDqa) delete fPIDqa;
    if(fOutput) delete fOutput;
    if(fQA) delete fQA;
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
  
  // Make lists for Output
  if(!fQA) fQA = new TList;
  fQA->SetOwner();
  if(!fOutput) fOutput = new TList;
  fOutput->SetOwner();
  
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
  if(IsAODanalysis()) {
    printf("AOD filter: %s \n", fUseFilterAOD ? "Yes" : "No");
    if(fUseFilterAOD) printf("AOD filter used: %d \n", fFilter);
    // First Part: Make QA histograms
    fQAAODCollection = new AliHFEcollection("TaskQAAOD", "QA histos from the AOD Electron Task");
    fQAAODCollection->CreateTH1F("Filterorigin", "AOD filter of tracks at the origin", 21, -1, 20);
    fQAAODCollection->CreateTH1F("Filterend", "AOD filter of tracks after all cuts", 21, -1, 20);
    fQA->Add(fQAAODCollection);
  }
  printf("MC Data available %s\n", HasMCData() ? "Yes" : "No");

  // Enable Trigger Analysis
  fTriggerAnalysis = new AliTriggerAnalysis;
  fTriggerAnalysis->EnableHistograms();
  fTriggerAnalysis->SetAnalyzeMC(HasMCData());

  // First Part: Make QA histograms
  fQACollection = new AliHFEcollection("TaskQA", "QA histos from the Electron Task");
  fQACollection->CreateTH1F("nElectronTracksEvent", "Number of Electron Candidates", 100, 0, 100);
  fQACollection->CreateTH1F("nElectron", "Number of electrons", 100, 0, 100);
  fQACollection->CreateTH2F("radius", "Production Vertex", 100, 0.0, 5.0, 100, 0.0, 5.0);
  fQACollection->CreateTH2F("TPCdEdxBeforePID", "TPC dE/dx; p (GeV/c); TPC dE/dx (a.u.)", 1000, 0., 10., 200, 0., 200.); 
  fQACollection->CreateTH2F("TPCnSigmaBeforePID", "TPC dE/dx; p (GeV/c); TPC dE/dx - <TPC dE/dx>|_{el} (#sigma)", 1000, 0., 10., 1000, -10., 10.);
 
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

  // Background subtraction-------------------------------------------------------------------
  if (GetPlugin(kNonPhotonicElectron)) {
    if(!fBackgroundSubtraction) fBackgroundSubtraction = new AliHFENonPhotonicElectron();
    fBackgroundSubtraction->Init();
    fOutput->Add(fBackgroundSubtraction->GetListOutput());
  }
  //------------------------------------------------------------------------------------------


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
    if(TestBit(kTreeStream)){
      fMCQA->EnableDebugStreamer();
    }
    fMCQA->CreatDefaultHistograms(fHistMCQA);
    fMCQA->SetBackgroundWeightFactor(fElecBackgroundFactor[0][0][0],fBinLimit);
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

  PrintStatus();
  // Done!!!
  PostData(1, fOutput);
  PostData(2, fQA);
}

//____________________________________________________________
void AliAnalysisTaskHFE::UserExec(Option_t *){
  //
  // Run the analysis
  // 

  //printf("test00\n");

  AliDebug(3, "Starting Single Event Analysis");
  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    //printf("Reconstructed Event not available");
    return;
  }
  if(HasMCData() && IsESDanalysis()){
    AliDebug(4, Form("MC Event: %p", fMCEvent));
    if(!fMCEvent){
      AliError("No MC Event, but MC Data required");
      //printf("No MC Event, but MC Data required");
      return;
    }
  }
  if(!fCuts){
    AliError("HFE cuts not available");
    //printf("HFE cuts not available");
    return;
  }
  if(!fPID->IsInitialized()){
    // Initialize PID with the given run number
    fPID->InitializePID(fInputEvent->GetRunNumber());
  }

  //printf("test0\n");

  // Initialize hadronic background from OADB Container
  AliDebug(2, Form("Apply background factors: %s, OADB Container %p", fBackGroundFactorApply ? "Yes" : "No", fHadronBackgroundOADB));
  if(fBackGroundFactorApply && !TestBit(kBackgroundInitialized)){ 
    AliDebug(2, "Initializing Background from OADB");
    if(!InitializeHadronBackground(fInputEvent->GetRunNumber())) AliError("Failed initializing hadronic background parameterization from OADB");
    else AliDebug(2, "Successfully loaded Background from OADB");
    SetBit(kBackgroundInitialized); 
  }

  //printf("test1\n");

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
    if(GetPlugin(kNonPhotonicElectron) && fBackgroundSubtraction) fBackgroundSubtraction->SetMCEvent(fMCEvent);
  }

  if(IsAODanalysis() && HasMCData()){
    // take MC info
    AliAODEvent *aodE = dynamic_cast<AliAODEvent *>(fInputEvent);
    if(!aodE){ 
      AliError("No AOD Event");
      return;
    }
    fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fInputEvent->FindListObject(AliAODMCHeader::StdBranchName()));
    if(!fAODMCHeader){ 
      AliError("No AliAODMCHeader");
      //printf("No AliAODMCHeader");
      return;
    }
    fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fAODArrayMCInfo){ 
      AliError("No AOD MC particles");
      //printf("No AOD MC particles");
      return;
    }
    fSignalCuts->SetMCAODInfo(fAODArrayMCInfo);
    // Background subtraction-------------------------------------------------------------------
    if (GetPlugin(kNonPhotonicElectron)) {
      fBackgroundSubtraction->SetAODArrayMCInfo(fAODArrayMCInfo);
    }
    //------------------------------------------------------------------------------------------
  }

  //printf("test2\n");

  // need the centrality for everything (MC also)
  fCentralityF = -1;
  if(!ReadCentrality()) fCentralityF = -1;
  //printf("pass centrality\n");
  //printf("Reading fCentralityF %d\n",fCentralityF);
  
  // See if pile up and z in the range
  RejectionPileUpVertexRangeEventCut();

  //printf("test3\n");

  // Protect agains missing 
  if(HasMCData()){
    //printf("Has MC data\n");
    fSignalCuts->SetMCEvent(fMCEvent);
    ProcessMC();  // Run the MC loop + MC QA in case MC Data are available
  }
  
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
    pidResponse = AliHFEtools::GetDefaultPID(HasMCData(), fInputEvent->IsA() == AliESDEvent::Class());
  }
  fPID->SetPIDResponse(pidResponse);
  if(fPIDpreselect) fPIDpreselect->SetPIDResponse(pidResponse);

  // Background subtraction-------------------------------------------------------------------
  if(GetPlugin(kNonPhotonicElectron)) fBackgroundSubtraction->InitRun(fInputEvent,pidResponse);
  //------------------------------------------------------------------------------------------

  // Event loop
  if(IsAODanalysis()){
    //printf("test4\n");
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
    //printf("Running post analysis\n");
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
  if(IsESDanalysis()) eventContainer[0] = fMCEvent->GetPrimaryVertex()->GetZ();
  else eventContainer[0] = fAODMCHeader->GetVtxZ();
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
        fMCQA->SetPercentrality(fCentralityPercent);
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
          fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG); // no accept cut
          fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG); // no accept cut
          fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kOthers, AliHFEmcQA::kElectronPDG); // no accept cut
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
  Int_t numberofmctracks = 0;
  if(IsESDanalysis()){
    numberofmctracks = fMCEvent->GetNumberOfTracks();
  }
  else {
    numberofmctracks = fAODArrayMCInfo->GetEntriesFast();
  }
  AliDebug(3, Form("Number of Tracks: %d",numberofmctracks));
  //printf("Number of MC track %d\n",numberofmctracks);
  for(Int_t imc = 0; imc <numberofmctracks; imc++){
    if(IsESDanalysis()) {
      if(!(mctrack = fMCEvent->GetTrack(imc))) continue;
    }
    else {
      if(!(mctrack = (AliVParticle *) fAODArrayMCInfo->At(imc))) continue;
    }
    //printf("Test in ProcessMC\n");
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
  AliDebug(1, Form("Task %s", GetName()));
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
  Double_t dataDca[6]; // [source, pT, dca, centrality]
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
  

  // Background subtraction-------------------------------------------------------------------
  if (GetPlugin(kNonPhotonicElectron)) {
    fBackgroundSubtraction->FillPoolAssociatedTracks(fInputEvent,fCentralityF);
  }
  //------------------------------------------------------------------------------------------

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
      AliDebug(1, "V0 PID done");
    }
 

    //Fill non-HFE source containers at reconstructed events cut step
    AliDebug(3, Form("Doing track %d, %p", itrack, track));


    //////////////////////////////////////
    // preselect
    //////////////////////////////////////
    if(fPIDpreselect || fCutspreselect) {
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
        
    if(fisNonHFEsystematics && IsPbPb())  {
      //FillProductionVertex(track);     
      if(fMCQA && signal){ 
	fMCQA->SetCentrality(fCentralityF);
	if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) == 11)){
	  Double_t weightElecBgV0[kBgLevels] = {0.,0.,0.};
          
	  fMCQA->SetTrkKine(track->Pt(),track->Eta(), track->Phi());
	  fMCQA->SetContainerStep(4);
          
	  weightElecBgV0[0] = fMCQA->GetWeightFactor(mctrack, 0); // positive:conversion e, negative: nonHFE 
          
	  //Fill additional containers for electron source distinction
	  Int_t elecSource = 0;
	  elecSource = fMCQA->GetElecSource(mctrack->Particle());
	  const Char_t *sourceName[kElecBgSpecies]={"Pion","Eta","Omega","Phi","EtaPrime","Rho"};
	  const Char_t *levelName[kBgLevels]={"Best","Lower","Upper"};
	  Int_t iName = 0;
	  for(Int_t iSource = AliHFEmcQA::kPi0; iSource <=  AliHFEmcQA::kGammaRho0; iSource++){
	    if((iSource == AliHFEmcQA::kElse)||(iSource == AliHFEmcQA::kMisID)) continue;
	    if(elecSource == iSource){
	      
	      if(weightElecBgV0[0]>0){ 
		fVarManager->FillContainer(fContainer, Form("conversionElecs%s%s",sourceName[iName], levelName[0]), 4, kTRUE, weightElecBgV0[0]);
	      } 
	      else if(weightElecBgV0[0]<0){ 
		fVarManager->FillContainer(fContainer, Form("mesonElecs%s%s",sourceName[iName], levelName[0]), 4, kTRUE, -1*weightElecBgV0[0]);
	      }           
	      break;
	    }
	    iName++;
	    if(iName == kElecBgSpecies)iName = 0;
	  }
	}
      }
    }
  
    // RecKine: ITSTPC cuts  
    if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    
    
    // RecPrim
    if(fRejectKinkMother) { 
      if(track->GetKinkIndex(0) != 0) continue; } // Quick and dirty fix to reject both kink mothers and daughters
    if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;

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

      if(fMCQA && signal){
        fMCQA->SetCentrality(fCentralityF);
        if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) == 11)){
         Double_t weightElecBgV0[kBgLevels] = {0.,0.,0.};
         Double_t hfeimpactRtmp=0., hfeimpactnsigmaRtmp=0.;
         fExtraCuts->GetHFEImpactParameters(track, hfeimpactRtmp, hfeimpactnsigmaRtmp);
         UChar_t itsPixel = track->GetITSClusterMap();
         Double_t ilyrhit=0, ilyrstat=0;
         for(Int_t ilyr=0; ilyr<6; ilyr++){
           if(TESTBIT(itsPixel, ilyr)) ilyrhit += TMath::Power(2,ilyr);
           if(fExtraCuts->CheckITSstatus(fExtraCuts->GetITSstatus(track,ilyr))) ilyrstat += TMath::Power(2,ilyr);
         }
         fMCQA->SetITSInfo(ilyrhit,ilyrstat);
         fMCQA->SetHFEImpactParameters(hfeimpactRtmp, hfeimpactnsigmaRtmp);
         fMCQA->SetTrkKine(track->Pt(),track->Eta(), track->Phi());
         fMCQA->SetContainerStep(3);
         for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
           weightElecBgV0[iLevel] = fMCQA->GetWeightFactor(mctrack, iLevel); // positive:conversion e, negative: nonHFE 
           if(!fisNonHFEsystematics || IsPbPb())break;   
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
                 if(weightElecBgV0[iLevel]>0){ 
                   fVarManager->FillContainer(fContainer, Form("conversionElecs%s%s",sourceName[iName], levelName[iLevel]), 3, kFALSE, weightElecBgV0[iLevel]);
                 } 
                 else if(weightElecBgV0[iLevel]<0){ 
                   fVarManager->FillContainer(fContainer, Form("mesonElecs%s%s",sourceName[iName], levelName[iLevel]), 3, kFALSE, -1*weightElecBgV0[iLevel]);
                 }
                 if(IsPbPb())break;
               }
               break;
             }
             iName++;
             if(iName == kElecBgSpecies)iName = 0;
           }
         }
         //else{
           if(weightElecBgV0[0]>0) {
	     fVarManager->FillContainer(fContainer, "conversionElecs", 3, kFALSE, weightElecBgV0[0]);
	     fVarManager->FillContainer(fContainer, "conversionElecs", 4, kTRUE, weightElecBgV0[0]);
	   }
           else if(weightElecBgV0[0]<0) {
	     fVarManager->FillContainer(fContainer, "mesonElecs", 3, kFALSE, -1*weightElecBgV0[0]);
	     fVarManager->FillContainer(fContainer, "mesonElecs", 4, kTRUE, -1*weightElecBgV0[0]);
	   }
           //}
        }
      }

      Double_t hfeimpactR4all=0., hfeimpactnsigmaR4all=0.;
      Int_t sourceDca =-1;
      if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) == 211)){
        if(track->Pt()>4.){
          fExtraCuts->GetHFEImpactParameters(track, hfeimpactR4all, hfeimpactnsigmaR4all);
          dataDca[0]=0; //pion
          dataDca[1]=track->Pt();
          dataDca[2]=hfeimpactR4all;
          dataDca[3]=fCentralityF;
          dataDca[4] = v0pid;
          dataDca[5] = double(track->Charge());
          fQACollection->Fill("Dca", dataDca);
        }
      }
      else if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) == 11)){ // to increas statistics for Martin
        if(signal){
          fExtraCuts->GetHFEImpactParameters(track, hfeimpactR4all, hfeimpactnsigmaR4all);
          if(fSignalCuts->IsCharmElectron(track)){
            sourceDca=1;
          }
          else if(fSignalCuts->IsBeautyElectron(track)){
            sourceDca=2;
          }
          else if(fSignalCuts->IsGammaElectron(track)){
            sourceDca=3;
          }
          else if(fSignalCuts->IsNonHFElectron(track)){
            sourceDca=4;
          }
          else if(fSignalCuts->IsJpsiElectron(track)){
            sourceDca=5;
          }
          else {
            sourceDca=6;
          }
          dataDca[0]=sourceDca;
          dataDca[1]=track->Pt();
          dataDca[2]=hfeimpactR4all;
          dataDca[3]=fCentralityF;
          dataDca[4] = v0pid;
          dataDca[5] = double(track->Charge());
          if(signal) fQACollection->Fill("Dca", dataDca);
        }
      }
    }

    if(TMath::Abs(track->Eta()) < 0.5){
      if(track->GetInnerParam())
        fQACollection->Fill("TPCdEdxBeforePID", track->GetInnerParam()->P(), track->GetTPCsignal());
      fQACollection->Fill("TPCnSigmaBeforePID", track->P(), fInputHandler->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron));
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
    
    // Background subtraction------------------------------------------------------------------------------------------
    if (GetPlugin(kNonPhotonicElectron)) {
      Int_t indexmother = -1;
      Int_t mcsource = 1;
      if(HasMCData()) mcsource = fBackgroundSubtraction->FindMother(mctrack->GetLabel(),indexmother);
      fBackgroundSubtraction->LookAtNonHFE(itrack, track, fInputEvent, 1, fCentralityF, -1, mcsource, indexmother);
    }
    //-----------------------------------------------------------------------------------------------------------------

    // Temporary histogram for chi2/ITS cluster
    if(IsPbPb()) {
            TBits shared = track->GetTPCSharedMap();
	    Int_t sharebit=0;
            if(shared.CountBits() >= 2) sharebit=1;

	    Double_t itschi2percluster = 0.0;
	    Double_t itsnbcls = static_cast<Double_t>(track->GetNcls(0));
	    if(itsnbcls > 0) itschi2percluster = track->GetITSchi2()/itsnbcls;

            Double_t itsChi2[7] = {track->Pt(),track->Eta(), track->Phi(),
				   fCentralityF,track->GetTPCsignalN(), sharebit,itschi2percluster};
            fQACollection->Fill("fChi2perITScluster", itsChi2);
    }
    else{

      Double_t itschi2percluster = 0.0;
      Double_t itsnbcls = static_cast<Double_t>(track->GetNcls(0));
      if(itsnbcls > 0) itschi2percluster = track->GetITSchi2()/itsnbcls;

      Double_t itsChi2[3] = {track->Pt(), fCentralityF, itschi2percluster};
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
            fQACollection->Fill("Ke3Kecorr",mctrack->Pt(),mctrackmother->Pt());
          else if(TMath::Abs(mctrackmother->Particle()->GetPdgCode())==130)
            fQACollection->Fill("Ke3K0Lecorr",mctrack->Pt(),mctrackmother->Pt());
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

    if (GetPlugin(kDEstep)) { 
      Double_t weightElecBgV0[kBgLevels] = {0.,0.,0.,};
      Int_t elecSource = 0;
      Double_t hfeimpactR=0., hfeimpactnsigmaR=0.;
      fExtraCuts->GetHFEImpactParameters(track, hfeimpactR, hfeimpactnsigmaR);
      if(HasMCData())
      {
        if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) != 11)){
            fQACollection->Fill("hadronsBeforeIPcut",track->Pt());
        } 
        if(fMCQA && signal) {
          
          fMCQA->SetContainerStep(0);
          for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
            weightElecBgV0[iLevel] = fMCQA->GetWeightFactor(mctrack, iLevel); // positive:conversion e, negative: nonHFE 
            if(!fisNonHFEsystematics || IsPbPb())break;        
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
                  if(IsPbPb())break;
                }
                break;
              }
              iName++;
              if(iName == kElecBgSpecies)iName = 0;
            }
          }
          //else{
          if(weightElecBgV0[0]>0) {
	    fVarManager->FillContainer(fContainer, "conversionElecs", 0, kFALSE, weightElecBgV0[0]);
	    fVarManager->FillContainer(fContainer, "conversionElecs", 5, kTRUE, weightElecBgV0[0]);
	  }
          else if(weightElecBgV0[0]<0) {
	    fVarManager->FillContainer(fContainer, "mesonElecs", 0, kFALSE, -1*weightElecBgV0[0]);
	    fVarManager->FillContainer(fContainer, "mesonElecs", 5, kTRUE, -1*weightElecBgV0[0]);
	  }  
          //}
          if(bTagged){ // bg estimation for the secondary vertex tagged signals
            if(weightElecBgV0[0]>0) fVarManager->FillContainer(fContainer, "conversionElecs", 2, kFALSE, weightElecBgV0[0]);
            else if(weightElecBgV0[0]<0) fVarManager->FillContainer(fContainer, "mesonElecs", 2, kFALSE, -1*weightElecBgV0[0]);
          }
        }
      } // end of MC

      dataDca[0]=-1; //for data, don't know the origin
      dataDca[1]=track->Pt();
      dataDca[2]=hfeimpactR;
      dataDca[3]=fCentralityF;
      dataDca[4] = v0pid;
      dataDca[5] = double(track->Charge());
      if (!HasMCData()) fQACollection->Fill("Dca", dataDca);

      // Fill Containers for impact parameter analysis
      if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsDca + AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack,track)) continue;
      if(signal) {
        // Apply weight for background contamination after ip cut
        if(fBackGroundFactorApply) {
              fWeightBackGround =  fkBackGroundFactorArray[0]->Eval(TMath::Abs(track->P())); // pp case
              if(fWeightBackGround < 0.0) fWeightBackGround = 0.0;
              else if(fWeightBackGround > 1.0) fWeightBackGround = 1.0;
              // weightBackGround as special weight
              fVarManager->FillContainer(fContainer, "hadronicBackground", 2, kFALSE, fWeightBackGround);
        }
      }

      if(HasMCData()){
        if(fMCQA && signal) {
          fMCQA->SetContainerStep(1);
          for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
            weightElecBgV0[iLevel] = fMCQA->GetWeightFactor(mctrack, iLevel); // positive:conversion e, negative: nonHFE 
            if(!fisNonHFEsystematics || IsPbPb())break;        
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
                  if(IsPbPb())break;
                }
                break;
              }
              iName++;
              if(iName == kElecBgSpecies)iName = 0;
            }
          }
          // else{
            if(weightElecBgV0[0]>0) {
	      fVarManager->FillContainer(fContainer, "conversionElecs", 1, kFALSE, weightElecBgV0[0]);
	      fVarManager->FillContainer(fContainer, "conversionElecs", 6, kTRUE, weightElecBgV0[0]);
	    }
            else if(weightElecBgV0[0]<0) {
	      fVarManager->FillContainer(fContainer, "mesonElecs", 1, kFALSE, -1*weightElecBgV0[0]);
	      fVarManager->FillContainer(fContainer, "mesonElecs", 6, kTRUE, -1*weightElecBgV0[0]);
            }
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
  //printf("Process AOD\n");
  AliDebug(3, "Processing AOD Event");
  Double_t eventContainer[4];
  eventContainer[0] = 0.0;
  if(HasMCData()) eventContainer[0] = fVz;
  else {
    if(fInputEvent->GetPrimaryVertex()) eventContainer[0] = fInputEvent->GetPrimaryVertex()->GetZ();
  }
  eventContainer[1] = 1.; // No Information available in AOD analysis, assume all events have V0AND
  eventContainer[2] = fCentralityF; 
  eventContainer[3] = fContributors; 
  
  //printf("value event container %f, %f, %f, %f\n",eventContainer[0],eventContainer[1],eventContainer[2],eventContainer[3]);

  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
  if(!fAOD){
    AliError("AOD Event required for AOD Analysis");
      return;
  }
  
  //printf("Will fill\n");
  //
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoCut);
  //printf("Fill\n");
  //
  if(fIdentifiedAsPileUp) return; 
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepRecNoPileUp);

  //
  if(fIdentifiedAsOutInz) return;
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepZRange);  

  //
  if(!fPassTheEventCut) return;
  fCFM->GetEventContainer()->Fill(eventContainer, AliHFEcuts::kEventStepReconstructed);
  //printf("pass\n");

  fContainer->NewEvent();

  fCFM->SetRecEventInfo(fAOD);

  // Look for kink mother
  Int_t numberofvertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
  for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
    if(!aodvertex) continue;
    if(aodvertex->GetType()==AliAODVertex::kKink) {
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      listofmotherkink[numberofmotherkink] = idmother;
      //printf("ID %d\n",idmother);
      numberofmotherkink++;
    }
  }
  //printf("Number of kink mother in the events %d\n",numberofmotherkink);

  // Background subtraction-------------------------------------------------------------------
  if (GetPlugin(kNonPhotonicElectron)) {
    fBackgroundSubtraction->FillPoolAssociatedTracks(fInputEvent,fCentralityF);
  }
  //------------------------------------------------------------------------------------------

  // Loop over tracks
  AliAODTrack *track = NULL;
  AliAODMCParticle *mctrack = NULL;
  Double_t dataE[6]; // [pT, eta, Phi, Charge, type, 'C' or 'B']
  Int_t nElectronCandidates = 0;
  Int_t pid;
  Bool_t signal;
  //printf("Number of track %d\n",(Int_t) fAOD->GetNumberOfTracks());
  for(Int_t itrack = 0; itrack < fAOD->GetNumberOfTracks(); itrack++){
    track = fAOD->GetTrack(itrack); mctrack = NULL;
    if(!track) continue;
    // Begining
    Bool_t passone = kFALSE;  
    fQAAODCollection->Fill("Filterorigin", -1);
    for(Int_t k=0; k<20; k++) {
      Int_t u = 1<<k;     
      if((track->TestFilterBit(u))) {
	fQAAODCollection->Fill("Filterorigin", k);
      	passone = kTRUE;
      }
    }
    //if(!passone) printf("what is the filter %d\n",track->GetFilterMap());

    if(fUseFilterAOD){
      //printf("Filter of the track  %d\n",track->GetFilterMap());
      if(!(track->TestFilterBit(fFilter))) continue;  // Only process AOD tracks where the HFE is set
    }
    //printf("Pass the flag\n");
    
    signal = kTRUE;
    if(HasMCData()){

      Int_t label = TMath::Abs(track->GetLabel());
      if(label && label < fAODArrayMCInfo->GetEntriesFast())
        mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
        if(fFillSignalOnly && !fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE;
    }
    
    fVarManager->NewTrack(track, mctrack, fCentralityF, -1, signal);
    
    if(fFillNoCuts) {
      if(signal || !fFillSignalOnly){
        fVarManager->FillContainer(fContainer, "recTrackContReco", AliHFEcuts::kStepRecNoCut, kFALSE);
        fVarManager->FillContainer(fContainer, "recTrackContMC", AliHFEcuts::kStepRecNoCut, kTRUE);
      }
    }

    if(fApplyCutAOD) {
      //printf("Apply cuts\n");
      // RecKine: ITSTPC cuts  
      if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;

      // Reject kink mother
      if(fRejectKinkMother) {
	Bool_t kinkmotherpass = kTRUE;
	for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
	  if(track->GetID() == listofmotherkink[kinkmother]) {
	    kinkmotherpass = kFALSE;
	    continue;
	  }
	}
	if(!kinkmotherpass) continue;
      }
      
      // RecPrim
      if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
      
      // HFEcuts: ITS layers cuts
      if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
      
      // HFE cuts: TOF PID and mismatch flag
      if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTOF, track)) continue;
      
      // HFE cuts: TPC PID cleanup
      if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
      
      // HFEcuts: Nb of tracklets TRD0
      if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTRD, track)) continue;
    }

    // Fill correlation maps before PID
    if(signal && fContainer->GetCorrelationMatrix("correlationstepbeforePID")) {
      //printf("Fill correlation maps before PID\n");
      fVarManager->FillCorrelationMatrix(fContainer->GetCorrelationMatrix("correlationstepbeforePID"));
    }

    //printf("Will process to PID\n");

    // track accepted, do PID
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
    hfetrack.SetRecTrack(track);
    if(HasMCData()) hfetrack.SetMCTrack(mctrack);
    hfetrack.SetCentrality(fCentralityF);
    if(IsPbPb()) hfetrack.SetPbPb();
    else hfetrack.SetPP();
    fPID->SetVarManager(fVarManager);
    if(!fPID->IsSelected(&hfetrack, fContainer, "recTrackCont", fPIDqa)) continue;   
    // we will do PID here as soon as possible

    // Background subtraction----------------------------------------------------------------------------------------------
    if (GetPlugin(kNonPhotonicElectron)) {
      Int_t indexmother = -1;
      Int_t mcsource = 1;
      if(HasMCData()) mcsource = fBackgroundSubtraction->FindMother(mctrack->GetLabel(),indexmother);
      fBackgroundSubtraction->LookAtNonHFE(itrack, track, fInputEvent, 1, fCentralityF, -1,mcsource, indexmother);
    }
    //---------------------------------------------------------------------------------------------------------------------

    // end AOD QA
    fQAAODCollection->Fill("Filterend", -1);  
    for(Int_t k=0; k<20; k++) {
      Int_t u = 1<<k;
      if((track->TestFilterBit(u))) {
	fQAAODCollection->Fill("Filterend", k);
      }
    }
       
    // Apply weight for background contamination
    //Double_t weightBackGround = 1.0;
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

  Double_t vertex[3] = {0.,0.,0.}; // Production vertex cut to mask gammas which are NOT supposed to have hits in the first ITS layer(s)
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

  //printf("MC Generated\n");
  if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, track)) return kFALSE;
  //printf("MC Generated pass\n");
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

  if(fCutspreselect) {
    //printf("test preselect\n");
    if(!fCutspreselect->IsSelected(track)) survived=kFALSE;
  }
  //printf("survived %d\n",(Int_t)survived);
  
  if(survived && fPIDpreselect){
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
  
  fContainer->CreateContainer("hadronicBackground", "Container for Hadronic Background", 3);
  fContainer->CreateContainer("recTrackContDEReco", "Container for displaced electron analysis with Reco information", 1);
  fContainer->CreateContainer("recTrackContDEMC", "Container for displaced electron analysis with MC information", 1);
  fContainer->CreateContainer("recTrackContSecvtxReco", "Container for secondary vertexing analysis with Reco information", 1);
  fContainer->CreateContainer("recTrackContSecvtxMC", "Container for secondary vertexing analysis with MC information", 1);

  if(HasMCData()){
    fContainer->CreateContainer("conversionElecs", "Container for weighted conversion electrons",7);
    fContainer->CreateContainer("mesonElecs", "Container for weighted electrons from meson decays",7);
    fContainer->Sumw2("conversionElecs");
    fContainer->Sumw2("mesonElecs");
   
    if(fisNonHFEsystematics){
      const Char_t *sourceName[kElecBgSpecies]={"Pion","Eta","Omega","Phi","EtaPrime","Rho"};
      const Char_t *levelName[kBgLevels]={"Best","Lower","Upper"};
      for(Int_t iSource = 0; iSource < kElecBgSpecies; iSource++){
        for(Int_t iLevel = 0; iLevel < kBgLevels; iLevel++){
          fContainer->CreateContainer(Form("conversionElecs%s%s",sourceName[iSource],levelName[iLevel]), Form("Container for weighted conversion electrons from %s grandm., %s level",sourceName[iSource],levelName[iLevel]),5);
          fContainer->CreateContainer(Form("mesonElecs%s%s",sourceName[iSource],levelName[iLevel]), Form("Container for weighted electrons from %s decays, %s level",sourceName[iSource],levelName[iLevel]),5);
          fContainer->Sumw2(Form("conversionElecs%s%s",sourceName[iSource],levelName[iLevel]));
          fContainer->Sumw2(Form("mesonElecs%s%s",sourceName[iSource],levelName[iLevel]));
          if(IsPbPb())break;
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

  TObjArray *array = fVarManager->GetVariables();
  Int_t nvars = array->GetEntriesFast();
  for(Int_t v = 0; v < nvars; v++) {
    AliHFEvarManager::AliHFEvariable *variable = (AliHFEvarManager::AliHFEvariable *) array->At(v);
    if(!variable) continue;
    TString name(((AliHFEvarManager::AliHFEvariable *)variable)->GetName());
    if(!name.CompareTo("pt")) {
      const Int_t nBinPt  = variable->GetNumberOfBins();
      const Double_t *kPtRange = variable->GetBinning();

      fQACollection->CreateTH1Farray("hadronsBeforeIPcut", "Hadrons before IP cut", nBinPt, kPtRange);
      fQACollection->CreateTH1Farray("hadronsAfterIPcut", "Hadrons after IP cut", nBinPt, kPtRange);

      fQACollection->CreateTH2Farray("Ke3Kecorr", "Ke3 decay e and K correlation; Ke3K p_{t}; Ke3e p_{t}; ", nBinPt, kPtRange, 20,0.,20.);
      fQACollection->CreateTH2Farray("Ke3K0Lecorr", "Ke3 decay e and K0L correlation; Ke3K0L p_{t}; Ke3e p_{t}; ", nBinPt, kPtRange, 20,0.,20.);
      fQACollection->CreateTH1Farray("Kptspectra", "Charged Kaons: MC p_{t} ", nBinPt, kPtRange);
      fQACollection->CreateTH1Farray("K0Lptspectra", "K0L: MC p_{t} ", nBinPt, kPtRange);

      const Double_t kDCAbound[2] = {-0.2, 0.2};

      const Int_t nDimDca=6;
      const Int_t nBinDca[nDimDca] = { 8, nBinPt, 800, 12,  6, 2};
      Double_t minimaDca[nDimDca]  = { -1., 0., kDCAbound[0], -1., -1, -1.1};
      Double_t maximaDca[nDimDca]  = { 7., 20., kDCAbound[1], 11.,  5, 1.1};

      Double_t *sourceBins = AliHFEtools::MakeLinearBinning(nBinDca[0], minimaDca[0], maximaDca[0]);
      Double_t *dcaBins = AliHFEtools::MakeLinearBinning(nBinDca[2], minimaDca[2], maximaDca[2]);
      Double_t *centralityBins = AliHFEtools::MakeLinearBinning(nBinDca[3], minimaDca[3], maximaDca[3]);
      Double_t *v0PIDBins = AliHFEtools::MakeLinearBinning(nBinDca[4], minimaDca[4], maximaDca[4]);
      Double_t *chargeBins = AliHFEtools::MakeLinearBinning(nBinDca[5], minimaDca[5], maximaDca[5]);

      fQACollection->CreateTHnSparseNoLimits("Dca", "Dca; source (0-all, 1-charm,etc); pT [GeV/c]; dca; centrality bin; v0pid; charge", nDimDca, nBinDca);
      ((THnSparse*)(fQACollection->Get("Dca")))->SetBinEdges(0, sourceBins);
      ((THnSparse*)(fQACollection->Get("Dca")))->SetBinEdges(1, kPtRange);
      ((THnSparse*)(fQACollection->Get("Dca")))->SetBinEdges(2, dcaBins);
      ((THnSparse*)(fQACollection->Get("Dca")))->SetBinEdges(3, centralityBins);
      ((THnSparse*)(fQACollection->Get("Dca")))->SetBinEdges(4, v0PIDBins);
      ((THnSparse*)(fQACollection->Get("Dca")))->SetBinEdges(5, chargeBins);

      break;
    }  
  }

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
    case kNonPhotonicElectron: SETBIT(fPlugins, plug); break; 
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
  
  Float_t fCentralityLimitstemp[12];
  Float_t fCentralityLimitsdefault[12]= {0.,5.,10., 20., 30., 40., 50., 60.,70.,80., 90., 100.};
  if(!fPbPbUserCentralityBinning) memcpy(fCentralityLimitstemp,fCentralityLimitsdefault,sizeof(fCentralityLimitsdefault));
  else memcpy(fCentralityLimitstemp,fCentralityLimits,sizeof(fCentralityLimitsdefault));
  

  Int_t bin = -1;
  if(IsPbPb()) {
    // Centrality
    AliCentrality *centrality = fInputEvent->GetCentrality();
    fCentralityPercent = centrality->GetCentralityPercentile("V0M");
    //printf("centrality %f\n",fCentralityPercent);

    for(Int_t ibin = 0; ibin < 11; ibin++){
      if(fCentralityPercent >= fCentralityLimitstemp[ibin] && fCentralityPercent < fCentralityLimitstemp[ibin+1]){
        bin = ibin;
	//printf("test bin %f, low %f, high %f, %d\n",fCentralityPercent,fCentralityLimitstemp[ibin],fCentralityLimitstemp[ibin+1],ibin);
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
    if( contributorstemp <=  0) {
      fContributors =  0.5;
      //printf("Number of contributors %d and vz %f\n",contributorstemp,vtx->GetZ());
    }
    else fContributors = 1.5;
    //printf("Number of contributors %d\n",contributorstemp);
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
   fIdentifiedAsPileUp = kFALSE;
   if(fRemovePileUp && fAOD->IsPileupFromSPD()) fIdentifiedAsPileUp = kTRUE; 
   // Z vertex
   fIdentifiedAsOutInz = kFALSE;
   //printf("Z vertex %f and out %f\n",fAOD->GetPrimaryVertex()->GetZ(),fCuts->GetVertexRange());
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

