#include "TString.h"
class AliAnalysisDataContainer;
class TNamed;
AliAnalysisTaskMeanPtV2Corr* AddTaskMeanPtV2Corr(TString name = "name", Bool_t IsMC=kFALSE, TString stage = "weights",
                                                  TString weightPath = "", TString meanPtPath="", TString NUAPath="")
{
  Int_t StageSwitch = 0;
  if(stage.Contains("weights")) StageSwitch=1;
  if(stage.Contains("meanpt")) StageSwitch=2;
  if(stage.Contains("full")) StageSwitch=3;
  if(StageSwitch==0) return 0;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  if(IsMC) {
    if(!mgr->GetMCtruthEventHandler()) {
      Error("AddTaskMeanPtV2Corr","Could not get MC truth handler");
      return NULL;
    };
    AliMCEventHandler *handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    handler->SetReadTR(kTRUE);
  };
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskMeanPtV2Corr* task = new AliAnalysisTaskMeanPtV2Corr(name.Data(), IsMC, stage);
  if(!task)
    return 0x0;
  mgr->AddTask(task); // add your task to the manager
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0);
  //Producing weights
  if(StageSwitch==1) {
    AliAnalysisDataContainer *weightCont = mgr->CreateContainer("WeightList",TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,1,weightCont);
    return task;
  }
  //Load input mean pt
  if(StageSwitch==2) {
    TObjArray *AllContainers = mgr->GetContainers();
    if(!AllContainers->FindObject("Weights")) {
      if(weightPath.IsNull()) AliFatal("Weight path not provided!\n");
      if(weightPath.Contains("alien:")) TGrid::Connect("alien:");
      TFile *tfWeights = TFile::Open(weightPath.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MergedWeights.root"
      if(!tfWeights) AliFatal("Could not open weights file\n");
      if(tfWeights->IsZombie()) AliFatal("Weight file is a zombie\n");
      TList *fList = (TList*)tfWeights->Get("WeightList");
      if(!fList) { AliFatal("Could not fetch the weight list!\n"); return; };
      AliAnalysisDataContainer *cWeights = mgr->CreateContainer("Weights",TList::Class(), AliAnalysisManager::kInputContainer);
      cWeights->SetData(fList);
      mgr->ConnectInput(task,1,cWeights);
    };
    AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer("MPTProfileList", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputMPT);
    AliAnalysisDataContainer *multiDist  = mgr->CreateContainer("multiDist",TH1D::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,2,multiDist);
    return task;
  };
  //Full
  if(StageSwitch==3) {
    TObjArray *AllContainers = mgr->GetContainers();
    if(!AllContainers->FindObject("Weights")) {
      if(weightPath.IsNull()) AliFatal("Weight path not provided!\n");
      if(weightPath.Contains("alien:")) TGrid::Connect("alien:");
      TFile *tfWeights = TFile::Open(weightPath.Data());
      TList *fList = (TList*)tfWeights->Get("WeightList");
      AliAnalysisDataContainer *cWeights = mgr->CreateContainer("WeightList",TList::Class(), AliAnalysisManager::kInputContainer);
      cWeights->SetData(fList);
      mgr->ConnectInput(task,1,cWeights);
    };
    if(!AllContainers->FindObject("InputMeanPt")) {
      if(meanPtPath.IsNull()) AliFatal("Mean pT path not provided!\n");
      if((!weightPath.Contains("alien:")) && meanPtPath.Contains("alien:")) TGrid::Connect("alien:"); //Only connect if not connected yet
      TFile *tfMPT = TFile::Open(meanPtPath.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MeanPts_05_20.root"
      TList *fMPTList = (TList*)tfMPT->Get("MPTProfileList");
      if(!fMPTList) AliFatal("fMPT list from file not fetcehd!");
      AliAnalysisDataContainer *cInMPT = mgr->CreateContainer("InputMeanPt",TList::Class(), AliAnalysisManager::kInputContainer);
      cInMPT->SetData(fMPTList);
      mgr->ConnectInput(task,2,cInMPT);
    };
    if(!AllContainers->FindObject("AdHocNUA")) {
      if(NUAPath.IsNull()) AliFatal("AdHoc NUA path not provided!\n");
      if((!weightPath.Contains("alien:")) && !meanPtPath.Contains("alien:") && NUAPath.Contains("alien:") ) TGrid::Connect("alien:"); //Only connect if not connected yet
      TFile *tfNUA = TFile::Open(NUAPath.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MeanPts_05_20.root"
      TList *fNUAList = (TList*)tfNUA->Get("PIDWeights");
      AliAnalysisDataContainer *cInNUA = mgr->CreateContainer("AdHocNUA",TList::Class(), AliAnalysisManager::kInputContainer);
      cInNUA->SetData(fNUAList);
      mgr->ConnectInput(task,3,cInNUA);
    };

    AliAnalysisDataContainer *cOutputCOV = mgr->CreateContainer("MPTDiff",TProfile::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputCOV);
    AliAnalysisDataContainer *cOutputFC  = mgr->CreateContainer("FlowCont",AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,2,cOutputFC);
    AliAnalysisDataContainer *cOutputFC  = mgr->CreateContainer("Covariance",TProfile::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,3,cOutputFC);
    return task;
  };
  return 0;
}
