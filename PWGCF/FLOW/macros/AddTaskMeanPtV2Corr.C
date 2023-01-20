#include "TString.h"
#include "TGrid.h"
class AliAnalysisDataContainer;
class TNamed;
Bool_t ConnectToGrid() {
  if(!gGrid) TGrid::Connect("alien:");
  if(!gGrid) {printf("Task requires connection to grid, but it could not be established!\n"); return kFALSE; };
  return kTRUE;
}
AliAnalysisTaskMeanPtV2Corr* AddTaskMeanPtV2Corr(TString name, Bool_t IsMC, TString stage,
                                                  TString efficiencyPath, TString meanPtPath, TString NUAPath, TString subfix1, TString subfix2)
{
  Int_t StageSwitch = 0;
  printf("Stage switch name: %s\n",stage.Data());
  if(stage.Contains("weights")) StageSwitch=1;
  if(stage.Contains("Efficiency")) StageSwitch=7;
  if(stage.Contains("CovSkipMpt")) StageSwitch=9;
  if(stage.Contains("EfTest")) StageSwitch=10;
  if(StageSwitch==0) return 0;
  TString l_ContName(subfix1);
  if(!subfix2.IsNull()) l_ContName+="_"+subfix2;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  // if(IsMC) {
  //   if(!mgr->GetMCtruthEventHandler()) {
  //     Error("AddTaskMeanPtV2Corr","Could not get MC truth handler");
  //     return NULL;
  //   };
  //   AliMCEventHandler *handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  //   handler->SetReadTR(kTRUE);
  // };
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskMeanPtV2Corr* task = new AliAnalysisTaskMeanPtV2Corr(name.Data(), IsMC, stage, l_ContName);
  if(!task)
    return 0x0;
  mgr->AddTask(task); // add your task to the manager
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0);
  //Producing weights
  if(StageSwitch==1) {
    AliAnalysisDataContainer *weightCont = mgr->CreateContainer(Form("WeightList%s",l_ContName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,1,weightCont);
    return task;
  }
  if(StageSwitch==7) { //Producing Pt spectra with filter bit
    AliAnalysisDataContainer *spectraCont = mgr->CreateContainer(Form("SpectraList%s",l_ContName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,1,spectraCont);
    return task;
  }
  //Full
  if(StageSwitch==9) {
    TObjArray *AllContainers = mgr->GetContainers();
    // AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer(Form("MPTProfileList%s",l_ContName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    Bool_t gridConnected=kFALSE;
    if(!IsMC) {
      if(!AllContainers->FindObject("Efficiency")) {
        if(efficiencyPath.IsNull()) { printf("Efficiency path not provided!\n"); return 0; };
        if(efficiencyPath.Contains("alien:")) if(!ConnectToGrid()) return 0;//{ TGrid::Connect("alien:"); gridConnected = kTRUE; };
        TFile *tfEff = TFile::Open(efficiencyPath.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MergedWeights.root"
        if(!tfEff) { printf("Could not open efficiency file\n"); return 0; };
        if(tfEff->IsZombie()) { printf("Efficiency file is a zombie\n"); return 0; };
        TList *fEffList = new TList();
        TList *keyList = tfEff->GetListOfKeys();
        for(Int_t i=0;i<keyList->GetEntries();i++) fEffList->Add((TList*)tfEff->Get(keyList->At(i)->GetName()));
        // TList *fList = (TList*)tfEff->Get("EffAndFD");
        // if(!fList) { printf("Could not fetch the efficiency list!\n"); return 0; };
        AliAnalysisDataContainer *cEff = mgr->CreateContainer("Efficiency",TList::Class(), AliAnalysisManager::kInputContainer);
        cEff->SetData(fEffList);
        mgr->ConnectInput(task,1,cEff);
      } else mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("Efficiency"));
      if(!AllContainers->FindObject("Weights")) {
        if(NUAPath.IsNull()) { printf("Weight path not provided!\n"); return 0; };
        if(NUAPath.Contains("alien:")) if(!ConnectToGrid()) return 0;//{ TGrid::Connect("alien:"); gridConnected = kTRUE; };
        TFile *tfWeights = TFile::Open(NUAPath.Data());
        TList *fList = (TList*)tfWeights->Get("WeightList");
        AliAnalysisDataContainer *cWeights = mgr->CreateContainer("Weights",TList::Class(), AliAnalysisManager::kInputContainer);
        cWeights->SetData(fList);
        mgr->ConnectInput(task,2,cWeights);
      } else mgr->ConnectInput(task,2,(AliAnalysisDataContainer*)AllContainers->FindObject("Weights"));
    };
    AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer(Form("MPTDiff%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputMPT);
    AliAnalysisDataContainer *cOutputFC  = mgr->CreateContainer(Form("FlowCont%s",l_ContName.Data()),AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,2,cOutputFC);
    AliAnalysisDataContainer *cOutputCov  = mgr->CreateContainer(Form("Covariance%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,3,cOutputCov);
    AliAnalysisDataContainer *cOutputQA = mgr->CreateContainer(Form("QAContainer%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,4,cOutputQA); //For QA
    return task;
  };
  if(StageSwitch==10) {
    AliAnalysisDataContainer *effCont = mgr->CreateContainer(Form("EffList%s",l_ContName.Data()),AliEffFDContainer::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,1,effCont);
    return task;
  }

  return 0;
}
AliAnalysisTaskMeanPtV2Corr* AddTaskMeanPtV2Corr(TString name="name", Bool_t IsMC=kFALSE, TString stage="Efficiency",
                                                  TString efficiencyPath="", TString meanPtPath="", TString NUAPath="", TString subfix2="")
{
  return AddTaskMeanPtV2Corr(name,IsMC,stage,efficiencyPath,meanPtPath,NUAPath,"",subfix2);
}
