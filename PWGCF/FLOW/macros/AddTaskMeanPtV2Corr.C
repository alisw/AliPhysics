#include "TString.h"
#include "TGrid.h"
class AliAnalysisDataContainer;
class TNamed;
Bool_t ConnectToGrid() {
  if(!gGrid) TGrid::Connect("alien:");
  if(!gGrid) {printf("Task requires connection to grid, but it could not be established!\n"); return kFALSE; };
  return kTRUE;
}
AliAnalysisTaskMeanPtV2Corr* AddTaskMeanPtV2Corr(TString name = "name", Bool_t IsMC=kFALSE, TString stage = "weights",
                                                  TString efficiencyPath = "", TString meanPtPath="", TString NUAPath="", TString subfix="")
{
  Int_t StageSwitch = 0;
  if(stage.Contains("weights")) StageSwitch=1;
  if(stage.Contains("meanpt")) StageSwitch=2;
  if(stage.Contains("full")) StageSwitch=3;
  if(stage.Contains("ALICEMpt")) StageSwitch=4;
  if(stage.Contains("ALICECov")) StageSwitch=5;
  if(stage.Contains("FBSpectra")) StageSwitch=6;
  if(stage.Contains("Efficiency")) StageSwitch=7;
  if(StageSwitch==0) return 0;

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
  AliAnalysisTaskMeanPtV2Corr* task = new AliAnalysisTaskMeanPtV2Corr(name.Data(), IsMC, stage, subfix);
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
    if(!IsMC) {
      TObjArray *AllContainers = mgr->GetContainers();
      if(!AllContainers->FindObject("Efficiency")) {
        if(efficiencyPath.IsNull()) { printf("Efficiency path not provided!\n"); return 0; };
        if(efficiencyPath.Contains("alien:")) if(!ConnectToGrid()) return 0;//TGrid::Connect("alien:");
        TFile *tfWeights = TFile::Open(efficiencyPath.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MergedWeights.root"
        if(!tfWeights) { printf("Could not open efficiency file\n"); return 0; };
        if(tfWeights->IsZombie()) { printf("Efficiency file is a zombie\n"); return 0; };
        TList *fList = (TList*)tfWeights->Get("EffAndFD");
        if(!fList) { printf("Could not fetch the efficiency list!\n"); return 0; };
        AliAnalysisDataContainer *cEff = mgr->CreateContainer("Efficiency",TList::Class(), AliAnalysisManager::kInputContainer);
        cEff->SetData(fList);
        mgr->ConnectInput(task,1,cEff);
      } else mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("Efficiency"));
    };
    TString l_ContName=subfix.IsNull()?"":("_" + subfix);
    AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer(Form("MPTProfileList%s",l_ContName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputMPT);
    return task;
  };
  //Full
  if(StageSwitch==3) {
    TObjArray *AllContainers = mgr->GetContainers();
    TString l_ContName=subfix.IsNull()?"":("_" + subfix);
    // AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer(Form("MPTProfileList%s",l_ContName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    Bool_t gridConnected=kFALSE;
    if(!AllContainers->FindObject(Form("MPTProfileList%s",l_ContName.Data()))) {
      if(meanPtPath.IsNull()) { printf("Mean pT path not provided!\n"); return 0; };
      if(meanPtPath.Contains("alien:")) if(!ConnectToGrid()) return 0;//{ TGrid::Connect("alien:"); gridConnected = kTRUE; }; //Only connect if not connected yet
      TFile *tfMPT = TFile::Open(meanPtPath.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MeanPts_05_20.root"
      TList *fMPTList = (TList*)tfMPT->Get(Form("MPTProfileList%s",l_ContName.Data()));
      if(!fMPTList) { printf("fMPT list from file not fetched! I was looking for MPTProfileList%s; contents: \n",l_ContName.Data()); tfMPT->ls(); return 0; };
      AliAnalysisDataContainer *cInMPT = mgr->CreateContainer(Form("MPTProfileList%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kInputContainer);
      cInMPT->SetData(fMPTList);
      mgr->ConnectInput(task,1,cInMPT);
    } else mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject(Form("MPTProfileList%s",l_ContName.Data())));
    if(!IsMC) {
      if(!AllContainers->FindObject("Efficiency")) {
        if(efficiencyPath.IsNull()) { printf("Efficiency path not provided!\n"); return 0; };
        if(efficiencyPath.Contains("alien:")) if(!ConnectToGrid()) return 0;//{ TGrid::Connect("alien:"); gridConnected = kTRUE; };
        TFile *tfEff = TFile::Open(efficiencyPath.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MergedWeights.root"
        if(!tfEff) { printf("Could not open efficiency file\n"); return 0; };
        if(tfEff->IsZombie()) { printf("Efficiency file is a zombie\n"); return 0; };
        TList *fList = (TList*)tfEff->Get("EffAndFD");
        if(!fList) { printf("Could not fetch the efficiency list!\n"); return 0; };
        AliAnalysisDataContainer *cEff = mgr->CreateContainer("Efficiency",TList::Class(), AliAnalysisManager::kInputContainer);
        cEff->SetData(fList);
        mgr->ConnectInput(task,2,cEff);
      } else mgr->ConnectInput(task,2,(AliAnalysisDataContainer*)AllContainers->FindObject("Efficiency"));
      if(!AllContainers->FindObject("Weights")) {
        if(NUAPath.IsNull()) { printf("Weight path not provided!\n"); return 0; };
        if(NUAPath.Contains("alien:")) if(!ConnectToGrid()) return 0;//{ TGrid::Connect("alien:"); gridConnected = kTRUE; };
        TFile *tfWeights = TFile::Open(NUAPath.Data());
        TList *fList = (TList*)tfWeights->Get("WeightList");
        AliAnalysisDataContainer *cWeights = mgr->CreateContainer("Weights",TList::Class(), AliAnalysisManager::kInputContainer);
        cWeights->SetData(fList);
        mgr->ConnectInput(task,3,cWeights);
      } else mgr->ConnectInput(task,3,(AliAnalysisDataContainer*)AllContainers->FindObject("Weights"));
    };

    AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer(Form("MPTDiff%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputMPT);
    AliAnalysisDataContainer *cOutputFC  = mgr->CreateContainer(Form("FlowCont%s",l_ContName.Data()),AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,2,cOutputFC);
    AliAnalysisDataContainer *cOutputCov  = mgr->CreateContainer(Form("Covariance%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,3,cOutputCov);
    AliAnalysisDataContainer *cOutputV2dPt  = mgr->CreateContainer(Form("V2vsdPt%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,4,cOutputV2dPt);
    return task;
  };
  if(StageSwitch==4) {
    AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer("MPTProfileList",TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputMPT);
    return task;
  }
  if(StageSwitch==5) {
    TObjArray *AllContainers = mgr->GetContainers();
    if(meanPtPath.IsNull()) { printf("Input weights not set\n"); return 0; };
    if(!AllContainers->FindObject("InputMeanPt")) {
      if(meanPtPath.Contains("alien:")) if(!ConnectToGrid()) return 0;//TGrid::Connect("alien:"); //Only connect if not connected yet
      TFile *tfMPT = TFile::Open(meanPtPath.Data()); //"alien:///alice/cern.ch/user/v/vvislavi/MeanPts/MeanPts_05_20.root"
      TList *fMPTList = (TList*)tfMPT->Get("MPTProfileList");
      if(!fMPTList) { printf("fMPT list from file not fetcehd!"); return 0; };
      AliAnalysisDataContainer *cInMPT = mgr->CreateContainer("InputMeanPt",TList::Class(), AliAnalysisManager::kInputContainer);
      cInMPT->SetData(fMPTList);
      mgr->ConnectInput(task,1,cInMPT);
    } else mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("InputMeanPt"));
    AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer("OutputList",TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputMPT);
  }
  if(StageSwitch==6) { //Producing Pt spectra with filter bit
    AliAnalysisDataContainer *cOutputSpectra = mgr->CreateContainer("PtSpectra",TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputSpectra);
    return task;
  }
  if(StageSwitch==7) { //Producing Pt spectra with filter bit
    TString l_ContName=subfix.IsNull()?"":("_" + subfix);
    AliAnalysisDataContainer *spectraCont = mgr->CreateContainer(Form("SpectraList%s",l_ContName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,1,spectraCont);
    return task;
  }

  return 0;
}
