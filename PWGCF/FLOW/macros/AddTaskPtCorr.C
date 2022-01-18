#include "TString.h"
#include "TGrid.h"
class AliAnalysisDataContainer;
class TNamed;
Bool_t ConnectToGrid() {
  if(!gGrid) TGrid::Connect("alien:");
  if(!gGrid) {printf("Task requires connection to grid, but it could not be established!\n"); return kFALSE; };
  return kTRUE;
}
AliAnalysisTaskPtCorr* AddTaskPtCorr(TString name, bool IsMC, TString stage,
                                                  TString efficiencyPath, TString NUAPath, TString meanPtPath, TString subfix1, TString subfix2 = "")
{
  int AnalysisStage = 0;
  printf("Stage switch name: %s\n",stage.Data());
  if(stage.Contains("weights")) AnalysisStage=1;
  if(stage.Contains("ptcorr")) AnalysisStage=2;
  if(stage.Contains("meanpt")) AnalysisStage=3;
  if(AnalysisStage==0) { printf("Analysis stage not set\n"); return 0; }
  TString l_ContName(subfix1);
  if(!l_ContName.IsNull()) l_ContName.Prepend("_");
  if(!subfix2.IsNull()) l_ContName+="_"+subfix2;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  bool dodynamics = (meanPtPath.IsNull())?false:true;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskPtCorr* task = new AliAnalysisTaskPtCorr(name.Data(), IsMC, dodynamics, stage, l_ContName);
  if(!task)
    return 0x0;
  mgr->AddTask(task); // add your task to the manager
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0);
  //Producing weights
  if(AnalysisStage==1) {
    AliAnalysisDataContainer *weightCont = mgr->CreateContainer(Form("WeightList%s",l_ContName.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    mgr->ConnectOutput(task,1,weightCont);
    return task;
  }
  //Full analysis
  printf("Getting input...\n");
  if(AnalysisStage==3) {
    TObjArray *AllContainers = mgr->GetContainers();
    Bool_t gridConnected=kFALSE;
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
        mgr->ConnectInput(task,1,cEff);
      } else mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("Efficiency"));
      printf("Efficiencies connected!\n");
    }
    AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer(Form("MPTProfileList%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(task,1,cOutputMPT);
    printf("Output connected!\n");
    return task;
  }
  if(AnalysisStage==2) {
    TObjArray *AllContainers = mgr->GetContainers();
    Bool_t gridConnected=kFALSE;
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
        mgr->ConnectInput(task,1,cEff);
      } else mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("Efficiency"));
      printf("Efficiencies connected!\n");      
      if(!AllContainers->FindObject("Weights")) {
        if(NUAPath.IsNull()) { printf("Weight path not provided!\n"); return 0; };
        if(NUAPath.Contains("alien:")) if(!ConnectToGrid()) { printf("Cannot connect to grid\n"); return 0; }//{ TGrid::Connect("alien:"); gridConnected = kTRUE; };
        TFile *tfWeights = TFile::Open(NUAPath.Data());
        if(!tfWeights) { printf("Could not open %s\n",NUAPath.Data()); return 0; }
        TList* fList = (TList*)tfWeights->Get("WeightList");
        AliAnalysisDataContainer *cWeights = mgr->CreateContainer("Weights",TList::Class(), AliAnalysisManager::kInputContainer);
        cWeights->SetData(fList);
        mgr->ConnectInput(task,2,cWeights);
      } else mgr->ConnectInput(task,2,(AliAnalysisDataContainer*)AllContainers->FindObject("Weights"));
      printf("Acceptance corrections connected!\n");
    };
    if(!dodynamics) 
    { 
      printf("Input mean pt not set! Continuing without dynamic correlations! Make sure this is intended!\n"); 
      task->DoDynamicCorrelations(false); }
    else 
    {
      if(!AllContainers->FindObject("InputMeanPt")) {
        if(meanPtPath.Contains("alien:")) if(!ConnectToGrid()) return 0;
        TFile *tfMPT = TFile::Open(meanPtPath.Data()); 
        if(!tfMPT) { printf("File %s not found!\n",meanPtPath.Data()); return 0;}
        TList *fMPTList = (TList*)tfMPT->Get(Form("MPTProfileList%s",l_ContName.Data()));
        if(!fMPTList) { printf("fMPT list MPTProfileList%s from file %s not fetched!\n",l_ContName.Data(),meanPtPath.Data()); return 0; };
        AliAnalysisDataContainer *cInMPT = mgr->CreateContainer("InputMeanPt",TList::Class(), AliAnalysisManager::kInputContainer);
        cInMPT->SetData(fMPTList);
        (IsMC)?mgr->ConnectInput(task,1,cInMPT):mgr->ConnectInput(task,3,cInMPT);
      } else (IsMC)?mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("InputMeanPt")):mgr->ConnectInput(task,3,(AliAnalysisDataContainer*)AllContainers->FindObject("InputMeanPt"));
      printf("Mean pt connected!\n");
    }
    printf("Inputs connected!\n");
    AliAnalysisDataContainer *cPtcorr = mgr->CreateContainer(Form("Correlations%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    AliAnalysisDataContainer *cdyncorr = mgr->CreateContainer(Form("DynamicCorrelations%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    printf("Output containers created!\n");
    mgr->ConnectOutput(task,1,cPtcorr);
    mgr->ConnectOutput(task,2,cdyncorr);
    return task;
  };

  return 0;
}

