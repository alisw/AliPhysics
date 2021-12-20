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
                                                  TString efficiencyPath, TString NUAPath, TString subfix1, TString subfix2 = "")
{
  int AnalysisStage = 0;
  printf("Stage switch name: %s\n",stage.Data());
  if(stage.Contains("weights")) AnalysisStage=1;
  if(stage.Contains("ptcorr")) AnalysisStage=2;
  if(AnalysisStage==0) return 0;
  TString l_ContName(subfix1);
  if(!subfix2.IsNull()) l_ContName+="_"+subfix2;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;

  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisTaskPtCorr* task = new AliAnalysisTaskPtCorr(name.Data(), IsMC, stage, l_ContName);
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
    printf("Inputs connected!\n");
    AliAnalysisDataContainer *cOutput = mgr->CreateContainer(Form("Output_%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    AliAnalysisDataContainer *cPtcorr = mgr->CreateContainer(Form("PtCorr_%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    printf("Output containers created!\n");
    mgr->ConnectOutput(task,1,cOutput);
    mgr->ConnectOutput(task,2,cPtcorr);
    return task;
  };

  return 0;
}

