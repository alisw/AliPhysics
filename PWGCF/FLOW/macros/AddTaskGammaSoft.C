#include "TString.h"
#include "TGrid.h"
class AliAnalysisDataContainer;
class TNamed;
Bool_t ConnectToGrid() {
  if(!gGrid) TGrid::Connect("alien:");
  if(!gGrid) {printf("Task requires connection to grid, but it could not be established!\n"); return kFALSE; };
  return kTRUE;
}
AliAnalysisTaskGammaSoft* AddTaskGammaSoft(TString name, Bool_t IsMC, TString efficiencyPath, TString NUAPath, TString subfix)
{
  printf("***********************************\n Adding Gamma soft task to grid\n***********************************\n");
  TString l_ContName(subfix);
  l_ContName.Prepend("_");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  // if(IsMC) {
  //   if(!mgr->GetMCtruthEventHandler()) {
  //     Error("AddTaskGammaSoft","Could not get MC truth handler");
  //     return NULL;
  //   };
  //   AliMCEventHandler *handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  //   handler->SetReadTR(kTRUE);
  // };
  AliAnalysisTaskGammaSoft* task = new AliAnalysisTaskGammaSoft(name.Data(), IsMC, l_ContName);
  if(!task)
    return 0x0;
  mgr->AddTask(task); // add your task to the manager
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0);
  TObjArray *AllContainers = mgr->GetContainers();
  // AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer(Form("MPTProfileList%s",l_ContName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  Bool_t gridConnected=kFALSE;
  if(!IsMC) {
    if(!AllContainers->FindObject("Weights")) {
      if(NUAPath.IsNull()) { printf("Weight path not provided!\n"); return 0; };
      if(NUAPath.Contains("alien:")) if(!ConnectToGrid()) return 0;
      TFile *tfWeights = TFile::Open(NUAPath.Data());
      TList *fList = (TList*)tfWeights->Get("WeightList");
      AliAnalysisDataContainer *cWeights = mgr->CreateContainer("Weights",TList::Class(), AliAnalysisManager::kInputContainer);
      cWeights->SetData(fList);
      mgr->ConnectInput(task,1,cWeights);
    } else mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("Weights"));
    if(!AllContainers->FindObject("Efficiency")) {
      if(efficiencyPath.IsNull()) { printf("Efficiency path not provided!\n"); return 0; };
      if(efficiencyPath.Contains("alien:")) if(!ConnectToGrid()) return 0;
      TFile *tfEff = TFile::Open(efficiencyPath.Data());
      if(!tfEff) { printf("Could not open efficiency file\n"); return 0; };
      if(tfEff->IsZombie()) { printf("Efficiency file is a zombie\n"); return 0; };
      TList *fList = (TList*)tfEff->Get("EffAndFD");
      if(!fList) { printf("Could not fetch the efficiency list!\n"); return 0; };
      AliAnalysisDataContainer *cEff = mgr->CreateContainer("Efficiency",TList::Class(), AliAnalysisManager::kInputContainer);
      cEff->SetData(fList);
      mgr->ConnectInput(task,2,cEff);
    } else mgr->ConnectInput(task,2,(AliAnalysisDataContainer*)AllContainers->FindObject("Efficiency"));
  };
  AliAnalysisDataContainer *cOutputMPT = mgr->CreateContainer(Form("MPTDiff%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  mgr->ConnectOutput(task,1,cOutputMPT);
  AliAnalysisDataContainer *cOutputFC  = mgr->CreateContainer(Form("FlowCont%s",l_ContName.Data()),AliGFWFlowContainer::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  mgr->ConnectOutput(task,2,cOutputFC);
  AliAnalysisDataContainer *cOutputCov  = mgr->CreateContainer(Form("Covariance%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  mgr->ConnectOutput(task,3,cOutputCov);
  AliAnalysisDataContainer *cOutputQA = mgr->CreateContainer(Form("QAContainer%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  mgr->ConnectOutput(task,4,cOutputQA); //For QA
  printf("***********************************\n Input and output connected \n***********************************\n");
  return task;
  return 0;
}

