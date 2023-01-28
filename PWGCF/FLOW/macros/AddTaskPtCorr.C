#include "TString.h"
#include "TGrid.h"
#include "AliAnalysisPtCorr.h"
class AliAnalysisDataContainer;
class TNamed;
Bool_t ConnectToGrid() {
  if(!gGrid) TGrid::Connect("alien:");
  if(!gGrid) {printf("Task requires connection to grid, but it could not be established!\n"); return kFALSE; };
  return kTRUE;
}
AliAnalysisTaskPtCorr* AddTaskPtCorr(TString name, bool IsMC, bool isOnTheFly, unsigned int fl_eff, TString efficiencyPath, TString subfix1)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  TString l_ContName(subfix1);
  printf("----Initialising task %s ----\n",l_ContName.Data());
  printf("IsMC: %i\n",IsMC);
  printf("IsOnTheFly: %i\n",isOnTheFly);
  printf("Efficiency flags:\n");
  printf("No efficiency: %i\n",(fl_eff&PtCorrFlags::noeff)==PtCorrFlags::noeff);
  printf("Constant efficiency: %i\n",(fl_eff&PtCorrFlags::consteff)==PtCorrFlags::consteff);
  printf("Gaussian efficiency: %i\n",(fl_eff&PtCorrFlags::gausseff)==PtCorrFlags::gausseff);
  printf("Flat efficiency: %i\n",(fl_eff&PtCorrFlags::flateff)==PtCorrFlags::flateff);
  printf("Power efficiency: %i\n",(fl_eff&PtCorrFlags::powereff)==PtCorrFlags::powereff);
  printf("Real efficiency input: %i\n",(fl_eff&PtCorrFlags::realeffin)==PtCorrFlags::realeffin);
  if(!l_ContName.IsNull()) l_ContName.Prepend("_");
  AliAnalysisTaskPtCorr* task = new AliAnalysisTaskPtCorr(name.Data(), IsMC, isOnTheFly, fl_eff, l_ContName);
  if(!task)
    return 0x0;
  mgr->AddTask(task); // add your task to the manager
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0);
  //Full analysis
  TObjArray *AllContainers = mgr->GetContainers();
  Bool_t gridConnected=kFALSE;
  if(!(IsMC&&!((fl_eff&PtCorrFlags::realeffin)==PtCorrFlags::realeffin)) &&!isOnTheFly) {
    if(!AllContainers->FindObject("Efficiency")) {
      printf("Getting input...\n");
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
      printf("Inputs connected!\n");    
    } else { mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("Efficiency")); printf("Inputs already connected\n"); }
  };
  AliAnalysisDataContainer *cPtcorr = mgr->CreateContainer(Form("Correlations%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  AliAnalysisDataContainer *cQA = mgr->CreateContainer(Form("QA%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  mgr->ConnectOutput(task,1,cPtcorr);
  mgr->ConnectOutput(task,2,cQA);
  return task;
}

