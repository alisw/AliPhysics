#include "AliAnalysisManager.h"
#include "TString.h"
#include "TGrid.h"
#include "AliAnalysisTaskPtCorr.h"
#include "AliAnalysisDataContainer.h"

Bool_t ConnectToGrid();
AliAnalysisTaskPtCorr* AddTaskPtCorr(TString name, bool IsMC, TString efficiencyPath = "", TString centcalPath = "", TString subfix1 = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;
  TString fileName = AliAnalysisManager::GetCommonFileName();
  TString l_ContName(subfix1);

  if(!l_ContName.IsNull()) l_ContName.Prepend("_");
  AliAnalysisTaskPtCorr* task = new AliAnalysisTaskPtCorr(name.Data(), IsMC, l_ContName);
  if(!task)
    return 0x0;
  mgr->AddTask(task); // add your task to the manager
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,0,cInput0);
  //Full analysis
  TObjArray *AllContainers = mgr->GetContainers();
  Bool_t gridConnected=kFALSE;
  if(!IsMC) {
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
  if(IsMC){
    if(!AllContainers->FindObject("CentCalibration")) {
      if(centcalPath.Contains("alien:")) if(!ConnectToGrid()) return 0;
      TFile *tfcentcal = TFile::Open(centcalPath.Data());
      if(!tfcentcal) { printf("Could not open centrality calibration file\n"); return 0; };
      if(tfcentcal->IsZombie()) { printf("Centrality calibration file is a zombie\n"); return 0; };
      TH1 *fcentcal = (TH1*)tfcentcal->Get("centcal");
      if(!fcentcal) { printf("Could not fetch the centrality calibration histogram!\n"); return 0; };
      AliAnalysisDataContainer *cCentCal = mgr->CreateContainer("CentCalibration",TH1::Class(), AliAnalysisManager::kInputContainer);
      cCentCal->SetData(fcentcal);
      mgr->ConnectInput(task,1,cCentCal);
      printf("Centrality calibration input connected!\n");
    } else { mgr->ConnectInput(task,1,(AliAnalysisDataContainer*)AllContainers->FindObject("CentCalibration")); printf("Inputs already connected\n"); }
  }
  AliAnalysisDataContainer *cPtcorr = mgr->CreateContainer(Form("Correlations%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  AliAnalysisDataContainer *cQA = mgr->CreateContainer(Form("QA%s",l_ContName.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  mgr->ConnectOutput(task,1,cPtcorr);
  mgr->ConnectOutput(task,2,cQA);
  return task;
}
Bool_t ConnectToGrid() {
  if(!gGrid) TGrid::Connect("alien:");
  if(!gGrid) {printf("Task requires connection to grid, but it could not be established!\n"); return kFALSE; };
  return kTRUE;
}
