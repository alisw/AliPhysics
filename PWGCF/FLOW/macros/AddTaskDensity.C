#include "TGrid.h"
#include "TString.h"
class TNamed;
AliAnalysisTaskDensity* AddTaskDensity(TString name = "name", TString efficiencyFile = "", const char* suffix = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0; 
  if (!mgr->GetInputEventHandler()) return 0x0;

  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":DensityTask";     
  TString taskName = Form("%s_%s", name.Data(), suffix);

  Bool_t useEfficiency = kFALSE;
  if(!efficiencyFile.IsNull()) useEfficiency = kTRUE;

  // setup task
  printf("Setting task");
  AliAnalysisTaskDensity* task = new AliAnalysisTaskDensity(taskName.Data(), useEfficiency);   
  if(!task) return 0x0;
  task->SelectCollisionCandidates(AliVEvent::kINT7+AliVEvent::kCentral+AliVEvent::kSemiCentral);    
  mgr->AddTask(task);
  
  // input
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // output
  mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("PtSubEventSamples_%s", taskName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("QAAliEventCuts_%s", taskName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("QATrackCuts_%s", taskName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

  // get NUE files as input
  if(useEfficiency) {
    TObjArray* taskContainers = mgr->GetContainers();
    if(!taskContainers) { printf("Task containers does not exists!\n"); return NULL; }

    AliAnalysisDataContainer* efficiency = (AliAnalysisDataContainer*) taskContainers->FindObject("inputEfficiency");
    if(!efficiency){
      printf("Input file name: %s \n", efficiencyFile.Data());
      if(efficiencyFile.Contains("alien:")){ if(!gGrid) TGrid::Connect("alien:"); }; 

      TFile* efficiency_file = TFile::Open(efficiencyFile.Data(),"READ");
      if(!efficiency_file) { printf("Input file with efficiency not found!\n"); return NULL; }

      TList* efficiency_list = (TList*)efficiency_file->Get("EffAndFD");
      if(!efficiency_list) { printf("E-AddTask: Input list with efficiency not found!\n"); efficiency_file->ls(); return NULL; }

      AliAnalysisDataContainer* cInputEfficiency = mgr->CreateContainer("inputEfficiency", TList::Class(), AliAnalysisManager::kInputContainer);
      cInputEfficiency->SetData(efficiency_list);
      mgr->ConnectInput(task,1,cInputEfficiency);  
    }
  else {
    mgr->ConnectInput(task,1,efficiency);
  }
  }
  return task;
}
