#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisTaskTrackingEffPID.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

AliAnalysisTaskTrackingEffPID* AddTaskTrackingEffPID(TString suffix = "",
						     TString collSyst="pp",
						     bool useGeneratedKine=kTRUE,
						     TString cutObjFile = "",
						     TString cutObjNam = "",
						     Int_t filtBit=4) {

  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLFefficiencies", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskLFefficiencies", "This task requires an input event handler");
    return 0x0;
  }
  
  AliAnalysisTaskTrackingEffPID *taskeff = new AliAnalysisTaskTrackingEffPID();
  taskeff->SetUseGeneratedKine(useGeneratedKine);
  taskeff->SetCollisionSystem(collSyst);

  if(!cutObjFile.IsNull()){
    TFile *f=TFile::Open(cutObjFile.Data(),"READ");
    if(f){
      AliESDtrackCuts* esdc = (AliESDtrackCuts*)f->Get(cutObjNam.Data());
      if(esdc){
	taskeff->SetTrackCuts(esdc);
	taskeff->UseTrackCutObjectForAODTracks();
	taskeff->SetFilterBitCutForAODTracks(-1);
	printf("Use AliESDtrackCuts from file %s\n",cutObjFile.Data());
      }
    }
  }
  
  mgr->AddTask(taskeff);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":TrackEffPID";

  TString listname="listTrackEffPID";
  listname+=suffix.Data();

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(listname,
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   outputFileName.Data() );

  TString cutlistname="listTrackCuts";
  cutlistname+=suffix.Data();

  AliAnalysisDataContainer* ccuts = mgr->CreateContainer(cutlistname,
							 TList::Class(),
							 AliAnalysisManager::kParamContainer,
							 outputFileName.Data() );
  
  mgr->ConnectInput  (taskeff,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (taskeff,  1, coutput);
  mgr->ConnectOutput (taskeff,  2, ccuts);
  return taskeff;
}
