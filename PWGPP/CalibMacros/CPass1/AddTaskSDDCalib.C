AliAnalysisTaskITSAlignQA *AddTaskSDDCalib(Int_t nrun=0, Bool_t writeITSTP=kFALSE, Bool_t useTPCcrv=kFALSE) 
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskITSAlign", "No analysis manager to connect to.");
    return NULL;
  }   
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskITSAlign", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); 
  if(type.Contains("AOD")){
    ::Error("AddTaskITSAlign", "This task requires to run on ESD");
    return NULL;
  }
  
  // Create and configure the task
  AliAnalysisTaskITSAlignQA *taskali = new AliAnalysisTaskITSAlignQA();
  //  taskali->SelectCollisionCandidates();
  if (nrun>0) taskali->SetOCDBInfo(nrun,Form("raw://")) ; 
  taskali->SetLoadGeometryFromOCDB(kFALSE);
  mgr->AddTask(taskali);
  //  
  taskali->SetUseVertex(kTRUE);
  taskali->SetUseVertexForZOnly(kFALSE);
  taskali->SetDoSPDResiduals(kFALSE);
  taskali->SetDoSDDResiduals(kFALSE);
  taskali->SetDoSSDResiduals(kFALSE);
  taskali->SetDoSDDDriftTime(kFALSE);
  taskali->SetMinMaxMult(20.,1070.);
  taskali->SetUseTPCMomentum(useTPCcrv);
  //
  taskali->SetDoFillTPTree(writeITSTP);
  //
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistSDDCalib",
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName );
  
  mgr->ConnectInput(taskali, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskali, 1, coutput1);
  //
  if (writeITSTP) { // if the TPtree to be written, create separate container
    AliAnalysisDataContainer *coutputTP = mgr->CreateContainer("ITSSumTP",
							       TTree::Class(),
							       AliAnalysisManager::kOutputContainer,
							       "AODtpITS.root" );
    mgr->ConnectOutput(taskali, 2, coutputTP);
  }
  //

  return taskali;
}
