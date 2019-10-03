//=========================================================================//
//                                                                         //
//          c AliEbyE OnFLy QA Tasks for Charge and PID  V1.0              //
//              Author: Satyajit Jena || Deepika Rathee                    //
//                      sjena@cern.ch || drathee@cern.ch                   //
//                                                                         //
//=========================================================================//


void AddAliEbyEPidQATaskOnFlyKine( Double_t etacut=0.5,Double_t ptcut=20.,Double_t vz = 30.) {
   
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); 

  TString basefilename = AliAnalysisManager::GetCommonFileName();
    
  AliEbyEPidQATaskOnFlyKine *taskqa = new AliEbyEPidQATaskOnFlyKine("QA");
  taskqa->SetKinematicCut(etacut,ptcut,vz);
  AliAnalysisDataContainer *couttqa = mgr->CreateContainer("QA",TList::Class(), AliAnalysisManager::kOutputContainer,
							   Form("%s",basefilename.Data()));
  mgr->ConnectInput(taskqa, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskqa, 1, couttqa);
  
  return;
}
