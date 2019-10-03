//=========================================================================//
//                                                                         //
//           Analysis AddTask for Particle Ratio Fluctuation Study         //
//              Author: Deepika Rathee  || Satyajit Jenara                 //
//                      drathee@cern.ch || sjena@cern.ch 
//                        Thu Jun 19 11:44:51 CEST 2014
//                                                                         //
//=========================================================================//

void AddAliEbyEPidTTaskMC(Double_t ptl=0.5, 
			  Int_t AODfilterBit = 768, 
			  Int_t pidtype = 2, 
			  Int_t requestTofPid = 1, 
			  Double_t nSigmaCut = 3., 
			  TString centralityEstimator = "V0M", 
			  TString ctaskname = "2011") {
  
  Bool_t isMC = 1; 
  
  TString taskname = "EbyECF_";
  taskname += ctaskname;
  taskname += "_";
  taskname += Form("%d",isMC);
  taskname += "_";
  taskname += centralityEstimator;
  taskname += "_";
  taskname += Form("%d",AODfilterBit);
  
 
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
  
  
  AliHelperPID* help = new AliHelperPID();
  help->SetNSigmaCut(nSigmaCut);
  help->SetPIDType(pidtype);        
  help->SetfRequestTOFPID(requestTofPid);
  help->SetfPtTOFPID(ptl);
  help->SetisMC(isMC);
  
  AliEbyEPidTTaskMC *taskqa;
  taskqa = new AliEbyEPidTTaskMC(taskname.Data());
  taskqa->SetAODtrackCutBit(AODfilterBit);
  taskqa->SetCentralityEstimator(centralityEstimator.Data());
  taskqa->SetAnalysisType(isMC);
  taskqa->SetHelperPID(help);
  mgr->AddTask(taskqa);
  
  AliAnalysisDataContainer *coutqa = mgr->CreateContainer(Form("%s_QA",taskname.Data()),TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s",basefilename.Data()));
  AliAnalysisDataContainer *coutt = mgr->CreateContainer("fEventTree",TTree::Class(),AliAnalysisManager::kOutputContainer, Form("%s",basefilename.Data()));

  mgr->ConnectInput(taskqa, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskqa, 1, coutqa);
 mgr->ConnectOutput(taskqa, 2, coutt);
 
 
  return;
}
