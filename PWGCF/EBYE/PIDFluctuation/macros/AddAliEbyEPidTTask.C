//=========================================================================//
//                                                                         //
//           Analysis AddTask for Particle Ratio Fluctuation Study         //
//              Author: Deepika Rathee  || Satyajit Jenara                 //
//                      drathee@cern.ch || sjena@cern.ch 
//                        Thu Jun 19 11:44:51 CEST 2014
//                                                                         //
//=========================================================================//

AliEbyEPidTTask *AddAliEbyEPidTTask(Bool_t isModeAOD    = 1,
				    Int_t AODfilterBit = 768, 
				    
				    Int_t pidtype       = 2, 
				    Int_t requestTofPid = 1,
				    Double_t nSigmaCut  = 3.,
				    Double_t lptfortof  = 0.5,
				    
				    Double_t ptl        = 0.5, 
				    Double_t pth        = 5.0, 
				    Double_t gEta       = 0.8,
				    
				    Double_t dcaxy     = 2.4,
				    Double_t dcaz      = 3.2,
				    
				    Double_t vz         = 10.,
				    
				    TString ctaskname = "D2011") {
  
 
  TString taskname = "";
  taskname += ctaskname;
  taskname += Form("%d",isMC);
  
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
  
  AliEbyEPidTTask *taskqa;
  taskqa = new AliEbyEPidTTask(taskname.Data());
  taskqa->SetAODtrackCutBit(AODfilterBit);
  taskqa->SetCentralityEstimator(centralityEstimator.Data());
  taskqa->SetAnalysisType(isMC);
  taskqa->SetHelperPID(help);
  mgr->AddTask(taskqa);
  
  AliAnalysisDataContainer *coutqa = mgr->CreateContainer(Form("%s_QA",taskname.Data()),TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s",basefilename.Data()));
  AliAnalysisDataContainer *coutt = mgr->CreateContainer("Event",TTree::Class(),AliAnalysisManager::kOutputContainer, Form("%s",basefilename.Data()));

  mgr->ConnectInput(taskqa, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskqa, 1, coutqa);
 mgr->ConnectOutput(taskqa, 2, coutt);
 
 
  return;
}
