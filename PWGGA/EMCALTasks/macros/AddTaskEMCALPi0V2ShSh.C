//astahlle
AliAnalysisTaskEMCALPi0V2ShSh *AddTaskEMCALPi0V2ShSh()
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  if (!mgr)
    {
      ::Error("AddTaskEMCALPi0V2ShSh", "No analysis manager to connect to.");
      return NULL;
    }  
  
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEventplane", "This task requires an input event handler");
      return NULL;
    }		
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE; 
  cout<<"AddTaskEMCALPi0V2ShSh - MC config is: "<<ismc<<endl;
  
  if (ismc) return 0;
  
  UInt_t physMB = 0;
  UInt_t physJet = 0;
  UInt_t physGam = 0;
  UInt_t physEMC = 0; 
  
  TString sGeomName = AliEMCALGeometry::GetDefaultGeometryName();
  
  AliAnalysisTaskEMCALPi0V2ShSh *task = new AliAnalysisTaskEMCALPi0V2ShSh("EMCALPi0ShowerShV2");
  
  task->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  
  if(!ismc)
    {		
      mgr->AddTask(task);
      
      AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("hist", TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s",AliAnalysisManager::GetCommonFileName()));
      
      mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(task, 1, coutput1);
      
    }
  return task;
}

