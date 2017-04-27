class AliAnalysisDataContainer;

AliAnalysisTaskHFEMultiplicity* AddTaskHFEMultiplicity(TString suffixName = "",
						       Int_t PhysSel = AliVEvent::kINT7,
						       Bool_t useTender   =  kTRUE,
						       Bool_t ClsTypeEMC  =  kTRUE,
						       Bool_t ClsTypeDCAL =  kTRUE
						       
						       
						       )
  
{ // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }
  
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":MyTask";     
  
  if(PhysSel == AliVEvent::kINT7){
  AliAnalysisTaskHFEMultiplicity* HFEtaskINT7 = new AliAnalysisTaskHFEMultiplicity("");
  mgr->AddTask(HFEtaskINT7); //HFEtask is my task
  HFEtaskINT7->SelectCollisionCandidates(AliVEvent::kINT7);
  HFEtaskINT7->SetTenderSwitch(useTender);
  HFEtaskINT7->SetClusterTypeEMC(ClsTypeEMC);
  HFEtaskINT7->SetClusterTypeDCAL(ClsTypeDCAL); 
  // Create containers for input/output
  TString finDirname         = "_INT7";
  TString outBasicname       = "EID";
  
  finDirname 	      +=   suffixName.Data();
  outBasicname      +=   finDirname.Data();
  
  
  
  mgr->ConnectInput(HFEtaskINT7,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(HFEtaskINT7,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  }
  
  if(PhysSel == AliVEvent::kEMCEGA){
    if(ClsTypeEMC){
      
      AliAnalysisTaskHFEMultiplicity* HFEtaskEG1 = new AliAnalysisTaskHFEMultiplicity("");
      mgr->AddTask(HFEtaskEG1);
      HFEtaskEG1->SelectCollisionCandidates(AliVEvent::kEMCEGA);
      HFEtaskEG1->SetEMCalTriggerEG1(kTRUE);
      HFEtaskEG1->SetTenderSwitch(useTender);
      HFEtaskEG1->SetClusterTypeEMC(ClsTypeEMC);
      HFEtaskEG1->SetClusterTypeDCAL(ClsTypeDCAL);
      
      
      TString finDirname         = "_TrigGAEG1";
      TString outBasicname       = "EID";
      
      finDirname 	      +=   suffixName.Data();
      outBasicname      +=   finDirname.Data();
    
      
      
      mgr->ConnectInput(HFEtaskEG1,0,mgr->GetCommonInputContainer());
      mgr->ConnectOutput(HFEtaskEG1,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
      // EMCal EGA EG2
    
      AliAnalysisTaskHFEMultiplicity* HFEtaskEG2 = new AliAnalysisTaskHFEMultiplicity("");
      mgr->AddTask(HFEtaskEG2);
      HFEtaskEG2->SelectCollisionCandidates(AliVEvent::kEMCEGA);
      HFEtaskEG2->SetEMCalTriggerEG2(kTRUE);
      HFEtaskEG2->SetTenderSwitch(useTender);
      HFEtaskEG2->SetClusterTypeEMC(ClsTypeEMC);
      HFEtaskEG2->SetClusterTypeDCAL(ClsTypeDCAL);
        

      TString finDirname         = "_TrigGAEG2";
      TString outBasicname       = "EID";
  
      finDirname 	      +=   suffixName.Data();
      outBasicname      +=   finDirname.Data();
 
  
  
      mgr->ConnectInput(HFEtaskEG2,0,mgr->GetCommonInputContainer());
      mgr->ConnectOutput(HFEtaskEG2,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
      
    }
    
    if(ClsTypeDCAL){
      // DCal EGA DG1
      AliAnalysisTaskHFEMultiplicity* HFEtaskDG1 = new AliAnalysisTaskHFEMultiplicity("");
      mgr->AddTask(HFEtaskDG1);
      HFEtaskDG1->SelectCollisionCandidates(AliVEvent::kEMCEGA);
      HFEtaskDG1->SetEMCalTriggerDG1(kTRUE);
      HFEtaskDG1->SetTenderSwitch(useTender);
      HFEtaskDG1->SetClusterTypeEMC(ClsTypeEMC);
      HFEtaskDG1->SetClusterTypeDCAL(ClsTypeDCAL);
      
      TString finDirname         = "_TrigGADG1";
      TString outBasicname       = "EID";
      
      finDirname 	      +=   suffixName.Data();
      outBasicname      +=   finDirname.Data();
      
      
      
      mgr->ConnectInput(HFEtaskDG1,0,mgr->GetCommonInputContainer());
      mgr->ConnectOutput(HFEtaskDG1,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
      // DCal EGA DG2
      AliAnalysisTaskHFEMultiplicity* HFEtaskDG2 = new AliAnalysisTaskHFEMultiplicity("");
      mgr->AddTask(HFEtaskDG2);
      HFEtaskDG2->SelectCollisionCandidates(AliVEvent::kEMCEGA);
      HFEtaskDG2->SetEMCalTriggerDG2(kTRUE);
      HFEtaskDG2->SetTenderSwitch(useTender);
      HFEtaskDG2->SetClusterTypeEMC(ClsTypeEMC);
      HFEtaskDG2->SetClusterTypeDCAL(ClsTypeDCAL);
    
      TString finDirname         = "_TrigGADG2";
      TString outBasicname       = "EID";
      
      finDirname 	      +=   suffixName.Data();
      outBasicname      +=   finDirname.Data();
      
    
      
      mgr->ConnectInput(HFEtaskDG2,0,mgr->GetCommonInputContainer());
      mgr->ConnectOutput(HFEtaskDG2,1,mgr->CreateContainer(outBasicname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
      
    }
  }
  
  
  
  return NULL;
}
