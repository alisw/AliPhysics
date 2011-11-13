AliAnalysisTaskEMCALTriggerQA * AddTaskEMCALTriggerQA(TString outputFile = "", Int_t run = 0){
  
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALTriggerQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALTriggerQA", "This task requires an input event handler");
    return NULL;
  }
    
  AliAnalysisTaskEMCALTriggerQA * qatrigger = new AliAnalysisTaskEMCALTriggerQA("QATrigger");
  
  AliEMCALRecoUtils * reco = qatrigger->GetRecoUtils();
  reco->SwitchOnRejectExoticCluster();
  
  // Pass the bad channels, need to access run number
  TString fileName="$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root";
  AliOADBContainer *contBC=new AliOADBContainer("");
  contBC->InitFromFile((char*)fileName.Data(),"AliEMCALBadChannels"); 
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(run); 
  if(arrayBC){
    TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject("pass1");
    if(arrayBCpass){
      
      reco->SwitchOnBadChannelsRemoval();
      printf("trigger REMOVE bad cells \n");
      
      for (Int_t i=0; i<10; ++i) {
        TH2I *hbm = reco->GetEMCALChannelStatusMap(i);
        if (hbm)
          delete hbm;
        hbm=(TH2I*)arrayBCpass->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
        
        if (!hbm) {
          AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
          continue;
        }
        
        hbm->SetDirectory(0);
        reco->SetEMCALChannelStatusMap(i,hbm);
      }
    } else printf("trigger AliEMCALRecoUtils ---Do NOT remove bad channels 1\n");
  }  else  printf("trigger AliEMCALRecoUtils ---Do NOT remove bad channels 2\n");
    
  if(outputFile.Length()==0)outputFile = AliAnalysisManager::GetCommonFileName(); 

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput = 
  mgr->CreateContainer("EMCALQATrigger", TList::Class(), AliAnalysisManager::kOutputContainer,  Form("%s:EMCALQATrigger",outputFile.Data()));
  mgr->AddTask(qatrigger);
  mgr->ConnectInput  (qatrigger, 0, cinput1);
  mgr->ConnectOutput (qatrigger, 1, coutput);
  
  return qatrigger;
  
}
