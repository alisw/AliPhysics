AliAnalysisTaskCaloConv * AddTaskCaloConv(){
  //Macro to add class CaloConv (conversion+calorimeters pi0 analysis) to train
  //Argument is the path to the PHOS recalibration parameters (file with OCDB entry)
  //Default path to the file with unit recalibration == no recalibnration
  //If file does not exist, no recalibration histograms will be filled

  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCaloConv", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCaloConv", "This task requires an input event handler");
    return NULL;
  }

  // Add task
  AliAnalysisTaskCaloConv *task = new AliAnalysisTaskCaloConv("CaloConv");
  mgr->AddTask(task);

  TDirectory* saveDir = gDirectory;
  TGrid a;
  if(a.IsConnected()){
    TFile *fBadMap = TFile::Open("alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10b.root") ;
    if(fBadMap->IsOpen()){
      if (saveDir) saveDir->cd(); else gROOT->cd();

      printf("Adding PHOS and EMCAL bad maps \n") ;
      char key[55] ;
      for(Int_t mod=1;mod<4; mod++){
        sprintf(key,"PHOS_BadMap_mod%d",mod) ;
        TH2I * h = (TH2I*)fBadMap->Get(key) ;
        if(h)
          task->SetPHOSBadMap(mod,h) ;
      }
      for(Int_t sm=0; sm<5; sm++){
        sprintf(key,"EMCAL_BadMap_mod%d",sm) ;
        TH2I * h = (TH2I*)fBadMap->Get(key) ;
        if(h)
          task->SetEMCALBadMap(mod,h) ;
      }
      fBadMap->Close() ;
    }
  }
  else{
    printf("Can not open Bad Map file \n") ;
  }


  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("CaloConv", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PWGGA_CaloConv",outputfile.Data()));
  mgr->ConnectOutput(task, 1, coutput);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("CFCaloConv", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PWGGA_CFCaloConv",outputfile.Data()));
  mgr->ConnectOutput(task, 2, coutput2);

  return task ;

}
