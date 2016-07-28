AliAnalysisTaskSEF01710fromAODtracks *AddTaskF01710fromAODtracks(
    Bool_t theMCon=kFALSE,
    Bool_t writeVariableTree=kFALSE,
    Int_t nTour=0
    )

{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLc2V0YW", "No analysis manager to connect to.");
    return NULL;
  }  

  //CREATE THE TASK
  printf("CREATE TASK\n");
  AliAnalysisTaskSEF01710fromAODtracks *task = new AliAnalysisTaskSEF01710fromAODtracks("AliAnalysisTaskSEF01710fromAODtracks",writeVariableTree);
  task->SetMC(theMCon);
  task->SetDebugLevel(1);

  task->SetEventMixingWithPools();
  //task->SetEventMixingOff();

  Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
  Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
  task->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);

  Double_t cent_mult_binlimitspp[] = { 0,100};
  Int_t cent_mult_bin_numbpp = sizeof(cent_mult_binlimitspp)/sizeof(Double_t) - 1;
  task->SetPoolCentBinLimits(cent_mult_bin_numbpp,cent_mult_binlimitspp);
  task->SetNumberOfEventsForMixing(10);//pp

  mgr->AddTask(task);

  // Create and connect containers for input/output  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWGLF_RES_F01710_";
  outputfile += nTour;

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // ----- output data -----
  AliAnalysisDataContainer *coutput1   = mgr->CreateContainer(Form("f0hist%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutputL2 = mgr->CreateContainer(Form("F0variables%1d",nTour),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); 
  mgr->ConnectOutput(task,2,coutputL2);
  AliAnalysisDataContainer *coutputL3 = mgr->CreateContainer(Form("F0All%1d",nTour),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,3,coutputL3);

  return task;

}
