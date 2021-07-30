AliAnalysisTask *AddTaskLeadingMC(){
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLeadingMC", "No analysis manager found.");
    return 0;
  }

  // currently don't accept AOD input
  if (!mgr->GetInputEventHandler()->InheritsFrom(AliESDInputHandler::Class())) {
    Error("AddTaskLeadingMC","The analysis tender only works with ESD input!");
    return 0;
  }

  //========= Add tender to the ANALYSIS manager and set default storage =====
  AliTaskLeadingMC *task = new AliTaskLeadingMC();

  task->SetEtaThreshold(8.0);                          // above this threashold candidate for leading effect are kept
  task->SetEnergyThreshold(0.0);                       // above this threashold candidate for	leading	effect are kept
  task->SetEtaBarrel(0.5);                             // barrel quanity computed in this eta range

  task->SetZDCPGeo(9., 27., -7., 11., 10000., 13000.); // hit kept in this volume for ZDCP
  task->SetZDCNGeo(-4, 4., -4., 4., 11000., 12500.);   // hit kept in this volume for ZDCN

  mgr->AddTask(task);

  //Attach input to my tasks
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  AliAnalysisDataContainer *cOutputT= mgr->CreateContainer("treeZDCmc",TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectOutput(task, 0, cOutputT);

  return task;
}
