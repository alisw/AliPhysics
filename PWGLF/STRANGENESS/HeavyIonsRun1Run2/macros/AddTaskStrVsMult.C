AliAnalysisTaskStrVsMult *AddTaskStrVsMult(UInt_t triggerMask = AliVEvent::kINT7)
{
  // analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { 
    ::Error("AddTaskStrangenessVsMultiplicity", "No analysis manager to connect to."); 
    return NULL; 
  }
  if (!mgr->GetInputEventHandler()) { 
    ::Error("AddTaskStrangenessVsMultiplicity", "This task requires an input event handler"); 
    return NULL; 
  }

  // Create the task and add it to the manager
  AliAnalysisTaskStrVsMult *mytask = new AliAnalysisTaskStrVsMult("StrVsMult_Task");
  mytask->SelectCollisionCandidates(triggerMask);
  mgr->AddTask(mytask);

  // output file name
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGLF_StrVsMult";
  printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

  //output containers
  AliAnalysisDataContainer *coutput_0, *coutput_1, *coutput_2, *coutput_3, *coutput_4, *coutput_5, *coutput_6, *coutput_7;
  coutput_0 = mgr->CreateContainer("chists_eve", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
  coutput_1 = mgr->CreateContainer("chists_K0S", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
  coutput_2 = mgr->CreateContainer("chists_Lam", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
  coutput_3 = mgr->CreateContainer("chists_ALam",TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
  coutput_4 = mgr->CreateContainer("chists_Xim", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
  coutput_5 = mgr->CreateContainer("chists_Xip", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
  coutput_6 = mgr->CreateContainer("chists_Omm", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );
  coutput_7 = mgr->CreateContainer("chists_Omp", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName );

  //connecting input and output
  mgr->ConnectInput (mytask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(mytask, 1, coutput_0);
  mgr->ConnectOutput(mytask, 2, coutput_1);
  mgr->ConnectOutput(mytask, 3, coutput_2);
  mgr->ConnectOutput(mytask, 4, coutput_3);
  mgr->ConnectOutput(mytask, 5, coutput_4);
  mgr->ConnectOutput(mytask, 6, coutput_5);
  mgr->ConnectOutput(mytask, 7, coutput_6);
  mgr->ConnectOutput(mytask, 8, coutput_7);

  return mytask;

}
