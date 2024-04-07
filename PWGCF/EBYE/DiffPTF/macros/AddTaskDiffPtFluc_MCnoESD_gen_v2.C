AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2 *AddTaskDiffPtFluc_MCnoESD_gen_v2( TString OutFileName = "_default", Int_t fMCGenChoice=2, Double_t fEta_leftCut=0.0, Double_t fEta_rightCut=0.0)
{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCMWPU                             \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":ResultsMC";         // This directory contains result histograms


  TString TaskMeanpt;
  TaskMeanpt.Form("gTaskMeanpt_%s", " ");
                                                
  AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2 *task_Mpt = new AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2(TaskMeanpt);

  task_Mpt->SetMCGeneratorChoice(fMCGenChoice);
  task_Mpt->SetEtaLeftCut(fEta_leftCut);
  task_Mpt->SetEtaRightCut(fEta_rightCut);

  
  mgr->AddTask(task_Mpt);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_Mpt, 0, cinput);        


  AliAnalysisDataContainer  *cOutPut1;
  AliAnalysisDataContainer  *cOutPut2;
  //AliAnalysisDataContainer  *cOutPut3;
  TString                  sMyOutName1;
  TString                  sMyOutName2;
  TString                  sMyOutName3;
  sMyOutName1 += "SimpleTask_tree";
  sMyOutName1 += OutFileName;
  sMyOutName2 += "Histogram_TrackVariables";
  sMyOutName2 += OutFileName;
  sMyOutName3 += "QAPileupPlots";
  sMyOutName3 += OutFileName;
  
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName1,TTree::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  cOutPut2 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName2,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  
  mgr->ConnectOutput(task_Mpt, 1, cOutPut1);
  mgr->ConnectOutput(task_Mpt, 2, cOutPut2);
  

  return task_Mpt;

}//Task Ends
