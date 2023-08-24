
void AddTaskfn(const char *suffix = "CMWchrg")

{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCMWPU                             \n");
  printf("===================================================================================\n");

  TGrid::Connect("alien://");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  Int_t gCentMin = 0;
  Int_t gCentMax = 100;

  TString TaskCMWPID;
  TaskCMWPID.Form("gTaskCMWCent%d_%d_%s", gCentMin, gCentMax, suffix);

  AliAnalysisTaskfn *task_CMW = new AliAnalysisTaskfn(TaskCMWPID);

  //task_CMW->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
  //printf("\n =========> AddTaskCMW::Info() Trigger = kINT7 \n");


  ///---> Now Pass data and containers to Analysis Object ----
  
  mgr->AddTask(task_CMW);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CMW, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_CMW, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTaskCMW() Configured properly <==================\n\n");

  
}
