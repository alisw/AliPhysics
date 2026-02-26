
void AddTaskAliresonance(const char *suffix = "Resochrg")

{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskResoPU                             \n");
  printf("===================================================================================\n");

  TGrid::Connect("alien://");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  Int_t gCentMin = 0;
  Int_t gCentMax = 100;

  TString TaskResoPID;
  TaskResoPID.Form("gTaskResoCent%d_%d_%s", gCentMin, gCentMax, suffix);

  Aliresonance *task_Reso = new Aliresonance(TaskResoPID);

  //task_Reso->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
  //printf("\n =========> AddTaskReso::Info() Trigger = kINT7 \n");


  ///---> Now Pass data and containers to Analysis Object ----
  
  mgr->AddTask(task_Reso);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_Reso, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_Reso, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTaskReso() Configured properly <==================\n\n");

  
}
