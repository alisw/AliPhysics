
AliAnalysisTaskMeanptFluctMC_effcheck *AddTaskMeanptFluctMC_effcheck(Int_t fCentralityMin=0, Int_t fCentralityMax=90,/* TString sTrigger="kINT7"*/ Double_t fVzMax=10, Double_t fdcaxy=0.1, Double_t fdcaz=1, Double_t fchi2tpc=2.5, Double_t fchi2its=36, Double_t fnCrossedRows=70, TString OutFileName = "_default", Double_t fEta=0.8, Int_t fMCGenChoice=2, Float_t fPtMin=0.2, Float_t fPtMax=5.0)
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

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TString TaskMeanpt;
  TaskMeanpt.Form("gTaskMeanpt%d_%d_%s", gCentMin, gCentMax, " ");
                                                
  AliAnalysisTaskMeanptFluctMC_effcheck *task_Mpt = new AliAnalysisTaskMeanptFluctMC_effcheck(TaskMeanpt);
  task_Mpt->SetMCGeneratorChoice(fMCGenChoice);
  task_Mpt->SetPtMaxMin(fPtMax,fPtMin);
  

  // TString OutTreeName;
  // OutTreeName += "fTreeEvent";
  // OutTreeName += OutFileName;
  // task_Mpt->SetTreeName(OutTreeName);
 
  
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
  //cOutPut3 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName3,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  
  mgr->ConnectOutput(task_Mpt, 1, cOutPut1);
  mgr->ConnectOutput(task_Mpt, 2, cOutPut2);
  

  return task_Mpt;

}//Task Ends
