
void AddTaskResonanceVsMultiplicity(Int_t fCentralityMin=0, Int_t fCentralityMax=90,/* TString sTrigger="kINT7"*/ Double_t fVzMax=10, Double_t fdcaxy=0.1, Double_t fdcaz=1, Double_t fchi2tpc=2.5, Double_t fchi2its=36, Double_t fnCrossedRows=70)
{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCMWPU                             \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TString TaskMeanpt;
  TaskMeanpt.Form("gTaskMeanpt%d_%d_%s", gCentMin, gCentMax, " ");
  gROOT->LoadMacro("AliAnalysisTaskResonanceVsMultiplicity.cxx++g");                                                           
  AliAnalysisTaskResonanceVsMultiplicity *task_Mpt = new AliAnalysisTaskResonanceVsMultiplicity(TaskMeanpt);

  /*
  ///-------> Analysis Object Created, now pass the arguments
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    task_Mpt->SelectCollisionCandidates(AliVEvent::kMB);
    printf("\n =========> AddTaskResonancevsMultiplicity::Info() Trigger = kMB  \n");
  }
  else if(sTrigger=="kSemiCentral" || sTrigger=="SemiCentral" || sTrigger=="semicentral"){
    task_Mpt->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kSemiCentral \n");
  }
  else if(sTrigger=="kCentral" || sTrigger=="Central" || sTrigger=="central"){
    task_Mpt->SelectCollisionCandidates(AliVEvent::kCentral);
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kCentral \n");
  }
  else{//if trigger==kINT7 or no trigger provided:
    task_Mpt->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kINT7 \n");
  }
  */
  ///swati: add event and track cuts
  
  ///Event cuts:
  /*
  task_Mpt->SetCentralityPercentileMin(fCentralityMin);
  task_Mpt->SetCentralityPercentileMax(fCentralityMax);
  task_Mpt->SetVzRangeMin(fVzMin);
  */
  task_Mpt->SetVzRangeMax(fVzMax);
       

  //Track cuts:
  task_Mpt->SetDCAXYRangeMax(fdcaxy);
  task_Mpt->SetDCAZRangeMax(fdcaz);
  task_Mpt->SetMaxChi2PerTPCClusterRange(fchi2tpc);
  task_Mpt->SetMaxChi2PerITSClusterRange(fchi2its);
  task_Mpt->SetMinNCrossedRowsTPCRange(fnCrossedRows);
  
  mgr->AddTask(task_Mpt);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_Mpt, 0, cinput);        


  AliAnalysisDataContainer  *cOutPut1;
  AliAnalysisDataContainer  *cOutPut2;
  AliAnalysisDataContainer  *cOutPut3;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s"," ");
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TTree::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  cOutPut2 = (AliAnalysisDataContainer *) mgr->CreateContainer("PhiSpectraPPb1",TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  cOutPut3 = (AliAnalysisDataContainer *) mgr->CreateContainer("QAPileupPlots",TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  
  mgr->ConnectOutput(task_Mpt, 1, cOutPut1);
  mgr->ConnectOutput(task_Mpt, 2, cOutPut2);
  mgr->ConnectOutput(task_Mpt, 3, cOutPut3);
  
 
  printf("\n\n ================> AddTaskRESONANCEVSMULTIPLICTY() Configured properly <==================\n\n",);

  //return task_Mpt;

}//Task Ends
