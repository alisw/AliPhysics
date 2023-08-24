
void AddTaskCMXQA(Int_t gFilterBit = 96, Double_t fPtMin=0.2, Double_t fPtMax=5.0, Double_t fEtaMin=-0.8, Double_t fEtaMax=0.8,  Int_t gNclustTPC=70, TString sCentEstimator="V0M", Double_t fCentralityMin=0., Double_t fCentralityMax=90., Float_t fVzMin = -10.0, Float_t fVzMax = 10.0, TString sTrigger="kINT7", const char *suffix = "")
{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCMXQA                          \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TString TaskCMEV0PID;
  TaskCMEV0PID.Form("SimpleTaskCent%d_%d_%s", gCentMin, gCentMax, suffix);
  //cout<<"Add taskname = "<<TaskCMEV0PID.Data()<<endl;

  AliAnalysisTaskCMXQA *task_CME = new AliAnalysisTaskCMXQA(TaskCMEV0PID);

  ///-------> Analysis Object Created, now pass the arguments
  
  task_CME->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    task_CME->SelectCollisionCandidates(AliVEvent::kMB);
  }
  ///Event cuts:
  task_CME->SetCentralityPercentileMin(fCentralityMin);
  task_CME->SetCentralityPercentileMax(fCentralityMax);
  task_CME->SetVzRangeMin(fVzMin);
  task_CME->SetVzRangeMax(fVzMax);
  if(sCentEstimator=="V0" || sCentEstimator=="V0M"){ 
    task_CME->SetCentralityEstimator("V0M");    
  }
  else{
    task_CME->SetCentralityEstimator(sCentEstimator);          //  use the Estimator provided in AddTask.
  }
     

  //Track cuts:
  task_CME->SetTrackCutNclusterMin(gNclustTPC);
  task_CME->SetFilterBit(gFilterBit);
  task_CME->SetEtaRangeMin(fEtaMin);
  task_CME->SetEtaRangeMax(fEtaMax);
  task_CME->SetPtRangeMin(fPtMin);
  task_CME->SetPtRangeMax(fPtMax);
  /// HardCoded at the moment. Can also be passed as AddTask Argument.
  task_CME->SetNSigmaCutTPC(2.0);               // For PID only; Does not apply to Inclusive Charged Tracks
  task_CME->SetTrackCutChi2Min(0.1);     
  task_CME->SetTrackCutdEdxMin(10.0);           
  task_CME->SetFlagUseKinkTracks(kFALSE);
  




  ///---> Now Pass data and containers to Analysis Object ----
  
  mgr->AddTask(task_CME);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CME, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_CME, 1, cOutPut1);
  
  printf("\n ================> AliAnalysisTaskCMXQA() Configured properly <==================\n");

  
}//Task Ends
