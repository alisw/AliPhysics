
void AddTaskCMW(Int_t gFilterBit = 768, Double_t fPtMin=0.2, Double_t fPtMax=10.0, Double_t fEtaMin=-0.8, Double_t fEtaMax=0.8,
		Int_t gNclustTPC=70, TString sCentEstimator="V0M", Double_t fCentralityMin=0., Double_t fCentralityMax=90.,
		Float_t fVzMin = -10.0, Float_t fVzMax = 10.0, TString sTrigger="kINT7", Int_t fparticle=3,
		Double_t nSigTPC = 3.0, Double_t nSigTOF = 3.0, Int_t vnHarmonic=2,
		TString sMCfilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/HijingMC_LHC18q_FB768_DeftCut.root",
		TString sNUAFilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/wgtCharge_NUAFB768NoPUcutRun296244.root",
		const char *suffix = "")
{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCMW                             \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TString TaskCMWPID;
  TaskCMWPID.Form("gTaskCMWCent%d_%d_%s", gCentMin, gCentMax, suffix);

  AliAnalysisTaskCMW *task_CMW = new AliAnalysisTaskCMW(TaskCMWPID);

  ///-------> Analysis Object Created, now pass the arguments
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    task_CMW->SelectCollisionCandidates(AliVEvent::kMB);
    printf("\n =========> AddTaskCMW::Info() Trigger = kMB  \n");
  }
  else if(sTrigger=="kSemiCentral" || sTrigger=="SemiCentral" || sTrigger=="semicentral"){
    task_CMW->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    printf("\n =========> AddTaskCMW::Info() Trigger = kSemiCentral \n");
  }
  else if(sTrigger=="kCentral" || sTrigger=="Central" || sTrigger=="central"){
    task_CMW->SelectCollisionCandidates(AliVEvent::kCentral);
    printf("\n =========> AddTaskCMW::Info() Trigger = kCentral \n");
  }
  else{//if trigger==kINT7 or no trigger provided:
    task_CMW->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskCMW::Info() Trigger = kINT7 \n");
  }
  
  
  ///Event cuts:
  task_CMW->SetCentralityPercentileMin(fCentralityMin);
  task_CMW->SetCentralityPercentileMax(fCentralityMax);
  task_CMW->SetVzRangeMin(fVzMin);
  task_CMW->SetVzRangeMax(fVzMax);
  task_CMW->SetParticle(fparticle);
  
  if(sCentEstimator=="V0" || sCentEstimator=="V0M"){ 
    task_CMW->SetCentralityEstimator("V0M");    
  }
  else{
    task_CMW->SetCentralityEstimator(sCentEstimator);          //  use the Estimator provided in AddTask.
  }
     

  //Track cuts:
  task_CMW->SetTrackCutNclusterMin(gNclustTPC);
  task_CMW->SetFilterBit(gFilterBit);
  task_CMW->SetEtaRangeMin(fEtaMin);
  task_CMW->SetEtaRangeMax(fEtaMax);
  task_CMW->SetPtRangeMin(fPtMin);
  task_CMW->SetPtRangeMax(fPtMax);
  /// HardCoded at the moment. Can also be passed as AddTask Argument.
  task_CMW->SetNSigmaCutTPC(nSigTPC);               // For PID only; Does not apply to Inclusive Charged Tracks
  task_CMW->SetNSigmaCutTPC(nSigTOF);       
  task_CMW->SetTrackCutChi2Min(0.1);     
  task_CMW->SetTrackCutdEdxMin(10.0);           
  task_CMW->SetFlagUseKinkTracks(kFALSE);
  task_CMW->SetCumulantHarmonic(vnHarmonic);
  





  //--------------------------------------------------------------------------
  TFile *fMCFile = TFile::Open(sMCfilePath,"READ");
  TList *fListMC=NULL;
  
  if(fMCFile) {

    fListMC = dynamic_cast <TList*> (fMCFile->FindObjectAny("fMcEffiHij"));

    if(fListMC) {
      task_CMW->SetListForTrkCorr(fListMC); 
    }
    else{
      printf("\n\n *** AddTask::WARNING \n => MC file Exist, But TList Not Found!!! \n AddTask::Info() ===> NO MC Correction!! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING \n => no MC file!!! \n AddTask::Info() ===> NO MC Correction!! \n\n");
  }
  


 //--------------------------------------------------------------------------
  TFile* fNUAFile = TFile::Open(sNUAFilePath,"READ");
  TList* fListNUA=NULL;

  if(fNUAFile) {    
    fListNUA = dynamic_cast <TList*> (fNUAFile->FindObjectAny("fNUA_ChPosChNeg"));
    //std::cout<<" \n ==============> List found for NUA, here is all the histograms : ";
    //fListNUA->ls();

    if(fListNUA) {
      task_CMW->SetListForNUACorr(fListNUA);
    }
    else{
      printf("\n\n *** AddTask::WARNING => NUA file Exist,But TList Not Found!!\n AddTask::Info() ===> NO NUA Correction!! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING => NUA file not Found!!\n AddTask::Info() ===> NO NUA Correction!! \n\n");
  }
  






  

  ///---> Now Pass data and containers to Analysis Object ----
  
  mgr->AddTask(task_CMW);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CMW, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_CMW, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTaskCMW() Configured properly <==================\n\n",);

  //return task_CMW;

}//Task Ends
