
void AddTaskCVEUNI(Int_t gFilterBit = 768, Double_t fPtMin=0.2, Double_t fPtMax=10.0, Double_t fEtaMin=-0.8, Double_t fEtaMax=0.8,
		Int_t gNclustTPC=70, TString sCentEstimator="V0M", Double_t fCentralityMin=0., Double_t fCentralityMax=90.,
		Float_t fVzMin = -10.0, Float_t fVzMax = 10.0, TString sTrigger="kINT7", Int_t fparticle=3,
		Double_t nSigTPC = 3.0, Double_t nSigTOF = 3.0, Int_t vnHarmonic=2,Double_t fEtaGapNeg=-0.1,Double_t fEtaGapPos=0.1,
		TString sMCfilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/HijingMC_LHC18q_FB768_DeftCut.root",
		TString sNUAFilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/wgtCharge_NUAFB768NoPUcutRun296244.root",
		const char *suffix = "")
{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCVE                             \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TString TaskCVEPID;
  TaskCVEPID.Form("gTaskCVECent%d_%d_%s", gCentMin, gCentMax, suffix);

  AliAnalysisTaskCVEUNI *task_CVE = new AliAnalysisTaskCVEUNI(TaskCVEPID);

  ///-------> Analysis Object Created, now pass the arguments
  
  task_CVE->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    task_CVE->SelectCollisionCandidates(AliVEvent::kMB);
  }
  ///Event cuts:
  task_CVE->SetCentralityPercentileMin(fCentralityMin);
  task_CVE->SetCentralityPercentileMax(fCentralityMax);
  task_CVE->SetVzRangeMin(fVzMin);
  task_CVE->SetVzRangeMax(fVzMax);
  task_CVE->SetParticle(fparticle);
  
  if(sCentEstimator=="V0" || sCentEstimator=="V0M"){ 
    task_CVE->SetCentralityEstimator("V0M");    
  }
  else{
    task_CVE->SetCentralityEstimator(sCentEstimator);          //  use the Estimator provided in AddTask.
  }
     

  //Track cuts:
  task_CVE->SetTrackCutNclusterMin(gNclustTPC);
  task_CVE->SetFilterBit(gFilterBit);
  task_CVE->SetEtaRangeMin(fEtaMin);
  task_CVE->SetEtaRangeMax(fEtaMax);
  task_CVE->SetPtRangeMin(fPtMin);
  task_CVE->SetPtRangeMax(fPtMax);

  task_CVE->SetEtaNeg(fEtaGapNeg);
  task_CVE->SetEtaPos(fEtaGapPos);
  
  /// HardCoded at the moment. Can also be passed as AddTask Argument.
  task_CVE->SetNSigmaCutTPC(nSigTPC);               // For PID only; Does not apply to Inclusive Charged Tracks
  task_CVE->SetNSigmaCutTPC(nSigTOF);       
  task_CVE->SetTrackCutChi2Min(0.1);     
  task_CVE->SetTrackCutdEdxMin(10.0);           
  task_CVE->SetFlagUseKinkTracks(kFALSE);
  task_CVE->SetCumulantHarmonic(vnHarmonic);
  





  //--------------------------------------------------------------------------
  TFile *fMCFile = TFile::Open(sMCfilePath,"READ");
  TList *fListMC=NULL;
  
  if(fMCFile) {

    fListMC = dynamic_cast <TList*> (fMCFile->FindObjectAny("fMcEffiHij"));

    if(fListMC) {
      task_CVE->SetListForTrkCorr(fListMC); 
    }
    else{
      printf("\n\n *** ERROR: MC file Exist, But TList Not Found!!! \n please check name \n\n");
    }
  }
  else{
    printf("\n\n *** AddTastSimple::ERROR => MC file path/name Wrong!!! \n please check name \n\n");
  }
  


 //--------------------------------------------------------------------------
  TFile* fNUAFile = TFile::Open(sNUAFilePath,"READ");
  TList* fListNUA=NULL;

  if(fNUAFile) {
    
    fListNUA = dynamic_cast <TList*> (fNUAFile->FindObjectAny("fNUA_ChPosChNeg"));
    //std::cout<<" \n ==============> List found for NUA, here is all the histograms : ";
    //fListNUA->ls();

    if(fListNUA) {
      task_CVE->SetListForNUACorr(fListNUA);
    }
    else{
      printf("\n\n *** AddTastSimple::ERROR => NUA file Exist, But fList name is wrong!!\n please check name \n\n");
    }
  }
  






  

  ///---> Now Pass data and containers to Analysis Object ----
  
  mgr->AddTask(task_CVE);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CVE, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_CVE, 1, cOutPut1);
  
 
  printf("\n ================> AddTaskCVE() Configured properly <==================\n",);

  //return task_CVE;

}//Task Ends
