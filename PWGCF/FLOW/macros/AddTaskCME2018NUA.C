
void AddTaskCME2018NUA(Int_t whichData = 2018, TString period = "q", Int_t gFilterBit = 768, Double_t fPtMin=0.2, Double_t fPtMax=10.0,Double_t fdcaxy=2.4, Double_t fdcaz=3.2, Double_t fChi2=4.0,Double_t fSlope=3.45, Float_t fConst=100, Bool_t bSkipPileUp=kFALSE, Double_t fEtaMin=-0.8, Double_t fEtaMax=0.8, Int_t gNclustTPC=70, TString sCentEstimator="V0M", Double_t fCentralityMin=0., Double_t fCentralityMax=90.,Float_t fVzMin = -10.0, Float_t fVzMax = 10.0, TString sTrigger="kINT7", Int_t fparticle=3,Double_t nSigTPC = 3.0, Double_t nSigTOF = 3.0, Int_t vnHarmonic=2,Double_t fEtaGapNeg=-0.1,Double_t fEtaGapPos=0.1,TString sMCfilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/HijingMC_LHC18q_FB768_DeftCut.root",TString sNUAFilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/wgtCharge_NUAFB768NoPUcutRun296244.root",TString sEVNTWGTFilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/wgtCharge_NUAFB768NoPUcutRun296244.root",const char *suffix = "")
{
  TGrid::Connect("alien://");
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AddTaskCME2018NUA                             \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TString TaskCMEPID;
  TaskCMEPID.Form("gTaskCMWCent%d_%d_%s", gCentMin, gCentMax, suffix);

  AliAnalysisTaskCME2018NUA *task_CMENUA = new AliAnalysisTaskCME2018NUA(TaskCMEPID);

  ///-------> Analysis Object Created, now pass the arguments
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    task_CMENUA->SelectCollisionCandidates(AliVEvent::kMB);
    printf("\n =========> AddTaskCMW::Info() Trigger = kMB  \n");
  }
  else if(sTrigger=="kSemiCentral" || sTrigger=="SemiCentral" || sTrigger=="semicentral"){
    task_CMENUA->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    printf("\n =========> AddTaskCMW::Info() Trigger = kSemiCentral \n");
  }
  else if(sTrigger=="kCentral" || sTrigger=="Central" || sTrigger=="central"){
    task_CMENUA->SelectCollisionCandidates(AliVEvent::kCentral);
    printf("\n =========> AddTaskCMW::Info() Trigger = kCentral \n");
  }
  else{//if trigger==kINT7 or no trigger provided:
    task_CMENUA->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskCMW::Info() Trigger = kINT7 \n");
  }
  
  
  ///Event cuts:
  task_CMENUA->SetCentralityPercentileMin(fCentralityMin);
  task_CMENUA->SetCentralityPercentileMax(fCentralityMax);
  task_CMENUA->SetVzRangeMin(fVzMin);
  task_CMENUA->SetVzRangeMax(fVzMax);
  task_CMENUA->SetParticle(fparticle);
  task_CMENUA->SetPileUpCutParam(fSlope,fConst);
  task_CMENUA->SetFlagSkipPileUpCuts(bSkipPileUp);  
  if (whichData == 2018) {
	  if (period == "q") {
		  task_CMENUA->SetupPileUpRemovalFunctions18qPass3();
	  } else if (period == "r") {
		  task_CMENUA->SetupPileUpRemovalFunctions18rPass3();
	  } else {
		  task_CMENUA->SetupPileUpRemovalFunctions18();
	  }
  }
  //

  if(sCentEstimator=="V0" || sCentEstimator=="V0M"){ 
    task_CMENUA->SetCentralityEstimator("V0M");    
  }
  else{
    task_CMENUA->SetCentralityEstimator(sCentEstimator);          //  use the Estimator provided in AddTask.
  }
     

  //Track cuts:
  task_CMENUA->SetTrackCutNclusterMin(gNclustTPC);
  task_CMENUA->SetFilterBit(gFilterBit);
  task_CMENUA->SetEtaRangeMin(fEtaMin);
  task_CMENUA->SetEtaRangeMax(fEtaMax);
  task_CMENUA->SetPtRangeMin(fPtMin);
  task_CMENUA->SetPtRangeMax(fPtMax);
  task_CMENUA->SetDCAXYRangeMax(fdcaxy);
  task_CMENUA->SetDCAZRangeMax(fdcaz);
  task_CMENUA->SetChi2Range(fChi2);
  task_CMENUA->SetEtaNeg(fEtaGapNeg);
  task_CMENUA->SetEtaPos(fEtaGapPos);
  /// HardCoded at the moment. Can also be passed as AddTask Argument.
  task_CMENUA->SetNSigmaCutTPC(nSigTPC);               // For PID only; Does not apply to Inclusive Charged Tracks
  task_CMENUA->SetNSigmaCutTPC(nSigTOF);       
  task_CMENUA->SetTrackCutChi2Min(0.1);     
  task_CMENUA->SetTrackCutdEdxMin(10.0);           
  task_CMENUA->SetFlagUseKinkTracks(kFALSE);
  task_CMENUA->SetCumulantHarmonic(vnHarmonic);
  





  //--------------------------------------------------------------------------
  TFile *fMCFile = TFile::Open(sMCfilePath,"READ");
  TList *fListMC=NULL;
  
  if(fMCFile) {

    fListMC = dynamic_cast <TList*> (fMCFile->FindObjectAny("fMcEffiHij"));

    if(fListMC) {
      task_CMENUA->SetListForTrkCorr(fListMC); 
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
    std::cout<<" \n ==============> List found for NUA, here is all the histograms : ";
    //fListNUA->ls();

    if(fListNUA) {
      task_CMENUA->SetListForNUACorr(fListNUA);
    }
    else{
      printf("\n\n *** AddTask::WARNING => NUA file Exist,But TList Not Found!!\n AddTask::Info() ===> NO NUA Correction!! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING => NUA file not Found!!\n AddTask::Info() ===> NO NUA Correction!! \n\n");
  }



  //-----------------------------------------------------------------------------


  TFile* fEVNTWGTFile = TFile::Open(sEVNTWGTFilePath,"READ");
  TList* fListEVNTWGT=NULL;

  if(fEVNTWGTFile) {    
    fListEVNTWGT = dynamic_cast <TList*> (fEVNTWGTFile->FindObjectAny("fNUA_ChPosChNeg"));
    std::cout<<" \n ==============> List found for EventWeight, here is all the histograms : ";
    //fListEVNTWGT->ls();

    if(fListEVNTWGT) {
      task_CMENUA->SetListForV0MCorr(fListEVNTWGT);
    }
    else{
      printf("\n\n *** AddTask::WARNING => EVNTWGT file Exist,But TList Not Found!!\n AddTask::Info() ===> NO EVNTWGT Correction!! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING => EVNTWGT file not Found!!\n AddTask::Info() ===> NO EVNTWGT Correction!! \n\n");
  }
  






  

  ///---> Now Pass data and containers to Analysis Object ----
  
  mgr->AddTask(task_CMENUA);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CMENUA, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_CMENUA, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTaskCMW() Configured properly <==================\n\n");

  //return task_CMENUA;

}//Task Ends
