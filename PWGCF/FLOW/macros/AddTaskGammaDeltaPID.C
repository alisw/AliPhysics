
void AddTaskGammaDeltaPID(Int_t gFilterBit = 768,Double_t fPtMin=0.2,Double_t fPtMax=2.0,Double_t fV0DautPtMax=2.0,Double_t fEtaMin=-0.8, Double_t fEtaMax=0.8,Double_t fChi2=4.0, Int_t gNclustTPC=70, Int_t fparticle=3,Double_t nSigTPC = 3.0, Double_t nSigTOF = 3.0, Bool_t bSkipPileUp=kFALSE, TString sCentEstimator="V0M",Float_t fVzMin = -10.0, Float_t fVzMax = 10.0,TString sTrigger="kINT7",Int_t vnHarmonic=2,TString sMCfilePath="alien:///alice/cern.ch/user/m/mhaque/nuanue18/HijingMC_LHC18q_FB768_DeftCut.root",TString sNUAFilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/wgtCharge_NUAFB768NoPUcutRun296244.root",TString sEvtWgtPath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/wgtCharge_NUAFB768NoPUcutRun296244.root",Bool_t bSkipAnalysis=kFALSE,Bool_t bFillLambda=kTRUE,Double_t fMassMean = 1.115683, Double_t fV0MassCut = 0.005, Double_t fV0CosPAmin = 0.995, Double_t fV0RapidityMax=0.5, Double_t fV0DecLenMin=3.0, Double_t fV0DecLenMax=100, Double_t fV0DCAToPrimVtx=1.5, Double_t fV0DcaDiffDautMax=1.0, Double_t fV0DautDCAToPrimVtxMin = 0.02, const char *suffix = "")
{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                Initialising Task: AliAnalysisTaskGammaDeltaPID                  \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":Results";         // This directory contains result histograms

  //Int_t gCentMin = fCentralityMin;
  //Int_t gCentMax = fCentralityMax;
  Int_t gCentMin = 0;
  Int_t gCentMax = 90.;  

  TString TaskCMWPID;
  TaskCMWPID.Form("gTaskCMWCent%d_%d_%s", gCentMin, gCentMax, suffix);

  AliAnalysisTaskGammaDeltaPID *task_CMW = new AliAnalysisTaskGammaDeltaPID(TaskCMWPID);

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
  else if(sTrigger=="kAny" ||sTrigger=="kAll"){
    task_CMW->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kSemiCentral | AliVEvent::kCentral);
  }
  else{//if trigger==kINT7 or no trigger provided:
    task_CMW->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskCMW::Info() Trigger = kINT7 \n");
  }
  

  Float_t fSlope=3.45;
  Float_t fConst=15000;
  ///Event cuts:
  task_CMW->SetCentralityPercentileMin(gCentMin);
  task_CMW->SetCentralityPercentileMax(gCentMax);
  task_CMW->SetVzRangeMin(fVzMin);
  task_CMW->SetVzRangeMax(fVzMax);
  task_CMW->SetParticle(fparticle);
  task_CMW->SetPileUpCutParam(fSlope,fConst);
  task_CMW->SetFlagSkipPileUpCuts(bSkipPileUp);  
  task_CMW->SetFlagSkipAnalysis(bSkipAnalysis);  
  //

  if(sCentEstimator=="V0" || sCentEstimator=="V0M"){ 
    task_CMW->SetCentralityEstimator("V0M");    
  }
  else{
    task_CMW->SetCentralityEstimator(sCentEstimator);          //  use the Estimator provided in AddTask.
  }
     

  //// Hardcode some cuts: Can be made variable if Needed.

  // Gap in btween the TPC event plane:
  Double_t fEtaGapNeg = -0.1;   //// |-0.8 ~ fEtaGapNeg| <-gap-> | fEtaGapPos ~ 0.8|
  Double_t fEtaGapPos = +0.1;
  Double_t fDCAcutxy = 2.4;     //// Redundant As this is superseeded by the FB. 
  Double_t fDCAcutz  = 3.2;     //// Redundant As this is superseeded by the FB.   

  
  //Track cuts:
  task_CMW->SetTrackCutNclusterMin(gNclustTPC);
  task_CMW->SetFilterBit(gFilterBit);
  task_CMW->SetEtaRangeMin(fEtaMin);
  task_CMW->SetEtaRangeMax(fEtaMax);
  task_CMW->SetPtRangeMin(fPtMin);
  task_CMW->SetPtRangeMax(fPtMax);
  task_CMW->SetDCAXYRangeMax(fDCAcutxy);
  task_CMW->SetDCAZRangeMax(fDCAcutz);
  task_CMW->SetChi2Range(fChi2);
  task_CMW->SetEtaNeg(fEtaGapNeg);
  task_CMW->SetEtaPos(fEtaGapPos);
  task_CMW->SetNSigmaCutTPC(nSigTPC);              // For PID only; Does not apply to Inclusive Charged Tracks
  task_CMW->SetNSigmaCutTPC(nSigTOF);              // For PID only; Does not apply to Inclusive Charged Tracks
  task_CMW->SetTrackCutChi2Min(0.1);               // Min. Chi2 per cluster.  Usually not varied!
  task_CMW->SetTrackCutdEdxMin(10.0);              // Min. dEdx signal per track. Usually not varied!
  task_CMW->SetFlagUseKinkTracks(kFALSE);
  task_CMW->SetCumulantHarmonic(vnHarmonic);
  
  ///--------- Lambda Related cuts:-----------
  
  /// Re-use Some AddTask Variables:
  Int_t    fV0DautTPCnclsMax = gNclustTPC;  
  Double_t fV0PtMin = fPtMin;  
  Double_t fV0DautEtaMax = fEtaMax;
  Double_t fV0DautNsigmaMax = nSigTPC;


  
  task_CMW->SetFlagAnalyseLambda(bFillLambda);
  task_CMW->SetV0PtMin(fV0PtMin);                        
  task_CMW->SetV0CPAMin(fV0CosPAmin);                             
  task_CMW->SetMassMean(fMassMean);                             
  task_CMW->SetLambdaMassCut(fV0MassCut);                   
  task_CMW->SetV0RapidityMax(fV0RapidityMax);                   
  task_CMW->SetV0DecayLengthMax(fV0DecLenMax);             
  task_CMW->SetV0DecayLengthMin(fV0DecLenMin);             
  task_CMW->SetV0DCAToPrimVtxMax(fV0DCAToPrimVtx);           
  task_CMW->SetV0DcaBetweenDaughtersMax(fV0DcaDiffDautMax);
  //V0 Daughter Cut
  task_CMW->SetDaughtersPtMax(fV0DautPtMax);                 
  task_CMW->SetDaughtersEtaMax(fV0DautEtaMax);                 
  task_CMW->SetDaughtersNsigma(fV0DautNsigmaMax);              
  task_CMW->SetDaughtersTPCNclsMin(fV0DautTPCnclsMax);       
  task_CMW->SetDaughtersDCAToPrimVtxMin(fV0DautDCAToPrimVtxMin);


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
  std::cout<<" NUA file Path "<<sNUAFilePath.Data()<<std::endl;
  TFile* fNUAFile = TFile::Open(sNUAFilePath,"READ");
  TList* fListNUA=NULL;

  //if(fNUAFile->IsOpen()) {
  if(fNUAFile){    
    fListNUA = dynamic_cast <TList*> (fNUAFile->FindObjectAny("fNUA_ChPosChNeg"));
    std::cout<<" \n ==============> TList found for NUA, here is all the histograms : "<<std::endl;
    fListNUA->ls();
    if(fListNUA) {
      task_CMW->SetListForNUACorr(fListNUA);
    }
    else{
      printf("\n\n *** AddTask::WARNING => NUA file Exist,But TList Not Found!!\n AddTask::Info() ===> NO NUA Correction!! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING => NUA file not Found or Wrong path Set in AddTask Macro!! \n\n");
  }



  //-----------------------------------------------------------------------------

  /*
  TFile* fEVNTWGTFile = TFile::Open(sEvtWgtPath,"READ");
  TList* fListEVNTWGT=NULL;

  if(fEVNTWGTFile) {    
    fListEVNTWGT = dynamic_cast <TList*> (fEVNTWGTFile->FindObjectAny("fNUA_ChPosChNeg"));
    std::cout<<" \n ==============> List found for EventWeight, here is all the histograms : ";
    //fListEVNTWGT->ls();

    if(fListEVNTWGT) {
      task_CMW->SetListForV0MCorr(fListEVNTWGT);
    }
    else{
      printf("\n\n *** AddTask::WARNING => EVNTWGT file Exist,But TList Not Found!!\n AddTask::Info() ===> NO EVNTWGT Correction!! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING => EVNTWGT file not Found!!\n AddTask::Info() ===> NO EVNTWGT Correction!! \n\n");
  }*/






  

  ///---> Now Pass data and containers to Analysis Object ----
  
  mgr->AddTask(task_CMW);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CMW, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_CMW, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTaskCMW() Configured properly <==================\n\n");

  //return task_CMW;

}//Task Ends
