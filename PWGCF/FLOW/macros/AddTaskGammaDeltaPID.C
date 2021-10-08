
void AddTaskGammaDeltaPID(Int_t gFilterBit = 768,Double_t fPtMin=0.2,Double_t fPtMax=2.0,Double_t fEtaMin=-0.8, Double_t fEtaMax=0.8,Double_t fChi2max=4.0,Int_t gNclustTPC=70, Int_t fparticle=0,Double_t nSigTPC = 3.0, Double_t nSigTOF = 3.0, Bool_t bSkipPileUp=kFALSE, TString sCentEstimator="V0M", Float_t fVzMin = -10.0, Float_t fVzMax = 10.0, TString sTrigger="kINT7", Int_t vnHarmonic=2, TString sDetForEP="TPC", TString sMCfilePath="alien:///alice/cern.ch/user/m/mhaque/nuanue18/HijingMC_LHC18q_FB768_DeftCut.root", TString sNUAFilePath = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/wgtCharge_NUAFB768NoPUcutRun296244.root", TString sDetWgtsFile = "alien:///alice/cern.ch/user/m/mhaque/nuanue18/wgtCharge_NUAFB768NoPUcutRun296244.root", Bool_t bSkipAnalysis=kFALSE, const char *suffix = "")
{

  printf("===================================================================================\n");
  printf("                   Initialising Task: AddTaskGammaDeltaPID                         \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName; // common outfile filename
  list1OutName += ":Results";         // This directory contains result histograms


  TString TaskCMWPID;
  TaskCMWPID.Form("gTaskBugTest%d_%d_%s", gFilterBit, gNclustTPC, suffix);

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
  else if(sTrigger=="kAny" || sTrigger=="kAll"){
    task_CMW->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kSemiCentral | AliVEvent::kCentral);
  }
  else{//if trigger==kINT7 or no trigger provided:
    task_CMW->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskCMW::Info() Trigger = kINT7 \n");
  }
  

  ///Set Event cuts:
  task_CMW->SetVzRangeMin(fVzMin);
  task_CMW->SetVzRangeMax(fVzMax);
  task_CMW->SetFlagSkipPileUpCuts(bSkipPileUp);  
  task_CMW->SetFlagSkipAnalysis(bSkipAnalysis);



  cout<<"=========> AddTaskCMW::Info() setting Event Plane Det: "<<sDetForEP<<endl;
  task_CMW->SetDetectorforEventPlane(sDetForEP);


  
  if(sCentEstimator=="V0" || sCentEstimator=="V0M"){ 
    task_CMW->SetCentralityEstimator("V0M");    
  }
  else{
    task_CMW->SetCentralityEstimator(sCentEstimator);  //  use the Estimator provided in AddTask.
  }



  
  //Set Track cuts:

  task_CMW->SetPtRangeMin(fPtMin);
  task_CMW->SetPtRangeMax(fPtMax);
  task_CMW->SetEtaRangeMin(fEtaMin);
  task_CMW->SetEtaRangeMax(fEtaMax);
  task_CMW->SetTrackCutChi2Min(0.1);
  task_CMW->SetTrackCutdEdxMin(10.0);  
  task_CMW->SetFilterBit(gFilterBit);
  task_CMW->SetNSigmaCutTPC(nSigTPC);    /// For PID only.Does not apply to Inclusive Charged Tracks
  task_CMW->SetNSigmaCutTOF(nSigTOF);
  task_CMW->SetParticlePID(fparticle);
  task_CMW->SetTrackCutChi2Max(fChi2max);
  task_CMW->SetFlagUseKinkTracks(kFALSE);
  task_CMW->SetCumulantHarmonic(vnHarmonic);
  task_CMW->SetTrackCutNclusterMin(gNclustTPC);  
 
  Bool_t bFillLambda=kFALSE;
  task_CMW->SetFlagAnalyseLambda(bFillLambda);
   ///  -----> Separate AddTask Added For Lambda-X correlation
   ///  AddTaskGammaDeltaPID.C  
 




  //========================= Setup Correction Files ======================> 
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
    //fListNUA->ls();
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
 
  TFile* fV0ZDCWgtsFile = TFile::Open(sDetWgtsFile,"READ");
  TList* fListDetWgts=NULL;

  if(fV0ZDCWgtsFile) {    
    fListDetWgts = dynamic_cast <TList*> (fV0ZDCWgtsFile->FindObjectAny("fWgtsV0ZDC"));
    std::cout<<" \n ==============> TList found for V0/ZDC wgts.. GOOD! ";
    // fListDetWgts->ls(); 

    if(fListDetWgts) {
      task_CMW->SetListForV0MCorr(fListDetWgts);
    }
    else{
      printf("\n\n *** AddTask::WARNING => V0/ZDC Weights file Exist, But TList Not Found!!");
      printf("\n May be wrong TList name? No Correction for V0/ZDC !! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING => NO File Found for V0/ZDC Wgts!!\n AddTask::Info() ===> No V0/ZDC Correction!! \n\n");
  }
  //=================================================================================



  

  ///---> Now Pass data and containers to Analysis Object ----
 
  mgr->AddTask(task_CMW);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CMW, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(task_CMW, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTask was Configured properly... <==================\n\n");

  //return task_CMW;

}//Task Ends
