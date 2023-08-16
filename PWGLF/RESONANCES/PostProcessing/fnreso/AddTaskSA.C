
void AddTaskSA(const char *suffix = "CMWchrg", Int_t frame =1, TString sMCfilePath = "a.root",TString sNUAFilePath = "b.root", TString sEVNTWGTFilePath = "c.root")

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

  AliAnalysisTaskSA *task_CMW = new AliAnalysisTaskSA(TaskCMWPID);
  task_CMW->Setframe(frame);




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
    std::cout<<" \n ==============> List found for NUA, here is all the histograms : ";
  
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



  //-----------------------------------------------------------------------------


  TFile* fEVNTWGTFile = TFile::Open(sEVNTWGTFilePath,"READ");
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
  }
  





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
