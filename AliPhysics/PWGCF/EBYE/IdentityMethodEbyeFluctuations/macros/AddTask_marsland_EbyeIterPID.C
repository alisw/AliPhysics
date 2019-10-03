AliAnalysisTask *AddTask_marsland_EbyeIterPID(Bool_t getFromAlien=kFALSE, TString configFileName = "Config_marsland_EbyeIterPID.C",Int_t settingType = 0,Int_t lhcPeriod = 1, Int_t lookUpTableIndex = 1, const char* suffix = "", Int_t containerNameMode=0)
{
  //
  //get the current analysis manager
  //
  /*
  =================================================================================================
  ======= eta Range --> [-0.8,0.8], momentum Range --> [0.2,3.2], dEdx Range --> [20,1020]* =======
  Parameter 1: getFromAlien     --> decide if the config file will be taken from alien or hera
  Parameter 2: configFileName   --> config file which should exist both in alicen and hera
  Parameter 3: settingType      --> an integer to decide which setting to use
  Parameter 4: lhcPeriod        --> an integer to decide which period 1,2 or 3
  Parameter 5: lookUpTableIndex --> an integer to decide which lookup table to use
  Parameter 6: suffix           --> used for naming conflicts in lego train --> each wagon has a different suffix
  Parameter 7: containerNameMode--> decide either dump output files in a TDirectoryFile or not. 0 without 1 and 2 with TDirectoryFile,
  =================================================================================================
  */
  std::cout << " Info::marsland: ===== In the AddTask_marsland_EbyeIterPID ===== " << std::endl;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_marsland_EbyeIterPID", "No analysis manager found.");
    return 0;
  }
  mgr->SetDebugLevel(0);
  //
  //==================================================
  //      Add task to the ANALYSIS manager
  //==================================================
  //
  //Get the current train configuration
  TString configBasePath;
  if (getFromAlien)
  {
    gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/m/marsland/PWGCF/EBYE/IdentityMethodEbyeFluctuations/macros/%s .",configFileName.Data()));
    configBasePath=Form("%s/",gSystem->pwd());
  } else {
    std::cout << " Info::marsland: Settings for local testing " << std::endl;
    configBasePath = "";
    // gSystem->Load("libANALYSIS");
    // gSystem->Load("libANALYSISalice");
    // gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    // gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    //
    // MC closure for higher moments
    settingType = 0;      // 4 for Real data 100 for full MC
    lhcPeriod   = 2;
    lookUpTableIndex =0;
    suffix = "test";
    containerNameMode=0;
  }
  TString configFilePath(configBasePath+configFileName);
  std::cout << " Info::marsland: Configpath:  " << configFilePath << std::endl;
  //
  gROOT->LoadMacro(configFilePath.Data());
  //
  TString combinedName;
  combinedName.Form("marslandEbyeIter_%s", suffix);
  AliAnalysisTaskEbyeIterPID* task = Config_marsland_EbyeIterPID(getFromAlien,settingType,lhcPeriod,lookUpTableIndex,combinedName);
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  task->SetIsMCtrue(hasMC);
  printf(" ========================= MC info %d ========================= \n",hasMC);
  mgr->AddTask(task);
  // mgr->Dump(); task->Dump();
  //
  //==================================================
  //      Add task to the PID Response
  //==================================================
  //
  // gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  // Bool_t isMC=kTRUE; // kTRUE in case of MC
  // AddTaskPIDResponse(isMC);
  //
  //================================================
  //              data containers
  //================================================
  //
  // ****** Do not forget to "DefineOutput(5, TTree::Class());" In the contructor of the task ******
  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *cinput, *coutput1, *coutput2, *coutput3, *coutput4;
  AliAnalysisDataContainer *coutput5, *coutput6, *coutput7, *coutput8, *coutput9;
  AliAnalysisDataContainer *coutput10, *coutput11, *coutput12, *coutput13,*coutput14;
  //
  //  find and connect input container // Output files --> File opening order is important
  cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);
  //
  // Get the output file name which is by default "AnalysisResults.root" --> Either write output into a TDirectoryFile or directly in "AnalysisResults.root"
  char* outputFileName = AliAnalysisManager::GetCommonFileName();
  TString dirName = "";
  TString fileDirStructure = "";
  TString listName = "";
  if (containerNameMode==0){
    fileDirStructure = Form("%s", outputFileName);  // TDirectoryFile name to put all containers; AnalysisResults.root --> trees and hists
    listName = "cleanHists";
  } else if (containerNameMode==1){
    dirName = "PWGCF_marsland";                     // TDirectoryFile name to put all containers; AnalysisResults.root --> PWGCF_marsland --> trees and hists
    fileDirStructure = Form("%s:%s", outputFileName, dirName.Data());
    listName = combinedName+"_cleanHists";
  } else if (containerNameMode==2){
    dirName = Form("PWGCF_marsland_%s",suffix);     // TDirectoryFile name to put all containers; AnalysisResults.root --> PWGCF_marsland_<setting> --> trees and hists
    fileDirStructure = Form("%s:%s", outputFileName, dirName.Data());
    listName = combinedName+"_cleanHists";
  }
  //
  // Output containers
  coutput1  = mgr->CreateContainer(listName,                      TList::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput2  = mgr->CreateContainer(combinedName+"_armPodTree",    TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput3  = mgr->CreateContainer(combinedName+"_mcFull",        TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput4  = mgr->CreateContainer(combinedName+"_mcGen",         TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput5  = mgr->CreateContainer(combinedName+"_fTreeMC",       TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput6  = mgr->CreateContainer(combinedName+"_fTreedEdxCheck",TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput7  = mgr->CreateContainer(combinedName+"_tracks",        TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput8  = mgr->CreateContainer(combinedName+"_dnchdeta",      TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput9  = mgr->CreateContainer(combinedName+"_fullacc",       TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput10 = mgr->CreateContainer(combinedName+"_resonance",     TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput11 = mgr->CreateContainer(combinedName+"_mcGenMoms",     TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput12 = mgr->CreateContainer(combinedName+"_events",        TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput13 = mgr->CreateContainer(combinedName+"_dscaled",       TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  coutput14 = mgr->CreateContainer(combinedName+"_mcMoms",        TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
  mgr->ConnectOutput (task,  1, coutput1);
  mgr->ConnectOutput (task,  2, coutput2);
  mgr->ConnectOutput (task,  3, coutput3);
  mgr->ConnectOutput (task,  4, coutput4);
  mgr->ConnectOutput (task,  5, coutput5);
  mgr->ConnectOutput (task,  6, coutput6);
  mgr->ConnectOutput (task,  7, coutput7);
  mgr->ConnectOutput (task,  8, coutput8);
  mgr->ConnectOutput (task,  9, coutput9);
  mgr->ConnectOutput (task,  10, coutput10);
  mgr->ConnectOutput (task,  11, coutput11);
  mgr->ConnectOutput (task,  12, coutput12);
  mgr->ConnectOutput (task,  13, coutput13);
  mgr->ConnectOutput (task,  14, coutput14);

  std::cout << " Info::marsland: === Containers are ready === " << std::endl;
  return task;
}
