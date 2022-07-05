// R__ADD_INCLUDE_PATH($PWD)
#include "AliAnalysisTaskEbyeIterPID.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom.h"

AliAnalysisTask *AddTask_marsland_EbyeIterPID(Bool_t getFromAlien=kFALSE, TString configFileName = "Config_marsland_EbyeIterPID.C",Int_t settingType = 0,Int_t year = 1, TString periodName="15o", Int_t passIndex = 0, Int_t lookUpTableIndex = 1, const char* suffix = "", Int_t containerNameMode=0)
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGUDbase");
  gSystem->Load("libTPCcalib");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGLFspectra");
  gRandom->SetSeed(0);
  //
  //get the current analysis manager
  //
  /*
  =================================================================================================
  ======= eta Range --> [-0.8,0.8], momentum Range --> [0.2,3.2], dEdx Range --> [20,1020]* =======
  Parameter 1: getFromAlien     --> decide if the config file will be taken from alien or hera
  Parameter 2: configFileName   --> config file which should exist both in alicen and hera
  Parameter 3: settingType      --> an integer to decide which setting to use
  Parameter 4: year             --> year
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
  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
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
    // gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/m/marsland/PWGCF/EBYE/IdentityMethodEbyeFluctuations/macros/%s .",configFileName.Data()));
    // TString configFileName = "Config_marsland_EbyeIterPID.C";
    gSystem->Exec(Form("alien_cp /alice/cern.ch/user/m/marsland/PWGCF/EBYE/IdentityMethodEbyeFluctuations/macros/%s file:%s",configFileName.Data(),configFileName.Data()));
    configBasePath=Form("%s/",gSystem->pwd());
  } else {
    std::cout << " Info::marsland: Settings for local testing " << std::endl;
    configBasePath = "";
    //
    // gSystem->Load("libANALYSIS");
    // gSystem->Load("libANALYSISalice");
    // gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    // gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    //
    // MC closure for higher moments
    // settingType = 64;   // 1 for Real data 50 for full MC
    // year   = 2; periodName="15o"; passIndex=2  // 1 for 10h, 2 for 15o, 3 for 18[q,r]
    lookUpTableIndex =0;
    suffix = "test";
    containerNameMode=0;
  }
  TString combinedName;
  combinedName.Form("marslandTIdentity_%s", suffix);
  TString configFilePath(configBasePath+configFileName);
  std::cout << " Info::marsland: Configpath:  " << configFilePath << " year = " << year << " --- period name = " << periodName << " --- pass = " << passIndex << " --- lookUpTableIndex = " << lookUpTableIndex << " --- settingType = " << settingType << std::endl;
  //

  AliAnalysisTaskEbyeIterPID* task(0x0);
  #ifdef __CLING__
      std::stringstream triggermakeradd;
      triggermakeradd << ".x " << configFilePath.Data() << "(";
      triggermakeradd << (getFromAlien ? "kTRUE" : "kFALSE") << ", ";
      triggermakeradd << settingType << ", ";
      triggermakeradd << year << ", ";
      triggermakeradd << "\"" << periodName.Data() << "\"" << ", ";
      triggermakeradd << passIndex << ", ";
      triggermakeradd << lookUpTableIndex << ", ";
      triggermakeradd << "\"" << combinedName.Data() << "\"";
      triggermakeradd << ")";
      std::string triggermakeraddstring = triggermakeradd.str();
      std::cout << triggermakeraddstring << std::endl;
      task = (AliAnalysisTaskEbyeIterPID*)gROOT->ProcessLine(triggermakeraddstring.c_str());
  #else
      gROOT->LoadMacro(configFilePath.Data());
      task = Config_marsland_EbyeIterPID(getFromAlien,settingType,year,periodName,passIndex,lookUpTableIndex,combinedName);
  #endif
  
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
  // gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  // AddTaskPIDResponse(hasMC, kTRUE, kFALSE, "", kFALSE, "TPC-OADB:COMMON/PID/data/TPCPIDResponseOADB_pileupCorr.root;TPC-Maps:$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCetaMaps_pileupCorr.root" );
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
  const char* outputFileName = AliAnalysisManager::GetCommonFileName();
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
  coutput12 = mgr->CreateContainer(combinedName+"_eventInfo",     TTree::Class(), AliAnalysisManager::kOutputContainer, fileDirStructure);
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
