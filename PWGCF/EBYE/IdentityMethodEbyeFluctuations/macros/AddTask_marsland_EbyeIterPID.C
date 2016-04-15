AliAnalysisTask *AddTask_marsland_EbyeIterPID(Bool_t getFromAlien=kFALSE,TString configFileName = "Config_marsland_EbyeIterPID.C",Int_t settingType = 0)
{
  //  
  //get the current analysis manager
  //
  /*
   *   =================================================================================================
   *   ======= eta Range --> [-0.8,0.8], momentum Range --> [0.2,3.2], dEdx Range --> [20,1020]* =======
   *   Parameter 1: getFromAlien   --> decide if the config file will be taken from alien or hera
   *   Parameter 2: configFileName --> config file which should exist both in alicen and hera
   *   Parameter 3: settingType    --> an integer to decide which setting to use 
   *   =================================================================================================
   * 
   *    Real data --> settingType = 
   *    0.)  THnSparse is used: StandardTPCITScuts              16EtaBin_150pBins_9centBins (REFERENCE settings)
   *    1.)  THnSparse is used: StandardTPCITScuts + TIGHTCUTS  16EtaBin_150pBins_9centBins (active length cut)
   *    2.)  THnSparse is used: StandardTPCITScuts + DCAxySmall 16EtaBin_150pBins_9centBins  
   *    3.)  THnSparse is used: StandardTPCITScuts + DCAxyLarge 16EtaBin_150pBins_9centBins  
   *    4.)  THnSparse is used: StandardTPCITScuts + cRows60    16EtaBin_150pBins_9centBins  
   *    5.)  THnSparse is used: StandardTPCITScuts + cRows100   16EtaBin_150pBins_9centBins  
   *    6.)  THnSparse is used: StandardTPCITScuts + centEstCL1 16EtaBin_150pBins_9centBins  
   *    7.)  THnSparse is used: StandardTPCITScuts + centEstTRK 16EtaBin_150pBins_9centBins  
   *    8.)  THnSparse is used: StandardTPCITScuts + Vz8        16EtaBin_150pBins_9centBins  
   *    9.)  THnSparse is used: StandardTPCITScuts + Vz12       16EtaBin_150pBins_9centBins  
   *    10.) THnSparse is used: StandardTPCITScuts + Chi2Small  16EtaBin_150pBins_9centBins  
   *    11.) THnSparse is used: StandardTPCITScuts + Chi2Large  16EtaBin_150pBins_9centBins  
   *    12.) THnSparse is used: (REFERENCE settings) + allCuts are filled  
   *    13.) THnSparse is used: (REFERENCE settings) + Bayesian Probabilities are filled  
   *    14.) THnSparse is used: (REFERENCE settings) + dEdxTree is filled
   *    15.) THnSparse is used: (REFERENCE settings) + 18centBins  
   *    16.) THnSparse is used: (REFERENCE settings) + centBin 10  
   *    17.) THnSparse is used: (REFERENCE settings) + centBin 5  
   *    18.) THnSparse is used: ITS OFF 
   *    19.) THnSparse is used: (REFERENCE settings) + THnSparse is used: number of eta bins = 32 
   * 
   *    MC data --> settingType = 
   *    20.) THnSparse is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings) MC CLOSURE
   *    21.) FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins (REFERENCE settings)
   *    22.) FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins ETA DEPENDENCE
   *    23.) FastSimul is used: StandardTPCITScuts 8EtaBin_150pBins_9centBins Momentum DEPENDENCE
   *    24.) FullSinul is used: StandardTPCITScuts 16EtaBin_150pBins_9centBins EffMatrix
   *    25.) FullSinul is used: Tight Cuts         16EtaBin_150pBins_9centBins EffMatrix
   */
  cout << " ===== In the AddTask_marsland_EbyeIterPID ===== " << endl;
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
    configBasePath = "/hera/alice/marsland/train/trunk/marsland_EbyeRatios/";
  }
  TString configFilePath(configBasePath+configFileName);
  std::cout << "Configpath:  " << configFilePath << std::endl;  
  gROOT->LoadMacro(configFilePath.Data());
  AliAnalysisTaskEbyeIterPID* task = Config_marsland_EbyeIterPID(settingType);
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  task->SetIsMCtrue(hasMC);
  printf(" ========================= MC info %d ========================= \n",hasMC);
  task->Dump();
  mgr->AddTask(task);
  //   
  //==================================================
  //      Add task to the PID Response 
  //==================================================
  //   
  //   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  //   Bool_t isMC=kTRUE; // kTRUE in case of MC
  //   AddTaskPIDResponse(isMC); 
  // 
  //================================================
  //              data containers
  //================================================
  //   
  // ****** Do not forget to "DefineOutput(5, TTree::Class());" In the contructor of the task ******
  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *cinput,   *coutput1, *coutput2, *coutput3, *coutput4;
  AliAnalysisDataContainer *coutput5, *coutput6, *coutput7, *coutput8, *coutput9, *coutput10, *coutput11, *coutput12, *coutput13;
  
  //  find input container // Output files --> File opening order is important
  TString results = "AnalysisResults.root";
  cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task,  0, cinput);
  // Output containers 
  coutput1  = mgr->CreateContainer("cleanHists",    TList::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput2  = mgr->CreateContainer("tidtree",       TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput3  = mgr->CreateContainer("MCtidtree",     TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput4  = mgr->CreateContainer("dataTree",      TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput5  = mgr->CreateContainer("armPodTree",    TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput6  = mgr->CreateContainer("mcRec",         TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput7  = mgr->CreateContainer("mcGen",         TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput8  = mgr->CreateContainer("fTreeMC",       TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput9  = mgr->CreateContainer("fTreedEdxCheck",TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput10 = mgr->CreateContainer("fTreeBayes",    TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput11 = mgr->CreateContainer("fTreeCuts",     TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput12 = mgr->CreateContainer("dnchdeta",      TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
  coutput13 = mgr->CreateContainer("fullacc",       TTree::Class(), AliAnalysisManager::kOutputContainer ,results.Data());
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
  cout << " === Containers are ready === " << endl;
  return task;
}
