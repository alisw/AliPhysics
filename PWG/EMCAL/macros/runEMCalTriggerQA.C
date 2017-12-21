
void LoadMacros();

AliAnalysisManager* runEMCalTriggerQA(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cRunPeriod     = "LHC16s",                                // set the run period
    const char   *cLocalFiles    = "files_LHC16s_AOD.txt",   // set the local list file
    const UInt_t  iNumEvents     = 5000,                                    // number of events to be analyzed
    UInt_t        iDebugLevel    = 0,                                       // debug level
    const char   *cTaskName      = "EMCalAna",                              // sets name of analysis manager
    const UInt_t  iNumFiles      = 100                                     // number of files analyzed locally
)
{
  TString sRunPeriod(cRunPeriod);
  sRunPeriod.ToLower();

  AliAnalysisTaskEmcalLight::EBeamType_t iBeamType = AliAnalysisTaskEmcalLight::kpp;

  if (sRunPeriod == "lhc10h" || sRunPeriod == "lhc11h" || sRunPeriod == "lhc15o") {
    iBeamType = AliAnalysisTaskEmcalLight::kAA;
  }
  else if (sRunPeriod == "lhc12g" || sRunPeriod == "lhc13b" || sRunPeriod == "lhc13c" ||
      sRunPeriod == "lhc13d" || sRunPeriod == "lhc13e" || sRunPeriod == "lhc13f" ||
      sRunPeriod == "lhc16q" || sRunPeriod == "lhc16r" || sRunPeriod == "lhc16s" ||
      sRunPeriod == "lhc16t") {
    iBeamType = AliAnalysisTaskEmcalLight::kpA;
  }

  enum eDataType { kAod, kEsd };

  eDataType iDataType;
  if (!strcmp(cDataType, "ESD")) {
    iDataType = kEsd;
  }
  else if (!strcmp(cDataType, "AOD")) {
    iDataType = kAod;
  }
  else {
    Printf("Incorrect data type option, check third argument of run macro.");
    Printf("datatype = AOD or ESD");
    return 0;
  }

  Printf("%s analysis chosen.", cDataType);

  TString sLocalFiles(cLocalFiles);
  if (sLocalFiles == "") {
    Printf("You need to provide the list of local files!");
    return 0;
  }
  Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

  LoadMacros();

  if (iDataType == kAod) {
    AliAODInputHandler* pAODHandler = AliAnalysisTaskEmcal::AddAODHandler();
  }
  else {  
    AliESDInputHandler* pESDHandler = AliAnalysisTaskEmcal::AddESDHandler();
  }

  // Physics selection task
  AliPhysicsSelectionTask *pPhysSelTask = AddTaskPhysicsSelection();

  // Multiplicity task
  AliMultSelectionTask *pMultTask = AddTaskMultSelection();

  // CDBconnect task
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
  taskCDB->SetFallBackToRaw(kTRUE);

  // Trigger maker
  AliEmcalTriggerMakerTask *pTriggerMakerTask = AddTaskEmcalTriggerMakerNew("EmcalTriggers");
  pTriggerMakerTask->GetTriggerMaker()->ConfigureForPP2015();
  pTriggerMakerTask->SetUseNewCentralityEstimation(kTRUE);

  // Trigger QA
  AliEmcalTriggerQATask* triggerQAtask = AliEmcalTriggerQATask::AddTaskEmcalTriggerQA_QAtrain(267110);
 
  if (!pMgr->InitAnalysis()) return 0;
  pMgr->PrintStatus();
    
  if (iDebugLevel == 0) {
    pMgr->SetUseProgressBar(kTRUE, 250);
  }
  else {
    pMgr->SetDebugLevel(iDebugLevel);
  }
  
  TFile *pOutFile = new TFile("train.root","RECREATE");
  pOutFile->cd();
  pMgr->Write();
  pOutFile->Close();
  delete pOutFile;

  TChain* pChain = 0;
  if (iDataType == kAod) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
    pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
  }
  else {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    pChain = CreateESDChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
  }

  // start analysis
  Printf("Starting Analysis...");
  pMgr->StartAnalysis("local", pChain, iNumEvents);

  return pMgr;
}

void LoadMacros()
{
  // Aliroot macros
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTriggerMakerNew.C");
}
