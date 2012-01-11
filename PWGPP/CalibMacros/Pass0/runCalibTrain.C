/*
  Template of calibration/filtering macro using ESD:
  - requires AliESDs.root and AliESDfriend.root
  - requires OCDB access (default set to "raw://")
  - requires run number as argument to init OCDB
  - calls LoadLibraries.C, ConfigCalibTrain.C and AddTaskTPCCalib.C macros
  - output AliESDfriends_v1.root with TPC and TRD calibration objects are created

  Example:
  .L $ALICE_ROOT/ANALYSIS/macros/runCalibTrain.C
  runCalibTrain("104892");
*/

void runCalibTrain(TString runNumberString, const char *inFileName = "AliESDs.root", const char *ocdb="raw://")
{
  //
  // macro to run TPC calibration train 
  //
  AliLog::SetGlobalLogLevel(AliLog::kError); 
  gROOT->Macro("LoadLibraries.C");
  gROOT->LoadMacro("ConfigCalibTrain.C");

  // detector tasks
  gROOT->LoadMacro("AddTaskTPCCalib.C");
  gROOT->LoadMacro("AddTaskTRDCalib.C");
  gROOT->LoadMacro("AddTOFAnalysisTaskCalibPass0.C");
  gROOT->LoadMacro("AddTaskT0Calib.C");
  gROOT->LoadMacro("AddTaskMeanVertexCalib.C");
  gROOT->LoadMacro("AddTaskSDDCalib.C"); 

  // switch off debug 
  AliLog::SetClassDebugLevel("AliESDEvent",0);
  
  // steering input chain
  TChain *chain = new TChain("esdTree");
  chain->Add(inFileName);

  // config calibration train
  // setting geometry and B-field from GRP
  Int_t runNumber = runNumberString.Atoi();
  printf("runNumber from runCalibTrain = %d\n",runNumber);
  printf("ocdb from runCalibTrain = %s\n",ocdb);
  ConfigCalibTrain(runNumber, ocdb);
  
  //
  // setup analysis
  //
  AliAnalysisManager *mgr  = new AliAnalysisManager("ESD to ESD", "Analysis Manager");
  // mgr->SetDebugLevel(3);
  
  // Input
  AliESDInputHandler* inpHandler = new AliESDInputHandler();
  inpHandler->SetReadFriends(1);
  mgr->SetInputEventHandler(inpHandler);
  
  // Output
  const char *outFile = "AliESDfriends_v1.root";
  AliESDHandler* esdHandler   = new AliESDHandler();
  mgr->SetOutputEventHandler(esdHandler);
  esdHandler->SetOutputFileName(outFile);
  mgr->SetCommonFileName(outFile);
  //  
  // Detector Tasks
  AliAnalysisTask* tTPC = AddTaskTPCCalib(runNumber);
  AliAnalysisTask* tTRD = AddTaskTRDCalib(runNumber);
  AliTOFAnalysisTaskCalibPass0 *thisTask = AddTOFAnalysisTaskCalibPass0();
  AliAnalysisTask* tT0 = AddTaskT0Calib(runNumber);
  AliMeanVertexCalibTask *tMeanVtx = AddTaskMeanVertexCalib();
  AliAnalysisTaskITSAlignQA *itsAlign = AddTaskSDDCalib();

  // Run the analysis
  if (!mgr->InitAnalysis()) {
    printf("Analysis cannot be started, returning\n");
    return;
  }
  
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain);
  
  return;
}
