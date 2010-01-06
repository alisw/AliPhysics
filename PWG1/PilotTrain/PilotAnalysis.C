void LoadLibraries();
void AddAnalysisTasks(); 
class AliAnalysisAlien;                                                                                                                    
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode);

Bool_t doQAsym        = 1;   // output ok
Bool_t doVZERO        = 1;   // output ok but there is a 2nd file
Bool_t doVertex       = 1;   // output ok
Bool_t doSPD          = 1;   // output ok, needs RP   
Bool_t doFMD          = 1;   // output ok
Bool_t doTPC          = 1;   // output ok
Bool_t doEventStat    = 1;   // output ok
Bool_t doSDD          = 1;   // outout ok needs RP

void PilotAnalysis(const char *plugin_mode = "full")
{
  TGrid::Connect("alien://");
  if (!gGrid || !gGrid->IsConnected()) {
    ::Error("PilotAnalysis", "No grid connection");
    return;
  }
  // Load libraries
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS");
  LoadLibraries();
  // Create manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("PilotAnalysis", "Production train");
  mgr->SetNSysInfo(1);
  // Input handler
  AliESDInputHandlerRP *esdHandler = new AliESDInputHandlerRP();
  esdHandler->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdHandler);

  mgr->SetDebugLevel(3);

  // AnalysisTasks
  AddAnalysisTasks();
  // Grid handler
  AliAnalysisGrid *alienHandler = CreateAlienHandler(plugin_mode);
  mgr->SetGridHandler(alienHandler);
  if (mgr->InitAnalysis()) {                                                                                                              
    mgr->PrintStatus(); 
    mgr->StartAnalysis("grid");
  }
}

void LoadLibraries()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
  gSystem->Load("libPWG1.so");
  gSystem->Load("libPWG2.so");
  gSystem->Load("libPWG2forward.so");

  if (doSPD) {   
    TFile::Cp(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/PilotTrain/AliAnalysisTaskSPD.cxx"), "AliAnalysisTaskSPD.cxx");
    TFile::Cp(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/PilotTrain/AliAnalysisTaskSPD.h"), "AliAnalysisTaskSPD.h");
    gROOT->LoadMacro("AliAnalysisTaskSPD.cxx++g");
  }
  if (doSDD) {  
    TFile::Cp(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/PilotTrain/AliAnalysisTaskSDDRP.cxx"), "AliAnalysisTaskSDDRP.cxx");
    TFile::Cp(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/PilotTrain/AliAnalysisTaskSDDRP.h"), "AliAnalysisTaskSDDRP.h");
    gROOT->LoadMacro("AliAnalysisTaskSDDRP.cxx++g");
  }  
}

void AddAnalysisTasks()
{
  // Vertexing (A. Dainese)
  // 
  if (doVertex) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskVertexESD.C");
    AliAnalysisTaskVertexESD* taskvertexesd =  AddTaskVertexESD();
    taskvertexesd->SelectCollisionCandidates();
  }  

  // TPC QA (E. Sicking)
  //
  if (doQAsym) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskQAsym.C");
    AliAnalysisTaskSE * taskqasim = AddTaskQAsym();
    taskqasim->SelectCollisionCandidates();
  }  
  //
  // VZERO QA  (C. Cheshkov)
  //
  if (doVZERO) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskVZEROQA.C");
    AliAnalysisTaskSE * taskv0qa = AddTaskVZEROQA(0);
//  taskv0qa->SelectCollisionCandidates();
  }
  // FMD (Hans Hjersing Dalsgaard)
  //
  if (doFMD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskFMD.C");
    AliAnalysisTaskSE* taskfmd = AddTaskFMD();
    taskfmd->SelectCollisionCandidates();
  }  
  //
  // TPC (Jacek Otwinowski)
  //
  if (doTPC) {
    gROOT->LoadMacro("$(ALICE_ROOT)/PWG1/TPC/macros/AddTaskPerformanceTPCQA.C");
    AliPerformanceTask *tpcQA = AddTaskPerformanceTPCQA(kFALSE, kTRUE);
  }  
  //
  // SPD (A. Mastroserio)
  //
  if (doSPD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskSPDQA.C");
    AliAnalysisTaskSE* taskspdqa = AddTaskSPDQA();
    taskspdqa->SelectCollisionCandidates();
  }  
  //
  // SDD (F. Prino)
  //
  if (doSDD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddSDDPoints.C");
    AliAnalysisTaskSE* tasksdd = AddSDDPoints();
    tasksdd->SelectCollisionCandidates();
  }
  //
  // Event Statistics (Jan Fiete)
  //

  if (doEventStat) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
      AliPhysicsSelection* physSel = physSelTask->GetPhysicsSelection();
      physSel->AddBackgroundIdentification(new AliBackgroundSelection());
  }    
}

//______________________________________________________________________________
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode)
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(plugin_mode);
//   if (useProductionMode) plugin->SetProductionMode();
   plugin->SetJobTag("Pilot analysis train");
   plugin->SetNtestFiles(1);
// Set versions of used packages
//   plugin->SetAPIVersion("V2.4");
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-25-04-3");
   plugin->SetAliROOTVersion("v4-18-14-AN-1");
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
   plugin->SetGridDataDir("/alice/data/2009/LHC09d");
// Set data search pattern
   plugin->SetDataPattern("*ESD.tag.root");
// ...then add run numbers to be considered
//   if (!iAODanalysis) plugin->SetRunRange(run_range[0], run_range[1]);
   plugin->SetRunPrefix("000");
   plugin->SetOutputSingleFolder("output");
   plugin->SetOutputToRunNo();
   Int_t run_numbers[30] = {104065, 104155, 104157, 104159, 104160, 104315, 104316, 104320, 104321, 104439, 
                            104792, 104793, 104799, 104800, 104801, 104802, 104803, 104821, 104824, 104825,
                            104841, 104845, 104849, 104852, 104865, 104867, 104876, 104892, 105143, 105160};
   for (Int_t i=0; i<30; i++) {
      plugin->AddRunNumber(run_numbers[i]);
   }   
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir("analysisDATA");
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("pilotAnalysis2");

//   plugin->EnablePackage("");

// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include  -I$ALICE_ROOT/ITS");
   plugin->SetAnalysisSource("AliAnalysisTaskSPD.cxx AliAnalysisTaskSDDRP.cxx");
   plugin->SetAdditionalLibs("libTENDER.so libPWG0base.so libPWG0dep.so libPWG0selectors.so libPWG1.so libPWG2.so libPWG2forward.so AliAnalysisTaskSPD.h AliAnalysisTaskSPD.cxx AliAnalysisTaskSDDRP.h AliAnalysisTaskSDDRP.cxx");
     
// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetDefaultOutputs();
//   plugin->SetMergeExcludes(mergeExclude);
   plugin->SetMaxMergeFiles(100);
   plugin->SetNrunsPerMaster(1);
// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:AliAOD.root,AOD.tag.root@ALICE::NIHAM::File");
   
   // Put default output files to archive
   TString listhists = "";
   TString listaods  = "";
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   TIter next(mgr->GetOutputs());
   AliAnalysisDataContainer *output;
   while ((output=(AliAnalysisDataContainer*)next())) {
      const char *filename = output->GetFileName();
      if (!(strcmp(filename, "default"))) {
         if (!mgr->GetOutputEventHandler()) continue;
         filename = mgr->GetOutputEventHandler()->GetOutputFileName();
         if (listaods.Length()) listaods += ",";
         listaods += filename;
      } else {   
         if (listhists.Contains(filename)) continue;
         if (listhists.Length()) listhists += ",";
         listhists += filename;
      }   
   }
   if (mgr->GetExtraFiles().Length()) {
      if (listaods.Length()) listaods += ",";
      listaods += mgr->GetExtraFiles();
      listaods.ReplaceAll(" ", ",");
   }
   if (listhists.Length()) listhists = Form("hist_archive.zip:%s", listhists.Data());
   if (listaods.Length())  listaods  = Form("aod_archive.zip:%s", listaods.Data());
   if (!listhists.Length()) {
      ::Fatal("AnalysisTrainNew", "No task output !");
   }
   TString outputArchive = "log_archive.zip:stdout,stderr";
   if (listaods.Length()) {
      outputArchive += " ";
      outputArchive += listaods;
   }   
   if (listhists.Length()) {
      outputArchive += " ";
      outputArchive += listhists;
   }   
// Set friends
//   if (iAODanalysis && iPWG3d2h) 
//      plugin->SetFriendChainName("AliAOD.VertexingHF.root");
   plugin->SetOutputArchive(outputArchive);
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro("pilotAnalysis002.C");
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(1000);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(70000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName("pilotAnalysis002.jdl");
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable("pilotAnalysis002.sh");
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   plugin->SetExecutableCommand("aliroot -b -q");
   return plugin;
}
