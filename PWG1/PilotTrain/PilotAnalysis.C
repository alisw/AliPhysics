#include "Riostream.h"
void LoadLibraries();
void AddAnalysisTasks(); 
class AliAnalysisAlien;                                                                                                                    
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode);

Int_t runNumbers[5] = {119934};

Bool_t doAOD          = 0;   
Bool_t doQAsym        = 1;   // output ok
Bool_t doVZERO        = 1;   // output ok but there is a 2nd file
Bool_t doVertex       = 1;   // output ok
Bool_t doSPD          = 1;   // output ok, needs RP   
Bool_t doFMD          = 1;   // output ok
Bool_t doTPC          = 1;   // output ok
Bool_t doEventStat    = 1;   // output ok
Bool_t doSDD          = 1;   // outout ok needs RP
Bool_t doSSDdEdx      = 1;   // testing
// new 
Bool_t doTRD          = 1;   // TRD 
Bool_t doITS          = 1;   // ITS
Bool_t doCALO         = 1;   // Calorimeter
Bool_t doMUONTrig     = 1;   // MUON trigger
Bool_t doMUONEff      = 0;   // MUON efficiency  NEEDS geometry
Bool_t doV0           = 0;   // V0 recosntruction performance NEEDS MCtruth 

TString     train_name         = "QAtest";
//TString     train_name         = "TR019_PASS6";
TString     job_tag            = "QA4_LHC10c: PWG1 QA train";
//TString     job_tag            = "TR019: LHC09d-Pass6 ESD filtering w. PhysSelection -> AOD (including muon deltas)";
TString     root_version       = "v5-26-00b-5";
TString     aliroot_version    = "v4-19-13-AN";
TString     grid_datadir       = "/alice/data/2010/LHC10c";
//TString     grid_datadir       = "/alice/data/2009/LHC09d";
TString     data_pattern       = "*ESDs/pass1/*ESDs.root";
TString     alien_outdir       = "";
//TString     alien_outdir       = "/alice/cern.ch/user/m/mgheata/analysisDATA/output_QA007_PASS1_7TeV/000114917";
//TString     alien_outdir       = "/alice/data/2009/LHC09d/analysis/PASS6/AOD";
TString     mergeExcludes;

Bool_t useProductionMode       = kTRUE;
Bool_t useMergeViaJDL          = kTRUE;
Bool_t useFastReadOption       = kTRUE;
Bool_t useOverwriteMode        = kTRUE;
Bool_t useDevelopmentVersion   = kFALSE;

void PilotAnalysis(const char *plugin_mode = "full")
{
  TGrid::Connect("alien://");
  if (!gGrid || !gGrid->IsConnected()) {
    ::Error("PilotAnalysis", "No grid connection");
    return;
  }
  // Write configuration
  TString cdir = gSystem->WorkingDirectory();
  gSystem->MakeDirectory(train_name);
  gSystem->ChangeDirectory(train_name);
  ofstream out;
  out.open(Form("%sConfig.C",train_name.Data()), ios::out);
  out << "{" << endl;
  out << "   train_name      = " << "\"" << train_name.Data() << "\";" << endl;
  out << "   root_version    = " << "\"" << root_version.Data() << "\";" << endl;
  out << "   aliroot_version = " << "\"" << aliroot_version.Data() << "\";" << endl;
  out << "   grid_datadir   = " << "\"" << grid_datadir.Data() << "\";" << endl;
  if (!alien_outdir.Length()) alien_outdir = Form("output_%s",train_name.Data());
  out << "   alien_outdir    = " << "\"" << alien_outdir.Data() << "\";" << endl;
  out << "   doAOD           = " << doAOD << ";" << endl;
  out << "   doQAsim         = " << doQAsym << ";" << endl;
  out << "   doVZERO         = " << doVZERO << ";" << endl;
  out << "   doVertex        = " << doVertex << ";" << endl;
  out << "   doSPD           = " << doSPD << ";" << endl;
  out << "   doSDD           = " << doSDD << ";" << endl;
  out << "   doSSDdEdx       = " << doSSDdEdx << ";" << endl;
  out << "   doFMD           = " << doFMD << ";" << endl;
  out << "   doTPC           = " << doTPC << ";" << endl;
  out << "   doEventStat     = " << doEventStat << ";" << endl;
  out << "}" << endl;
  out.close();
  
  // Load libraries
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD");
  LoadLibraries();
  // Create manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("PilotAnalysis", "Production train");
  mgr->SetNSysInfo(100);
  // Input handler
  AliESDInputHandlerRP *esdHandler = new AliESDInputHandlerRP();
  esdHandler->SetReadFriends(kTRUE);
  esdHandler->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdHandler);
  if (doAOD) {
     // AOD output handler
     AliAODHandler* aodHandler   = new AliAODHandler();
     aodHandler->SetOutputFileName("AliAOD.root");
     if (!mergeExcludes.IsNull()) mergeExcludes += " ";
     mergeExcludes += "AliAOD.root";
     mgr->SetOutputEventHandler(aodHandler);
  }   

  mgr->SetDebugLevel(1);
  mgr->SetSaveCanvases(kTRUE);
  
  // AnalysisTasks
  AddAnalysisTasks();
  // Grid handler
  AliAnalysisAlien *alienHandler = CreateAlienHandler(plugin_mode);
  mgr->SetGridHandler(alienHandler);
  if (mgr->InitAnalysis()) {                                                                                                              
    mgr->PrintStatus(); 
    mgr->StartAnalysis("grid");
    TString alien_workdir = gGrid->GetHomeDirectory();
    alien_workdir += "analysisDATA";
    TString configName = Form("%s/%sConfig.C", alien_workdir.Data(), train_name.Data());
    if (strcmp(plugin_mode, "test")) {
      printf("=== Registering configuration file <%s>===\n", configName.Data());
      if (AliAnalysisAlien::FileExists(configName.Data())) gGrid->Rm(configName.Data());                                                     
      TFile::Cp(Form("file:%sConfig.C",train_name.Data()), Form("alien://%s", configName.Data()));  
    }  
  }
}

void LoadLibraries()
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
  gSystem->Load("libPWG1.so");
  gSystem->Load("libPWG2.so");
  gSystem->Load("libPWG2forward.so");

  if (doCALO) {
     gSystem->Load("libEMCALUtils");
     gSystem->Load("libPWG4PartCorrBase");
     gSystem->Load("libPWG4PartCorrDep");
  }  
  if(doMUONTrig || doAOD) {
     gSystem->Load("libPWG3base");
     gSystem->Load("libPWG3muon");
     gSystem->Load("libPWG3muondep");
  }   
}

void AddAnalysisTasks()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetCommonFileName("QAresults.root");
  // AOD creation with collision events
  if (doAOD) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C");
    mgr->RegisterExtraFile("AliAOD.Muons.root");
    mgr->RegisterExtraFile("AliAOD.Dimuons.root");
    AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(kFALSE, kTRUE, kTRUE, doEventStat);
  }   
  //
  // Event Statistics (Jan Fiete)
  //

  if (doEventStat) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
      mgr->RegisterExtraFile("event_stat.root");
  }
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
    AliAnalysisTaskSE * taskqasim = AddTaskQAsym(0);
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
  // SSD dEdx (Marek Chojnacki)
  //
  if (doSSDdEdx) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskdEdxSSDQA.C");
    AliAnalysisTaskSE* taskssddedx = AddTaskdEdxSSDQA();
    taskssddedx->SelectCollisionCandidates();
  }

  //
  // ITS
  //
  if (doITS) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceITS.C");
      AliAnalysisTaskITSTrackingCheck *itsQA = AddTaskPerformanceITS(kFALSE);
  }
  //
  // TRD (Alex Bercuci, M. Fasel) 
  //
  if(doTRD) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTrainPerformanceTRD.C");
      AddTrainPerformanceTRD("ALL");
  }

  //
  // Calorimetry (Gustavo Conesa)
  //

  if(doCALO) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/QA/AddTaskCalorimeterQA.C");
      AliAnalysisTaskParticleCorrelation *taskCaloQA = AddTaskCalorimeterQA("ESD", kTRUE, kFALSE);
      taskCaloQA->SetDebugLevel(0);
  }

  //
  // Muon Trigger
  //
  
  if(doMUONTrig) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskMTRchamberEfficiency.C");
      AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
  }

  //
  // Muon Efficiency
  //

  if(doMUONEff) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMUONTrackingEfficiency.C");
      AliAnalysisTaskMuonTrackingEff *taskMuonTrackEff = AddTaskMUONTrackingEfficiency();
  }
  
  //
  // V0-Decay Reconstruction (Ana Marin)
  // 

  if (doV0) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskV0QA.C");
      AliAnalysisTaskV0QA *taskv0QA = AddTaskV0QA(kFALSE);
  }
  // FMD (Hans Hjersing Dalsgaard)
  //
  if (doFMD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskFMD.C");
    AliAnalysisTaskSE* taskfmd = AddTaskFMD();
    taskfmd->SelectCollisionCandidates();
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
   if (useProductionMode) plugin->SetProductionMode();
   plugin->SetJobTag(job_tag);
   plugin->SetNtestFiles(1);
// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion(root_version);
   plugin->SetAliROOTVersion(aliroot_version);
// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
   plugin->SetGridDataDir(grid_datadir);
// Set data search pattern
   plugin->SetDataPattern(data_pattern);
// ...then add run numbers to be considered
//   if (!iAODanalysis) plugin->SetRunRange(run_range[0], run_range[1]);
   plugin->SetRunPrefix("000");
//   plugin->SetOutputSingleFolder("output");
   plugin->SetOutputToRunNo();
//   Int_t run_numbers[30] = {104065, 104155, 104157, 104159, 104160, 104315, 104316, 104320, 104321, 104439, 
//                            104792, 104793, 104799, 104800, 104801, 104802, 104803, 104821, 104824, 104825,
//                            104841, 104845, 104849, 104852, 104865, 104867, 104876, 104892, 105143, 105160};
//   Int_t run_numbers[8] = {114785, 114778, 114757, 114753, 114745, 114744, 114743, 114737};
//   Int_t run_numbers[2] = {114785, 114917};
   for (Int_t i=0; i<2; i++) {
      if (!runNumbers[i]) break;
      plugin->AddRunNumber(runNumbers[i]);
   }   
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir(train_name);
// Declare alien output directory. Relative to working directory.
   if (alien_outdir.IsNull()) alien_outdir = Form("output_%s",train_name.Data());
   plugin->SetGridOutputDir(alien_outdir);

   if (useDevelopmentVersion) {
     plugin->EnablePackage("STEERBase");
     plugin->EnablePackage("ESD");
     plugin->EnablePackage("AOD");
     plugin->EnablePackage("ANALYSIS");
     plugin->EnablePackage("ANALYSISalice");
   }

// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD");
   
   plugin->SetAdditionalLibs("libTENDER.so libPWG0base.so libPWG0dep.so libPWG0selectors.so libPWG1.so libPWG2.so \
                              libPWG2forward.so libEMCALUtils.so libPWG4PartCorrBase.so libPWG4PartCorrDep.so \
                              libPWG3base.so libPWG3muon.so libPWG3muondep.so");
     
// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetDefaultOutputs();
//   plugin->SetMergeExcludes(mergeExclude);
   plugin->SetMaxMergeFiles(20);
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
      if (listhists.Length()) listhists += ",";
      listhists += mgr->GetExtraFiles();
      listhists.ReplaceAll(" ", ",");
   }
   if (listhists.Length()) listhists = Form("hist_archive.zip:%s", listhists.Data());
   if (listaods.Length())  listaods  = Form("aod_archive.zip:%s", listaods.Data());
   if (!listhists.Length()) {
      ::Fatal("AnalysisTrainNew", "No task output !");
   }
   TString outputArchive = "log_archive.zip:stdout,stderr@disk=4";
   if (listaods.Length()) {
      outputArchive += " ";
      outputArchive += listaods;
      outputArchive += "@disk=4";
   }   
   if (listhists.Length()) {
      outputArchive += " ";
      outputArchive += listhists;
      outputArchive += "@disk=4";
   }   
   if (!mergeExcludes.IsNull()) plugin->SetMergeExcludes(mergeExcludes);
// Set friends
//   if (iAODanalysis && iPWG3d2h) 
//      plugin->SetFriendChainName("AliAOD.VertexingHF.root");
   plugin->SetOutputArchive(outputArchive);
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro(Form("%s.C", train_name.Data()));
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(1);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally resubmit threshold.
//   plugin->SetMasterResubmitThreshold(90);
// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(70000);
// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");
// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName(Form("%s.jdl", train_name.Data()));
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable(Form("%s.sh", train_name.Data()));
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   plugin->SetExecutableCommand("aliroot -b -q");
// Merge via JDL
   plugin->SetMergeViaJDL(useMergeViaJDL);
// Use fastread option
   plugin->SetFastReadOption(useFastReadOption);
// UseOverwrite mode
   plugin->SetOverwriteMode(useOverwriteMode);   
   return plugin;
}
