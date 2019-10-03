#include "Riostream.h"
void LoadLibraries();
void AddAnalysisTasks(); 
class AliAnalysisAlien;                                                                                                                    
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode);

Int_t runNumbers[5] = {126437};

Bool_t doQAsym        = 1;   // output ok
Bool_t doVZERO        = 1;   // output ok but there is a 2nd file
Bool_t doVertex       = 1;   // output ok
Bool_t doSPD          = 1;   // output ok, needs RP   
Bool_t doTPC          = 1;   // output ok
Bool_t doEventStat    = 0;   // output ok
Bool_t doSDD          = 1;   // outout ok needs RP
Bool_t doSSDdEdx      = 1;   // testing
// new 
Bool_t doTRD          = 1;   // TRD 
Bool_t doITS          = 1;   // ITS
Bool_t doCALO         = 1;   // Calorimeter
Bool_t doMUONTrig     = 1;   // MUON trigger
Bool_t doImpParRes    = 1;   // Impact parameter resolution
Bool_t doMUON         = 1;   // MUON QA

Bool_t doMUONEff      = 0;   // MUON efficiency  NEEDS geometry
Bool_t doV0           = 0;   // V0 recosntruction performance NEEDS MCtruth 

TString     train_name         = "QA";      // QA local folder name
TString     train_tag          = "";        // Train special tag appended to 
                                            // visible name. ("sim", "pp", ...)
               // Name in train page (DON'T CHANGE)
TString     visible_name       = Form("QA$2_$3%s", train_tag.Data()); //# FIXED #
TString     job_comment        = "PWGPP QA train"; // Can add observations here
               // Job tag (DON'T CHANGE)
TString     job_tag            = Form("%s: %s", visible_name.Data(), job_comment.Data());
               // Package versions - Modify as needed
TString     root_version       = "v5-27-06-1";
TString     aliroot_version    = "v4-20-12-AN";
               // Production directory - change as needed for test mode
TString     grid_datadir       = "/alice/sim/LHC10f7";
               // Work directory in GRID (DON'T CHANGE)
TString     grid_workdir       = "/alice/cern.ch/user/a/alidaq/QA/QA$2";
               // Job splitting
Int_t       grid_split         = 20;       // Splitting
               // Debug level
Int_t       debug_level        = 1;        // Debugging
               // Data pattern - change as needed for test mode
TString     data_pattern       = "*ESDs.root";
               // Output directory (DON'T CHANGE)
TString     alien_outdir       = "$1/QA$2";
               // Input collection (production mode)
TString     data_collection    = "$1/qa1.xml";
TString     mergeExcludes      = ""; // Files to be excluded for merging
TString     terminateFiles     = "trending.root"; // Files produced during Terminate

Bool_t useProductionMode       = kTRUE;
Bool_t useMergeViaJDL          = kTRUE;
Bool_t useFastReadOption       = kTRUE;
Bool_t useOverwriteMode        = kFALSE;
Bool_t useDevelopmentVersion   = kFALSE;

void PilotAnalysis(const char *plugin_mode = "full")
{
  TString smode(plugin_mode);
  smode.ToLower();
  if (smode == "test") useProductionMode = kFALSE;
  if (!useProductionMode) {
     TGrid::Connect("alien://");
     if (!gGrid || !gGrid->IsConnected()) {
       ::Error("PilotAnalysis", "No grid connection");
       return;
     }
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
  out << "   doQAsim         = " << doQAsym << ";" << endl;
  out << "   doVZERO         = " << doVZERO << ";" << endl;
  out << "   doVertex        = " << doVertex << ";" << endl;
  out << "   doSPD           = " << doSPD << ";" << endl;
  out << "   doSDD           = " << doSDD << ";" << endl;
  out << "   doSSDdEdx       = " << doSSDdEdx << ";" << endl;
  out << "   doTPC           = " << doTPC << ";" << endl;
  out << "   doTRD           = " << doTRD << ";" << endl;
  out << "   doImpParRes     = " << doImpParRes << ";" << endl;
  out << "   doMUON          = " << doMUON << ";" << endl;
  out << "   doEventStat     = " << doEventStat << ";" << endl;
  out << "}" << endl;
  out.close();
  
  // Load libraries
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TRD");
  LoadLibraries();
  // Create manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("PilotAnalysis", "Production train");
  mgr->SetNSysInfo(100);
  // Input handler
  AliESDInputHandlerRP *esdHandler = new AliESDInputHandlerRP();
  esdHandler->SetReadFriends(kTRUE);
  esdHandler->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdHandler);
  mgr->SetDebugLevel(debug_level);
  mgr->SetSaveCanvases(kFALSE);

  // AnalysisTasks
  AddAnalysisTasks();
  // Grid handler
  AliAnalysisAlien *alienHandler = CreateAlienHandler(plugin_mode);
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
  gSystem->Load("libCORRFW");
  gSystem->Load("libTender");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG0selectors");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWG2");
  gSystem->Load("libPWG2forward");

  if (doCALO) {
     gSystem->Load("libEMCALUtils");
     gSystem->Load("libPWG4PartCorrBase");
     gSystem->Load("libPWG4PartCorrDep");
  }  
  if(doMUONTrig || doAOD) {
     gSystem->Load("libPWGHFbase");
     gSystem->Load("libPWGmuon");
     gSystem->Load("libPWGmuondep");
  }   
}

void AddAnalysisTasks()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetCommonFileName("QAresults.root");
  //
  // Event Statistics (Jan Fiete)
  //

  if (doEventStat) {
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
      if (!terminateFiles.IsNull()) terminateFiles += ",";
      terminateFiles += "event_stat.root";
  }
  // Vertexing (A. Dainese)
  // 
  if (doVertex) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskVertexESD.C");
    AliAnalysisTaskVertexESD* taskvertexesd =  AddTaskVertexESD();
    taskvertexesd->SelectCollisionCandidates();
  }  

  // TPC QA (E. Sicking)
  //
  if (doQAsym) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskQAsym.C");
    AliAnalysisTaskSE * taskqasim = AddTaskQAsym(0);
    taskqasim->SelectCollisionCandidates();
  }  
  //
  // VZERO QA  (C. Cheshkov)
  //
  if (doVZERO) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskVZEROQA.C");
    AliAnalysisTaskSE * taskv0qa = AddTaskVZEROQA(0);
//  taskv0qa->SelectCollisionCandidates();
  }
  //
  // TPC (Jacek Otwinowski & Michael Knichel)
  //
  if (doTPC) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
    AliPerformanceTask *tpcQA = AddTaskPerformanceTPCdEdxQA(kFALSE, kTRUE);   
    tpcQA->SelectCollisionCandidates();
  }  
  //
  // SPD (A. Mastroserio)
  //
  if (doSPD) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskSPDQA.C");
    AliAnalysisTaskSE* taskspdqa = AddTaskSPDQA();
    taskspdqa->SelectCollisionCandidates();
  }  
  //
  // SDD (F. Prino)
  //
  if (doSDD) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddSDDPoints.C");
    AliAnalysisTaskSE* tasksdd = AddSDDPoints();
    tasksdd->SelectCollisionCandidates();
  }
  //
  // SSD dEdx (Marek Chojnacki)
  //
  if (doSSDdEdx) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskdEdxSSDQA.C");
    AliAnalysisTaskSE* taskssddedx = AddTaskdEdxSSDQA();
    taskssddedx->SelectCollisionCandidates();
  }

  //
  // ITS
  //
  if (doITS) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskPerformanceITS.C");
      AliAnalysisTaskITSTrackingCheck *itsQA = AddTaskPerformanceITS(kFALSE);
  }
  //
  // TRD (Alex Bercuci, M. Fasel) 
  //
  if(doTRD) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTrainPerformanceTRD.C");
      // steer individual TRD tasks
      Bool_t 
      doCheckESD(kTRUE),  // AliTRDcheckESD
      doCheckDET(kTRUE),  // AliTRDcheckDET
      doEffic(kTRUE),     // AliTRDefficiency
      doResolution(kTRUE),// AliTRDresolution
      doCheckPID(kTRUE),  // AliTRDcheckPID
      doV0Monitor(kFALSE);// AliTRDv0Monitor
      AddTrainPerformanceTRD(Translate(doCheckESD, doCheckDET, doEffic, doResolution, doCheckPID, doV0Monitor));
  }

  //
  // Calorimetry (Gustavo Conesa)
  //

  if(doCALO) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG4/macros/QA/AddTaskCalorimeterQA.C");
      AliAnalysisTaskParticleCorrelation *taskCaloQA = AddTaskCalorimeterQA("ESD", kTRUE, kFALSE);
      taskCaloQA->SetDebugLevel(0);
  }

  //
  // Muon Trigger
  //
  
  if(doMUONTrig) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskMTRchamberEfficiency.C");
      AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
  }

  //
  // Muon Efficiency
  //

  if(doMUONEff) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG3/muondep/AddTaskMUONTrackingEfficiency.C");
      AliAnalysisTaskMuonTrackingEff *taskMuonTrackEff = AddTaskMUONTrackingEfficiency();
  }
  
  //
  // V0-Decay Reconstruction (Ana Marin)
  // 

  if (doV0) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskV0QA.C");
      AliAnalysisTaskV0QA *taskv0QA = AddTaskV0QA(kFALSE);
  }
  // Impact parameter resolution (xianbao.yuan@pd.infn.it, andrea.dainese@pd.infn.it)
  //
  if (doImpParRes) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskImpParRes.C");
    AliAnalysisTaskSE* taskimpparres= AddTaskImpParRes();
    taskimpparres->SelectCollisionCandidates();
  }  
  // MUON QA (Philippe Pillot)
  //
  if (doMUON) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG3/muon/AddTaskMuonQA.C");
    AliAnalysisTaskSE* taskmuonqa= AddTaskMuonQA(kFALSE, kFALSE);
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
   if (useProductionMode) {
      plugin->SetProductionMode();
      plugin->AddDataFile(data_collection);
   }   
   plugin->SetJobTag(job_tag);
   plugin->SetNtestFiles(1);
   plugin->SetCheckCopy(kFALSE);
   plugin->SetOneStageMerging(kTRUE);
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
   //plugin->SetRunPrefix("000");
//   plugin->SetOutputSingleFolder("output");
   if (!useProductionMode) {
      plugin->SetOutputToRunNo();
      for (Int_t i=0; i<2; i++) {
         if (!runNumbers[i]) break;
         plugin->AddRunNumber(runNumbers[i]);
      }   
   }
// Define alien work directory where all files will be copied. Relative to alien $HOME.
   plugin->SetGridWorkingDir(grid_workdir);
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
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS/TRD");
   
   plugin->SetAdditionalLibs("libTender.so libPWG0base.so libPWG0dep.so libPWG0selectors.so libPWGPP.so libPWG2.so \
                              libPWG2forward.so libEMCALUtils.so libPWG4PartCorrBase.so libPWG4PartCorrDep.so \
                              libPWGHFbase.so libPWGmuon.so libPWGmuondep.so");
     
// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
   plugin->SetDefaultOutputs();
   plugin->SetMaxMergeFiles(20);
   plugin->SetNrunsPerMaster(1);
   
   // Put default output files to archive
   TString listhists = "";
   TString listaods  = "";
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mergeExcludes.IsNull()) plugin->SetMergeExcludes(mergeExcludes);
   if (!terminateFiles.IsNull()) plugin->SetTerminateFiles(terminateFiles);
// Set friends
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro(Form("%s.C", train_name.Data()));
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(grid_split);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally modify the number of replicas
   plugin->SetNumberOfReplicas(4);
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
/*********************************************************
 ***     PROOF MODE SPECIFIC SETTINGS         ************
 *********************************************************/
// Proof cluster
//   plugin->SetProofCluster("alice-caf");
   plugin->SetProofCluster("skaf.saske.sk");
// Dataset to be used   
   plugin->SetProofDataSet("/alice/data/LHC10e_000128175_p1#esdTree");
// May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
   plugin->SetProofReset(0);
// May limit number of workers
   plugin->SetNproofWorkers(20);   
// May use a specific version of root installed in proof
   plugin->SetRootVersionForProof("current_dbg");
// May set the aliroot mode. Check http://aaf.cern.ch/node/83 
   plugin->SetAliRootMode("ALIROOT"); // Loads AF libs by default
// May request ClearPackages (individual ClearPackage not supported)
   plugin->SetClearPackages(kFALSE);
// Plugin test mode works only providing a file containing test file locations
   plugin->SetFileForTestMode(gSystem->ExpandPathName("$ALICE_PHYSICS/PWGPP/PilotTrain/files.txt"));
   return plugin;
}
