#include "Riostream.h"
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
// new 
Bool_t doTRD          = 1;   // TRD 
Bool_t doITS          = 1;   // ITS
Bool_t doCALO         = 1;   // Calorimeter
Bool_t doMUONTrig     = 1;   // MUON trigger
Bool_t doMUONEff      = 0;   // MUON efficiency  NEEDS geometry
Bool_t doV0           = 1;   // V0 recosntruction performance NEEDS MCtruth

TString     train_name         = "QA001_PASS4";
TString     job_tag            = "QA001: PWG1 QA train";
TString     root_version       = "v5-26-00b";
TString     aliroot_version    = "v4-19-04-AN";
TString     grid_datadir       = "/alice/data/2009/LHC09d";
TString     data_pattern       = "*ESDs/pass4/*ESDs.root";
//TString     alien_outdir       = "";
TString     alien_outdir       = "/alice/data/2009/LHC09d/analysis/QA001_PASS4";

Bool_t useProductionMode       = kTRUE;

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
  out << "   doQAsim         = " << doQAsym << ";" << endl;
  out << "   doVZERO         = " << doVZERO << ";" << endl;
  out << "   doVertex        = " << doVertex << ";" << endl;
  out << "   doSPD           = " << doSPD << ";" << endl;
  out << "   doSDD           = " << doSDD << ";" << endl;
  out << "   doFMD           = " << doFMD << ";" << endl;
  out << "   doTPC           = " << doTPC << ";" << endl;
  out << "   doEventStat     = " << doEventStat << ";" << endl;
  out << "}" << endl;
  out.close();
  
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
gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
  gSystem->Load("libPWG1.so");
  gSystem->Load("libPWG2.so");
  gSystem->Load("libPWG3muon.so");
  gSystem->Load("libPWG3muondep.so");
  gSystem->Load("libPWG2forward.so");
  gSystem->Load("libPWG4PartCorrBase.so");
  gSystem->Load("libPWG4PartCorrDep.so");
 
  if (doSPD) {   
    TFile::Cp(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/ITS/AliAnalysisTaskSPD.cxx"), "AliAnalysisTaskSPD.cxx");
    TFile::Cp(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/ITS/AliAnalysisTaskSPD.h"), "AliAnalysisTaskSPD.h");
    gROOT->LoadMacro("AliAnalysisTaskSPD.cxx++g");
  }
  if (doSDD) {  
    TFile::Cp(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/ITS/AliAnalysisTaskSDDRP.cxx"), "AliAnalysisTaskSDDRP.cxx");
    TFile::Cp(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/ITS/AliAnalysisTaskSDDRP.h"), "AliAnalysisTaskSDDRP.h");
//    gROOT->LoadMacro("AliAnalysisTaskSDDRP.cxx++g");
  }  
}

void AddAnalysisTasks()
{
  //
  // Vertexing (A. Dainese)
  // 

  if (doVertex) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskVertexESD.C");
    AliAnalysisTaskVertexESD* taskvertexesd =  AddTaskVertexESD();
    taskvertexesd->SelectCollisionCandidates();
  }  

  //
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
  }

  //
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
  // ITS
  // 
  if (doITS) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceITS.C");
      AliAnalysisTaskITSTrackingCheck *itsQA = AddTaskPerformanceITS(kFALSE);
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
      AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile("event_stat.root");
  }
   

  //
  // TRD (Alex Bercuci, M. Fasel) 
  //
  if(doTRD) {
      AliAnalysisDataContainer *ci[] = {0x0, 0x0, 0x0};
      //
      // Check ESD
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDcheckESD.C++");
      AddTRDcheckESD(mgr);
      //
      // Info top task
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDinfoGen.C++");
      AddTRDinfoGen(mgr, "ALL", 0x0, ci);
      //
      // check DET
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDcheckDET.C++");
      AddTRDcheckDET(mgr, "ALL", ci);
      //
      // check PID (ref maker ???)
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDcheckPID.C++");
      AddTRDcheckPID(mgr, "ALL", ci);
      //
      // Efficiency
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDefficiency.C++");
      AddTRDefficiency(mgr, "ALL", ci);
      //
      // Resolution
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDresolution.C++");      
      AddTRDresolution(mgr, "ALL", ci);
  }

  //
  // Calorimetry (Gustavo Conesa)
  //

  if(doCALO) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/QA/AddTaskCalorimeterQA.C");
      AliAnalysisTaskParticleCorrelation *taskCaloQA = AddTaskCalorimeterQA("ESD", kTRUE, kFALSE);
      taskCaloQA->Dump();
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
   if (alien_outdir.IsNull()) alien_outdir = Form("output_%s",train_name.Data());
   plugin->SetGridOutputDir(alien_outdir);

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
// Set friends
//   if (iAODanalysis && iPWG3d2h) 
//      plugin->SetFriendChainName("AliAOD.VertexingHF.root");
   plugin->SetOutputArchive(outputArchive);
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro(Form("%s.C", train_name.Data()));
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
   plugin->SetJDLName(Form("%s.jdl", train_name.Data()));
// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable(Form("%s.sh", train_name.Data()));
// Optionally modify job price (default 1)
   plugin->SetPrice(1);      
// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");
   plugin->SetExecutableCommand("aliroot -b -q");
   return plugin;
}
