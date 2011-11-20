#include "Riostream.h"
void LoadLibraries();
void AddAnalysisTasks(); 
class AliAnalysisAlien;                                                                                                                    
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode);

// Collision type: 0 = p-p   1 = Pb-Pb
Int_t  iCollisionType = 0;
// Trigger mask.
// Still need to change:
UInt_t kTriggerInt = AliVEvent::kAnyINT;
UInt_t kTriggerMuonAll = AliVEvent::kMUL7 | AliVEvent::kMUSH7 | AliVEvent::kMUU7 | AliVEvent::kMUS7;
UInt_t kTriggerMuonBarell = AliVEvent::kMUU7;
UInt_t kTriggerEMC   = AliVEvent::kEMC7;
UInt_t kTriggerHM   = AliVEvent::kHighMult;
// Main trigger mask used:
UInt_t kTriggerMask = kTriggerInt;

Int_t runNumbers[5] = {154780};

Bool_t doCDBconnect   = 1;
Bool_t doEventStat    = 1;
Bool_t doCentrality   = 0;
Bool_t doQAsym        = 1;
Bool_t doVZERO        = 1;   // there is a 2nd file
Bool_t doVertex       = 1;
Bool_t doSPD          = 1;   // needs RP   
Bool_t doTPC          = 1;
Bool_t doSDD          = 1;   // needs RP
Bool_t doSSDdEdx      = 1;

Bool_t doTRD          = 1;
Bool_t doITS          = 1;
Bool_t doITSsaTracks  = 1; 
Bool_t doITSalign     = 1;  
Bool_t doCALO         = 1;
Bool_t doMUONTrig     = 1;
Bool_t doImpParRes    = 1;
Bool_t doMUON         = 1;
Bool_t doTOF          = 1;
Bool_t doHMPID        = 1;
Bool_t doT0           = 1;
Bool_t doZDC          = 1;
Bool_t doPIDResponse  = 1;
Bool_t doPIDqa        = 1; //new
Bool_t doFMD          = 1; // new

Bool_t doMUONEff      = 0;   // NEEDS geometry
Bool_t doV0           = 0;   // NEEDS MCtruth 

TString     train_name         = "QA";      // QA local folder name
TString     train_tag          = (iCollisionType)?"_Pb-Pb":"_p-p";        // Train special tag appended to 
                                            // visible name. ("sim", "pp", ...)
               // Name in train page (DON'T CHANGE)
TString     visible_name       = Form("QA$2_$3%s", train_tag.Data()); //# FIXED #
TString     job_comment        = "PWG1 QA kAnyInt, QAsym(add kMUU7 and kEMC7) CALO (add kEMC7)  triggers"; // Can add observations here
               // Job tag (DON'T CHANGE)
TString     job_tag            = Form("%s: %s", visible_name.Data(), job_comment.Data());
               // Package versions - Modify as needed
TString     root_version       = "v5-28-00e-1";
TString     aliroot_version    = "v4-21-29-AN";
               // Production directory - change as needed for test mode
TString     grid_datadir       = "/alice/data/2011/LHC11c";
               // Work directory in GRID (DON'T CHANGE)
TString     grid_workdir       = "/alice/cern.ch/user/a/alidaq/QA/QA$2";
               // Job splitting
Int_t       grid_split         = 20;       // Splitting
               // Debug level
Int_t       debug_level        = 1;        // Debugging
               // File merging
Int_t       maxMergeFiles      = 10;       // Max files to merge in a chunk
               // Data pattern - change as needed for test mode
TString     data_pattern       = "*ESDs/pass1/*ESDs.root";
               // Output directory (DON'T CHANGE)
TString     alien_outdir       = "$1/QA$2";
               // Input collection (production mode)
TString     data_collection    = "$1/qa1.xml";
TString     mergeExcludes      = ""; // Files to be excluded for merging
TString     mergeDirName       = "QA$2";
TString     terminateFiles     = "trending.root"; // Files produced during Terminate

Bool_t useProductionMode       = kTRUE;
Bool_t useMergeViaJDL          = kTRUE;
Bool_t useFastReadOption       = kFALSE;
Bool_t useOverwriteMode        = kTRUE;
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
  out << "   doITS           = " << doITS << ";" << endl;
  out << "   doITSsaTracks   = " << doITSsaTracks << ";" << endl;
  out << "   doITSalign      = " << doITSalign << ";" << endl;
  out << "   doZDC           = " << doZDC << ";" << endl;
  out << "   doImpParRes     = " << doImpParRes << ";" << endl;
  out << "   doMUON          = " << doMUON << ";" << endl;
  out << "   doTOF           = " << doTOF << ";" << endl;
  out << "   doHMPID         = " << doHMPID << ";" << endl;
  out << "   doZDC           = " << doZDC << ";" << endl;
  out << "   doT0            = " << doT0 << ";" << endl;
  out << "   doPIDResponse   = " << doPIDResponse << ";" << endl;
  out << "   doPIDqa         = " << doPIDqa << ";" << endl;
  out << "   doFMD           = " << doFMD << ";" << endl;
  out << "   doEventStat     = " << doEventStat << ";" << endl;
  if (iCollisionType) out << "   doCentrality    = " << doCentrality << ";" << endl;
  out << "}" << endl;
  out.close();
  
  // Load libraries
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWG1");
  LoadLibraries();
  // Create manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("PilotAnalysis", "Production train");
  if (!strcmp(plugin_mode,"test")) mgr->SetNSysInfo(100);
  // Input handler
  AliESDInputHandlerRP *esdHandler = new AliESDInputHandlerRP();
  esdHandler->SetReadFriends(kTRUE);
  esdHandler->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdHandler);
  mgr->SetDebugLevel(debug_level);
  
  // AnalysisTasks
  AddAnalysisTasks();
  // Grid handler
  AliAnalysisAlien *alienHandler = CreateAlienHandler(plugin_mode);
  mgr->SetGridHandler(alienHandler);
  if (mgr->InitAnalysis()) {                                                                                                              
    mgr->PrintStatus(); 
    if (!strcmp(plugin_mode, "local")) mgr->StartAnalysis("local");
    else mgr->StartAnalysis("grid");
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

  if (doCALO) {
     gSystem->Load("libEMCALUtils");
     gSystem->Load("libPHOSUtils");
     gSystem->Load("libPWG4PartCorrBase");
     gSystem->Load("libPWG4PartCorrDep");
  }  
  if(doMUON || doMUONTrig) {
     gSystem->Load("libPWG3base");
     gSystem->Load("libPWG3muon");
     gSystem->Load("libPWG3muondep");
  }
  if (doFMD) {
     gSystem->Load("libPWG2forward2");
  }      
}

void AddAnalysisTasks()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetCommonFileName("QAresults.root");
  // Statistics task
  mgr->AddStatisticsTask(kTriggerMask);
  //
  // CDB connection
  //
  if (doCDBconnect) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
    if (!taskCDB) return;
    taskCDB->SetRunNumber(runNumbers[0]);
  }    
  
  //
  // Event Statistics (Jan Fiete)
  //
  if (doEventStat) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE /*MC*/);
//      mgr->RegisterExtraFile("event_stat.root");
      if (!terminateFiles.IsNull()) terminateFiles += ",";
      terminateFiles += "event_stat.root";
  }
  
  //
  // Centrality (A. Toia)
  //
  if (doCentrality) {
     if (!iCollisionType) {
        printf("Disabling centrality task for p-p\n");
        doCentrality = kFALSE;
     } else {           
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
        AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
     }   
  }   
  
  // Vertexing (A. Dainese)
  // 
  if (doVertex) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskVertexESD.C");
    AliAnalysisTaskVertexESD* taskvertexesd =  AddTaskVertexESD(kFALSE, kTriggerMask);
    taskvertexesd->SelectCollisionCandidates(kTriggerMask);
  }  

  // TPC QA (E. Sicking)
  //
  if (doQAsym) {
  // offline trigger in AddTask
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskQAsym.C");
    AliAnalysisTaskSE * taskqasim = AddTaskQAsym(0, kTriggerMask, kTriggerHM, kTriggerEMC, kTriggerMuonBarell);
    // taskqasim->SelectCollisionCandidates(); // Set by AddTask
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
  // TPC (Jacek Otwinowski & Michael Knichel & Weilin Yu)
  //
  //
  // Optionally MC information can be used by setting the 1st argument to true
  // Optionally friends information can be switched off by setting the 2st argument 
  // to false
  // Optionally highMult axis can be used by setting the 3st argument to true (for PbPb)
  if (doTPC) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
    // low multiplicity (pp) 
    //AliPerformanceTask *tpcQA = AddTaskPerformanceTPCdEdxQA(kFALSE, kTRUE, kFALSE);
    // high multiplicity (Pb-Pb)
    AliPerformanceTask *tpcQA = AddTaskPerformanceTPCdEdxQA(kFALSE, kTRUE, kTRUE);
    tpcQA->SelectCollisionCandidates(kTriggerMask);
  }  
  //
  // SPD (A. Mastroserio)
  //
  if (doSPD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskSPDQA.C");
    AliAnalysisTaskSE* taskspdqa = AddTaskSPDQA();
    taskspdqa->SelectCollisionCandidates(kTriggerMask);
  }  
  //
  // SDD (F. Prino)
  //
  if (doSDD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddSDDPoints.C");
    AliAnalysisTaskSE* tasksdd = AddSDDPoints();
    tasksdd->SelectCollisionCandidates(kTriggerMask);
  }
  //
  // SSD dEdx (Marek Chojnacki)
  //
  if (doSSDdEdx) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskdEdxSSDQA.C");
    AliAnalysisTaskSE* taskssddedx = AddTaskdEdxSSDQA();
    taskssddedx->SelectCollisionCandidates(kTriggerMask);
  }

  //
  // ITS
  //
  if (doITS) {
  // hardcoded non-zero trigger mask
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceITS.C");
      AliAnalysisTaskITSTrackingCheck *itsQA = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent0010 = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent3050 = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent6080 = 0;
      if(iCollisionType==0) {
        itsQA = AddTaskPerformanceITS(kFALSE);
      } else {
        itsQA = AddTaskPerformanceITS(kFALSE);
        itsQACent0010 = AddTaskPerformanceITS(kFALSE,kFALSE,kFALSE,3500,10000);
        itsQACent3050 = AddTaskPerformanceITS(kFALSE,kFALSE,kFALSE,590,1570);
        itsQACent6080 = AddTaskPerformanceITS(kFALSE,kFALSE,kFALSE,70,310);
      }
  }
  //
  // ITS saTracks, align (F.Prino)
  //
  if (doITSsaTracks) {
  // offline trigger in AddTask
     gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskITSsaTracks.C");
     AliAnalysisTaskITSsaTracks *itssaTracks = AddTaskITSsaTracks(kFALSE,kFALSE);
     itssaTracks->SelectCollisionCandidates(kTriggerMask);
  }   
  if (doITSalign) {
  // no offline trigger selection
     gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskITSAlign.C");
     AliAnalysisTaskITSAlignQA *itsAlign = AddTaskITSAlign(0,2011);
  }   
  //
  // TRD (Alex Bercuci, M. Fasel) 
  //
  if(doTRD) {
  // no offline trigger selection
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTrainPerformanceTRD.C");
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
  // ZDC (Chiara Oppedisano) 
  //
  if(doZDC) {
  // hardcoded kMB trigger mask
     gROOT->LoadMacro("$ALICE_ROOT/PWG1/ZDC/AddTaskZDCQA.C");
     AliAnalysisTaskSE *taskZDC = AddTaskZDCQA();
     taskZDC->SelectCollisionCandidates(kTriggerMask);
  }   
  //
  // Calorimetry (Gustavo Conesa)
  //

  if(doCALO) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/QA/AddTaskCalorimeterQA.C");
      AliAnalysisTaskParticleCorrelation *taskCaloQA = AddTaskCalorimeterQA("ESD", 2011, kFALSE, kFALSE);
      taskCaloQA->SetDebugLevel(0);
      // offline mask set in AddTask to kMB
      taskCaloQA->SelectCollisionCandidates(kTriggerMask);
      // Add a new calo task with EMC1 trigger only
      taskCaloQA = AddTaskCalorimeterQA("ESD", 2011, kFALSE, kFALSE, "", "EMC7");
      taskCaloQA->SetDebugLevel(0);
      taskCaloQA->SelectCollisionCandidates(kTriggerEMC);
  }

  //
  // Muon Trigger
  //
  
  if(doMUONTrig) {
  // no offline trigger selection
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskMTRchamberEfficiency.C");
      AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
  }

  //
  // Muon Efficiency (not used)
  //

  if(doMUONEff) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMUONTrackingEfficiency.C");
      AliAnalysisTaskMuonTrackingEff *taskMuonTrackEff = AddTaskMUONTrackingEfficiency();
  }
  
  //
  // V0-Decay Reconstruction (Ana Marin) (not used)
  // 

  if (doV0) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskV0QA.C");
      AliAnalysisTaskV0QA *taskv0QA = AddTaskV0QA(kFALSE);
  }
  //
  // Impact parameter resolution (xianbao.yuan@pd.infn.it, andrea.dainese@pd.infn.it)
  //
  if (doImpParRes) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskImpParRes.C");
    AliAnalysisTaskSE* taskimpparres=0;
    if(iCollisionType==0) {
       taskimpparres= AddTaskImpParRes();
    } else {
       taskimpparres= AddTaskImpParRes(kFALSE,-1,kFALSE,kFALSE);
    }
    taskimpparres->SelectCollisionCandidates(kTriggerMask);
  }  
  //
  // MUON QA (Philippe Pillot)
  //
  if (doMUON) {
  // trigger analysis internal
    gROOT->LoadMacro("$ALICE_ROOT/PWG3/muon/AddTaskMuonQA.C");
    AliAnalysisTaskSE* taskmuonqa= AddTaskMuonQA();
  }  
  //
  // TOF (Francesca Bellini)
  //
  if (doTOF) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/TOF/AddTaskTOFQA.C");
    AliAnalysisTaskTOFqa *tofQA = AddTaskTOFQA();
    tofQA->SelectCollisionCandidates(kTriggerMask);
  } 
   //
  // PIDResponse(JENS)
  //
  if (doPIDqa) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"); 
    AliAnalysisTaskPIDResponse *PIDResponse = AddTaskPIDResponse();
    PIDResponse->SelectCollisionCandidates(kTriggerMask);
  }  

  //
  // PIDqa(JENS)
  //
  if (doPIDqa) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa *PIDQA = AddTaskPIDqa();
    PIDQA->SelectCollisionCandidates(kTriggerMask);
  }  
 
  //
  // HMPID QA (Giacomo Volpe)
  //
  if (doHMPID) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/HMPID/AddTaskHmpidQA.C");
    AliAnalysisTaskSE* taskhmpidqa= AddTaskHmpidQA(kFALSE);
      // offline mask set in AddTask to kMB
    taskhmpidqa->SelectCollisionCandidates(kTriggerMask);
  }      
  // T0 QA (Alla Mayevskaya)
  if (doT0) {
  // no offline trigger selection
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/T0/AddTaskT0QA.C");
    AliT0AnalysisTaskQA* taskt0qa= AddTaskT0QA();
    taskt0qa->SelectCollisionCandidates(kTriggerMask);
  }      
  // FMD QA (Christian Holm Christiansen)
  if (doFMD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskForwardQA.C");
    // Parameters: usemc, usecentrality
    AliAnalysisTaskSE *forwardQA = (AliAnalysisTaskSE *)AddTaskForwardQA(kFALSE, kFALSE);
    // No offline trigger config. needed (see #84077)
  }
}

//______________________________________________________________________________
AliAnalysisAlien* CreateAlienHandler(const char *plugin_mode)
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   plugin->SetRunMode(plugin_mode);
   if (useProductionMode) {
      plugin->SetProductionMode();
      plugin->AddDataFile(data_collection);
   }   
   plugin->SetJobTag(job_tag);
   plugin->SetNtestFiles(2);
   plugin->SetCheckCopy(kFALSE);
   plugin->SetMergeDirName(mergeDirName);
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
//   plugin->SetOutputSingleFolder("outpu$ALICE_ROOT/PWG1/T0/Addt");
   if (!useProductionMode) {
      plugin->SetRunPrefix("000");
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
   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD");
   
   plugin->SetAdditionalLibs("libCORRFW.so libTENDER.so libPWG0base.so libPWG0dep.so libPWG0selectors.so libPWG1.so \
                              libEMCALUtils.so libPHOSUtils.so libPWG4PartCorrBase.so libPWG4PartCorrDep.so \
                              libPWG3base.so libPWG3muon.so libPWG3muondep.so libPWG2forward2.so");
     
// Declare the output file names separated by blancs.
   plugin->SetDefaultOutputs();
   plugin->SetMaxMergeFiles(maxMergeFiles);
   plugin->SetNrunsPerMaster(1);
   
   // Put default output files to archive
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   TIter next(mgr->GetOutputs());
   AliAnalysisDataContainer *output;
   if (!mergeExcludes.IsNull()) plugin->SetMergeExcludes(mergeExcludes);
   if (!terminateFiles.IsNull()) plugin->SetTerminateFiles(terminateFiles);
// Set friends
// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro(Form("%s.C", train_name.Data()));
// Optionally set a name for the generated validation script
   plugin->SetValidationScript("validation.sh");
// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(grid_split);
// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
//   plugin->SetMaxInitFailed(5);
// Optionally modify the number of replicas
   plugin->SetNumberOfReplicas(5);
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
   plugin->SetNproofWorkers(1);   
// May use a specific version of root installed in proof
   plugin->SetRootVersionForProof("current_dbg");
// May set the aliroot mode. Check http://aaf.cern.ch/node/83 
   plugin->SetAliRootMode("ALIROOT"); // Loads AF libs by default
// May request ClearPackages (individual ClearPackage not supported)
   plugin->SetClearPackages(kFALSE);
// Plugin test mode works only providing a file containing test file locations
   plugin->SetFileForTestMode(gSystem->ExpandPathName("$ALICE_ROOT/PWG1/PilotTrain/files.txt"));
   return plugin;
}
