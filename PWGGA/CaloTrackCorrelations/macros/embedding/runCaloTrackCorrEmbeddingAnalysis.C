/** @file runCaloTrackCorrEmbeddingAnalysis.C
 * @brief Example embedding analysis run macro
 *
 * @ingroup CaloTrackCorrMacros
 * Example embedding framework run macro with CaloTrackCorrelation example analysis. 
 *
 *http://alidoc.cern.ch/AliPhysics/master/READMEemcEmbedding.html
 *
 * Such collections should be available in the case of embedding MC productions. Embedding
 * a real pp data set (which would not have embedded particle level information available)
 * must be setup with more care. 
 *
 * Data files location should be in txt file "aodFiles.txt"
 * MC files to be embedded in data should be located in "aodFilesEmbed.txt"
 * Some options can be modified in the section under the comment "Configuration options".
 * Adaptation of macro runEmbeddingAnalysis.C from Raymond Ehlers.
 * 
 * @author Gustavo Conesa Balbastre <gustavo.conesa.balbastre@cern.ch>, LPSC-Grenoble
 * @date Mar 11, 2020
 */

#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "AliLog.h"
#include "TChain.h"
#include "TSystemDirectory.h"
#include "Riostream.h"

class AliESDInputHandler;
class AliAODInputHandler;
class AliVEvent;
class AliAnalysisGrid;
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAnalysisAlien.h"
#include "AliPhysicsSelectionTask.h"
#include "AliMultSelectionTask.h"
#include "AliCentralitySelectionTask.h"
#include "AliTaskCDBconnect.h"

class AliAnalysisDataContainer;
class AliClusterContainer;
class AliParticleContainer;

#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliEmcalCorrectionTask.h"
#include "AliAnalysisTaskCaloTrackCorrelation.h"
#include "AliAnaCaloTrackCorrMaker.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
// Tell ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)

#include "OADB/macros/AddTaskPhysicsSelection.C"
#include "OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"
#include "OADB/macros/AddTaskCentrality.C"
#include "PWGPP/PilotTrain/AddTaskCDBconnect.C"
//#include "PWG/EMCAL/macros/CreateAODChain.C"
//#include "PWG/EMCAL/macros/CreateESDChain.C"

#include "PWGGA/CaloTrackCorrelations/macros/AddTaskGammaHadronCorrelationSelectAnalysis.C" 

#endif

void LoadMacros();
void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode);
AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC);

//______________________________________________________________________________
AliAnalysisManager* runCaloTrackCorrEmbeddingAnalysis(
    const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
    const char   *cRunPeriod     = "LHC18q",                                // set the run period
    const char   *cLocalFiles    = "aodFiles.txt",                          // set the local list file
    const UInt_t  iNumEvents     = 1000,                                    // number of events to be analyzed
    const UInt_t  kPhysSel       = AliVEvent::kAnyINT |
                    AliVEvent::kCentral | AliVEvent::kSemiCentral, // physics selection
    const char   *cTaskName      = "EMCalEmbeddingAnalysis",                     // sets name of analysis manager
    // 0 = only prepare the analysis manager but do not start the analysis
    // 1 = prepare the analysis manager and start the analysis
    // 2 = launch a grid analysis
    Int_t         iStartAnalysis = 1,
    const UInt_t  iNumFiles      = 5,                                     // number of files analyzed locally
    const char   *cGridMode      = "test"
)
{
  // Setup period
  TString sRunPeriod(cRunPeriod);
  sRunPeriod.ToLower();

  // Set Run 2
  Bool_t bIsRun2 = kFALSE;
  if (sRunPeriod.Length() == 6 && (sRunPeriod.BeginsWith("lhc15") || sRunPeriod.BeginsWith("lhc16"))) bIsRun2 = kTRUE;

  // Set beam type
  AliAnalysisTaskEmcal::BeamType iBeamType = AliAnalysisTaskEmcal::kpp;
  if (sRunPeriod == "lhc10h" || sRunPeriod == "lhc11h" || sRunPeriod == "lhc15o") {
    iBeamType = AliAnalysisTaskEmcal::kAA;
  }
  else if (sRunPeriod == "lhc12g" || sRunPeriod == "lhc13b" || sRunPeriod == "lhc13c" ||
      sRunPeriod == "lhc13d" || sRunPeriod == "lhc13e" || sRunPeriod == "lhc13f" ||
      sRunPeriod == "lhc16q" || sRunPeriod == "lhc16r" || sRunPeriod == "lhc16s" ||
      sRunPeriod == "lhc16t") {
    iBeamType = AliAnalysisTaskEmcal::kpA;
  }

  // Setup track container
  AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);
  std::cout << "Default track cut period set to: " << AliTrackContainer::GetDefTrackCutsPeriod().Data() << "\n";

  // Set data file type
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
  if (iStartAnalysis == 1) {
    if (sLocalFiles == "") {
      Printf("You need to provide the list of local files!");
      return 0;
    }
    Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);
  }

  // Load macros needed for the analysis
  #ifndef __CLING__
  LoadMacros();
  #endif

  ////////////////////////
  /// Configuration options
  ////////////////////////

  // Embedding files list
  const std::string embeddedFilesList = "aodFilesEmbed.txt";

  // Debug options
  //AliLog::SetClassDebugLevel("AliAnalysisTaskEmcalEmbeddingHelper", AliLog::kDebug+0);

  // Determine track, cluster, and cell names
  const bool IsEsd = (iDataType == kEsd);
  // These names correspond to the _uncombined_ input objects that we are interestd in the external (embedded) event
  TString externalEventParticlesName = "mcparticles";
  // Empty because there are no particle level clusters
  TString externalEventClustersName = "";

  // General input object names
  TString tracksName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kTrack, IsEsd);
  TString emcalCellsName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCaloCells, IsEsd);
  TString emcalCellsNameCombined = emcalCellsName + "Combined";
  TString clustersName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCluster, IsEsd);
  // Combined (PbPb + embedded det level) clusters
  TString clustersNameCombined = clustersName + "Combined";

 
  emcalCellsName = "";
  emcalCellsNameCombined = "";
  clustersName = "";
  clustersNameCombined = "";
  
  ///////////////////////////////
  /// Setup and Configure Analysis
  ///////////////////////////////

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

  // Create Input Handler
  if (iDataType == kAod) {
    AliAODInputHandler * pESDHandler = AliAnalysisTaskEmcal::AddAODHandler();
  }
  else {  
    AliESDInputHandler * pESDHandler = AliAnalysisTaskEmcal::AddESDHandler();
  }

  // Physics selection task
  if (iDataType == kEsd) {
    AliPhysicsSelectionTask * pPhysSelTask = AddTaskPhysicsSelection();
  }

  // Centrality task
  // The Run 2 condition is too restrictive, but until the switch to MultSelection is complete, it is the best we can do
  if (iDataType == kEsd && iBeamType != AliAnalysisTaskEmcal::kpp && bIsRun2 == kFALSE) {
    AliCentralitySelectionTask * pCentralityTask = AddTaskCentrality(kTRUE);
    pCentralityTask->SelectCollisionCandidates(AliVEvent::kAny);
  }
  // AliMultSelection
  // Works for pp, pPb, and PbPb for the periods that it is calibrated
  if (bIsRun2 == kTRUE) {
    AliMultSelectionTask * pMultSelectionTask = AddTaskMultSelection(kFALSE);
    pMultSelectionTask->SelectCollisionCandidates(AliVEvent::kAny);
  }

  // CDBconnect task
  AliTaskCDBconnect * taskCDB = AddTaskCDBconnect();
  taskCDB->SetFallBackToRaw(kTRUE);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Setup embedding task
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
  embeddingHelper->SelectCollisionCandidates(kPhysSel);
  // The pt hard bin should be set via the filenames in this file
  // If using a file pattern, it could be configured via embeddingHelper->SetPtHardBin(ptHardBin);
  embeddingHelper->SetFileListFilename(embeddedFilesList.c_str());

  // Some example settings for LHC12a15e_fix (anchored to LHC11h)
  embeddingHelper->SetNPtHardBins(11);
  embeddingHelper->SetMCRejectOutliers();

  // Setup internal event selection and additional configuration options
  embeddingHelper->SetConfigurationPath("./EmbeddingConfigurationCaloTrackCorrExample.yaml");

  // Initialize the task to complete the setup.
  embeddingHelper->Initialize();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// EMCal corrections
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TObjArray correctionTasks;
  
  // Create the Correction Tasks
  // "data" corresponds to the PbPb level
  // "embed" corresponds to the embedded detector level
  // "combined" corresponds to the hybrid (PbPb + embedded detector) level
  correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("data"));
  correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("embed"));
  // It is important that combined is last!
  correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("combined"));
  
  // Loop over all of the correction tasks to configure them
  AliEmcalCorrectionTask * tempCorrectionTask = 0;
  TIter next(&correctionTasks);
  while (( tempCorrectionTask = static_cast<AliEmcalCorrectionTask *>(next())))
  {
    tempCorrectionTask->SelectCollisionCandidates(kPhysSel);
    // Configure centrality
    tempCorrectionTask->SetNCentBins(5);
    if (bIsRun2) {
      tempCorrectionTask->SetUseNewCentralityEstimation(kTRUE);
    }
    tempCorrectionTask->SetUserConfigurationFilename("./EMCalCorrConfig_EmbeddMCData_ClV1_xTalk.yaml");
    
    tempCorrectionTask->Initialize();
  }

  TObjArray *pTopTasks = pMgr->GetTasks();
   for (Int_t i = 0; i < pTopTasks->GetEntries(); ++i) {
     AliAnalysisTaskSE *pTask = dynamic_cast<AliAnalysisTaskSE*>(pTopTasks->At(i));
     if (!pTask) continue;

     if (pTask->InheritsFrom("AliEmcalCorrectionTask")) {
       AliEmcalCorrectionTask * pTaskEmcalCorrection = static_cast<AliEmcalCorrectionTask*>(pTask);
       Printf("Setting beam type %d for task %s", iBeamType, pTaskEmcalCorrection->GetName());
       pTaskEmcalCorrection->SetForceBeamType(static_cast<AliEmcalCorrectionTask::BeamType>(iBeamType));
     }
   }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// //////////////////// CaloTrackCorr ////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  #if defined(__CINT__)
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/AddTaskGammaHadronCorrelationSelectAnalysis.C");
  #endif
  
  
  // Common settings for Correlation and QA tasks
  Bool_t   calibrate     = kFALSE;
  Int_t    minCen        = -1;
  Int_t    maxCen        = -1;
  Int_t    debug         = -1;
  
  Int_t    rejectEMCTrig = 0;
  Bool_t   nonLinOn      = kFALSE;
  Float_t  shshMax       = 0.3;
  Float_t  isoCone       = 0.4;
  Float_t  isoConeMin    = -1;
  Float_t  isoPtTh       = 2;
  Int_t    isoMethod     = AliIsolationCut::kSumBkgSubIC;//kSumPtIC;//kSumBkgSubIC;//kSumPtIC; //kSumBkgSubPhiBandIC
  Int_t    isoContent    = AliIsolationCut::kNeutralAndCharged;//kOnlyCharged;//kNeutralAndCharged;//
  Int_t    leading       = 0;
  Int_t    tm            = 2;
  Bool_t   mixOn         = kFALSE;
  TString  outputfile    = "";
  Bool_t   printSettings = kFALSE;
  
  Bool_t  kMC = kFALSE;
  TString kCollision= "PbPb";
  Int_t   kYear = 2018;
  TString kPeriod = "LHC18q";
  
  TString  cutSelected      = "SPDPileUp";  
  TString  analysisSelected = "QA_Charged";
  
  TString trigger       = "default";// MB
  TString clustersArray = "";//"caloClustersCombined"; // "caloClusters"
  TString cellsArray    = "";//"emcalCellsCombined"; // "emcalCells"
  
  // EMCal analysis
  //
  // *** Data ***
  AliAnalysisTaskCaloTrackCorrelation * emcData = AddTaskGammaHadronCorrelationSelectAnalysis
  ("EMCAL",kMC,kYear,kCollision,kPeriod,rejectEMCTrig,"",cutSelected,calibrate,nonLinOn, analysisSelected,
   shshMax,isoCone,isoConeMin,isoPtTh,isoMethod ,isoContent,leading,
   tm,minCen,maxCen,mixOn,outputfile,printSettings,debug,trigger);
  
  emcData->GetAnalysisMaker()->GetReader()->SetEMCALClusterListName(clustersArray);
  emcData->GetAnalysisMaker()->GetReader()->SetEMCALCellsListName(cellsArray);
  //emcData->GetAnalysisMaker()->GetReader()->SwitchOffRejectNoTrackEvents();
  
  emcData->GetAnalysisMaker()->GetReader()->SwitchOffUseEMCALTimeCut();
  emcData->GetAnalysisMaker()->GetReader()->SetEMCALTimeCut(-1e10,1e10); // Open time cut
  
  //printf("***EMCal FINAL Reader SETTINGS\n");
  
  emcData->GetAnalysisMaker()->GetReader()->Print("");
  
  // TList * anaListEmc = emcData->GetAnalysisMaker()->GetListOfAnalysisContainers();
  
  // *** Combined ***
  kMC = kTRUE;
  clustersArray = "caloClustersCombined";//"caloClustersCombined"; // "caloClusters"
  cellsArray    = "emcalCellsCombined";//"emcalCellsCombined"; // "emcalCells"
  AliAnalysisTaskCaloTrackCorrelation * emcComb = AddTaskGammaHadronCorrelationSelectAnalysis
  ("EMCAL",kMC,kYear,kCollision,kPeriod,rejectEMCTrig,"",cutSelected+"_EmbedMC",
   calibrate,nonLinOn, analysisSelected,
   shshMax,isoCone,isoConeMin,isoPtTh,isoMethod ,isoContent,leading,
   tm,minCen,maxCen,mixOn,outputfile,printSettings,debug,trigger);
    
  emcComb->GetAnalysisMaker()->GetReader()->SetEMCALClusterListName(clustersArray);
  emcComb->GetAnalysisMaker()->GetReader()->SetEMCALCellsListName(cellsArray);
  
  //emcComb->GetAnalysisMaker()->GetReader()->SwitchOffRejectNoTrackEvents();
  
  emcComb->GetAnalysisMaker()->GetReader()->SwitchOffUseEMCALTimeCut();
  emcComb->GetAnalysisMaker()->GetReader()->SetEMCALTimeCut(-1e10,1e10); // Open time cut
  
  //printf("***EMCal FINAL Reader SETTINGS\n");
  
  emcComb->GetAnalysisMaker()->GetReader()->Print("");
  
  // TList * anaListEmc = emcComb->GetAnalysisMaker()->GetListOfAnalysisContainers();
  
  // *** External MC ***
  clustersArray = "";//"caloClustersCombined"; // "caloClusters"
  cellsArray    = "";//"emcalCellsCombined"; // "emcalCells"
  AliAnalysisTaskCaloTrackCorrelation * emcExt = AddTaskGammaHadronCorrelationSelectAnalysis
  ("EMCAL",kMC,kYear,kCollision,kPeriod,rejectEMCTrig,"",cutSelected+"_EmbedMCInput",
   calibrate,nonLinOn, analysisSelected,
   shshMax,isoCone,isoConeMin,isoPtTh,isoMethod ,isoContent,leading,
   tm,minCen,maxCen,mixOn,outputfile,printSettings,debug,trigger);
  
  //emcExt->UseEmbeddedEvent(kTRUE,kTRUE); 
  
  emcExt->GetAnalysisMaker()->GetReader()->SetEMCALClusterListName(clustersArray);
  emcExt->GetAnalysisMaker()->GetReader()->SetEMCALCellsListName(cellsArray);
  //emcComb->GetAnalysisMaker()->GetReader()->SwitchOffRejectNoTrackEvents();
  
  emcExt->GetAnalysisMaker()->GetReader()->SwitchOffUseEMCALTimeCut();
  emcExt->GetAnalysisMaker()->GetReader()->SetEMCALTimeCut(-1e10,1e10); // Open time cut
  
  //printf("***EMCal FINAL Reader SETTINGS\n");
  
  emcExt->GetAnalysisMaker()->GetReader()->Print("");
  
  // TList * anaListEmc = emcComb->GetAnalysisMaker()->GetListOfAnalysisContainers();
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
 
  if (!pMgr->InitAnalysis()) return 0;
  pMgr->PrintStatus();
    
  //pMgr->SetDebugLevel(10);
  pMgr->SetUseProgressBar(kTRUE, 100);
  
  // Commented, it crashes, why??
//  TFile *pOutFile = new TFile("train.root","RECREATE");
//  pOutFile->cd();
//  pMgr->Write();
//  pOutFile->Close();
//  delete pOutFile;
  
  if (iStartAnalysis == 1) { // start local analysis
    TChain* pChain = 0;
    if (iDataType == kAod) {
      printf("AOD\n");
      #ifdef __CLING__
      std::stringstream aodChain;
      aodChain << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/CreateAODChain.C(";
      aodChain << "\"" << sLocalFiles.Data() << "\", ";
      aodChain << iNumEvents << ", ";
      aodChain << 0 << ", ";
      aodChain << std::boolalpha << kFALSE << ");";
      pChain = reinterpret_cast<TChain *>(gROOT->ProcessLine(aodChain.str().c_str()));
      #else
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
      #endif
      printf("AOD END\n");

    }
    else {
      #ifdef __CLING__
      std::stringstream esdChain;
      esdChain << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<  "/PWG/EMCAL/macros/CreateESDChain.C(";
      esdChain << "\"" << sLocalFiles.Data() << "\", ";
      esdChain << iNumEvents << ", ";
      esdChain << 0 << ", ";
      esdChain << std::boolalpha << kFALSE << ");";
      pChain = reinterpret_cast<TChain *>(gROOT->ProcessLine(esdChain.str().c_str()));
      #else
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      pChain = CreateESDChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
      #endif
    }

    // start analysis
    Printf("Starting Analysis...");
    pMgr->StartAnalysis("local", pChain, iNumEvents);
  }
  else if (iStartAnalysis == 2) {  // start grid analysis
    StartGridAnalysis(pMgr, cTaskName, cGridMode);
  }

  return pMgr;
}

void LoadMacros()
{
  // Aliroot macros
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
}

void StartGridAnalysis(AliAnalysisManager* pMgr, const char* uniqueName, const char* cGridMode)
{
  Int_t maxFilesPerWorker = 4;
  Int_t workerTTL = 7200;
  const char* runNumbers = "180720";
  const char* pattern = "pass2/AOD/*/AliAOD.root";
  const char* gridDir = "/alice/data/2012/LHC12c";
  const char* additionalCXXs = "";
  const char* additionalHs = "";

  AliAnalysisGrid *plugin = CreateAlienHandler(uniqueName, gridDir, cGridMode, runNumbers, pattern, additionalCXXs, additionalHs, maxFilesPerWorker, workerTTL, kFALSE);
  pMgr->SetGridHandler(plugin);

  // start analysis
  Printf("Starting GRID Analysis...");
  pMgr->SetDebugLevel(0);
  pMgr->StartAnalysis("grid");
}

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir, const char* gridMode, const char* runNumbers,
    const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker, Int_t workerTTL, Bool_t isMC)
{
  TDatime currentTime;
  TString tmpName(uniqueName);

  // Only add current date and time when not in terminate mode! In this case the exact name has to be supplied by the user
  if (strcmp(gridMode, "terminate")) {
    tmpName += "_";
    tmpName += currentTime.GetDate();
    tmpName += "_";
    tmpName += currentTime.GetTime();
  }

  TString macroName("");
  TString execName("");
  TString jdlName("");
  macroName = Form("%s.C", tmpName.Data());
  execName = Form("%s.sh", tmpName.Data());
  jdlName = Form("%s.jdl", tmpName.Data());

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetRunMode(gridMode);

  // Here you can set the (Ali)PHYSICS version you want to use
  plugin->SetAliPhysicsVersion("vAN-20160203-1");

  plugin->SetGridDataDir(gridDir); // e.g. "/alice/sim/LHC10a6"
  plugin->SetDataPattern(pattern); //dir structure in run directory

  if (!isMC) plugin->SetRunPrefix("000");

  plugin->AddRunList(runNumbers);

  plugin->SetGridWorkingDir(Form("work/%s",tmpName.Data()));
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

  plugin->SetAnalysisSource(additionalCode.Data());

  plugin->SetDefaultOutputs(kTRUE);
  plugin->SetAnalysisMacro(macroName.Data());
  plugin->SetSplitMaxInputFileNumber(maxFilesPerWorker);
  plugin->SetExecutable(execName.Data());
  plugin->SetTTL(workerTTL);
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName(jdlName.Data());
  plugin->SetPrice(1);
  plugin->SetSplitMode("se");

  // merging via jdl
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  plugin->SetMaxMergeStages(2);

  return plugin;
}
