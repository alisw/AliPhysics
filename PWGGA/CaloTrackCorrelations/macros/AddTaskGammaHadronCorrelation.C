/// \file AddTaskGammaHadronCorrelation.C
/// \ingroup CaloTrackCorrMacros
/// \brief Configuration of gamma-hadron and pi0-hadron + isolation, correlation analysis
///
/// Configuration macro for analysis of gamma-hadron and pi0-hadron correlation analysis
/// both pi0 and gamma isolated or not.
/// It does:
///   * Photon selection in calorimeter with AliAnaPhoton: Track matching, mild shower shape, NLM ... cuts
///   * Tagging of photon candidate clusters as decay from pi0 or eta or from their side bands with AliAnaPi0EbE (4 calls, one pi0, one eta, one pi0 side band and one eta side band)
///   * Tagging of photon candidate as isolated with AliAnaParticleIsolation
///   * Correlation of the photon and charged tracks, twice, one without isolation condition, and other with isolation condition
///   * Pi0 selection with the identification of merged clusters with AliAnaPi0EbE
///   * Tagging of pi0 as isolated with AliAnaParticleIsolation
///   * Correlation of the pi0 and charged tracks, twice, one without isolation condition, and other with isolation condition
///   * Optionally, the QA tasks AliAnaCalorimeterQA and AliAnaChargedParticle are executed
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

//// Set includes for compilation
//
//#if !defined(__CINT__) || defined(__MAKECINT__)
//
//#include <TString.h>
//#include <TROOT.h>
//
//#include "AliLog.h"
//#include "AliAnalysisTaskCaloTrackCorrelation.h"
//#include "AliCaloTrackESDReader.h"
//#include "AliCaloTrackAODReader.h"
//#include "AliCalorimeterUtils.h"
//#include "AliAnaPhoton.h"
//#include "AliAnaPi0EbE.h"
//#include "AliHistogramRanges.h"
//#include "AliAnaParticleIsolation.h"
//#include "AliAnaParticleHadronCorrelation.h"
//#include "AliAnaChargedParticles.h"
//#include "AliAnaCalorimeterQA.h"
//#include "AliAnaGeneratorKine.h"
//#include "AliAnalysisTaskCaloTrackCorrelation.h"
//#include "AliAnaCaloTrackCorrMaker.h"
//#include "AliAnalysisManager.h"
//#include "AliInputEventHandler.h"
//#include "AliVTrack.h"
//#include "ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C"
//#include "AliESDtrackCuts.h"
//#include "CreateTrackCutsPWGJE.C"
//#include "ConfigureEMCALRecoUtils.C"
//#endif
//
//// Declare methods for compilation
//
//AliCaloTrackReader  * ConfigureReader        (TString col,           Bool_t simulation,
//                                              TString clustersArray, Bool_t tender,
//                                              TString calorimeter,   Bool_t nonLinOn,
//                                              TString trigger,       Bool_t rejectEMCTrig,
//                                              Int_t   minCen,        Int_t  maxCen,
//                                              Bool_t  printSettings, Int_t  debug       );
//AliCalorimeterUtils * ConfigureCaloUtils     (TString col,           Bool_t simulation,
//                                                                     Bool_t tender,
//                                              Bool_t  nonLinOn,      Int_t  year,
//                                              Bool_t  printSettings, Int_t  debug            );
//AliAnaPhoton        * ConfigurePhotonAnalysis(TString col,           Bool_t simulation,
//                                              TString calorimeter,   Int_t  year,  Int_t tm,
//                                              Bool_t  printSettings, Int_t  debug            );
//AliAnaPi0EbE        * ConfigurePi0EbEAnalysis(TString particle,      Int_t  analysis,
//                                              Bool_t useSSIso,       Bool_t useAsy,
//                                              TString col,           Bool_t simulation,
//                                              TString calorimeter,   Int_t  year, Int_t tm,
//                                              Bool_t  printSettings, Int_t  debug            );
//AliAnaParticleIsolation* ConfigureIsolationAnalysis
//                                             (TString particle,      Int_t   leading,
//                                              Int_t partInCone,      Int_t   thresType,
//                                              Float_t cone,          Float_t pth,      Bool_t multi,
//                                              TString col,           Bool_t  simulation,
//                                              TString calorimeter,   Int_t   year,      Int_t tm,
//                                              Bool_t  printSettings, Int_t   debug                    );
//AliAnaParticleHadronCorrelation * ConfigureHadronCorrelationAnalysis
//                                             (TString particle,      Int_t   leading,
//                                              Bool_t  bIsolated,     Float_t shshMax,
//                                              Int_t partInCone,      Int_t   thresType,
//                                              Float_t cone,          Float_t pth,      Bool_t mixOn,
//                                              TString col,           Bool_t  simulation,
//                                              TString calorimeter,   Int_t   year,     Int_t tm,
//                                              Bool_t  printSettings, Int_t   debug                     );
//AliAnaChargedParticles* ConfigureChargedAnalysis
//                                             (Bool_t simulation,     Bool_t  printSettings, Int_t   debug);
//AliAnaCalorimeterQA * ConfigureQAAnalysis    (TString col,           Bool_t  simulation,
//                                              TString calorimeter,   Int_t   year,
//                                              Bool_t  printSettings, Int_t   debug                     );
//AliAnaGeneratorKine* ConfigureGenKineAnalysis(Int_t   thresType,     Float_t cone,      Float_t pth,
//                                              TString col,           Bool_t  simulation,
//                                              TString calorimeter,   Int_t   year,
//                                              Bool_t  printSettings, Int_t   debug                     );
//
//void SetAnalysisCommonParameters             (AliAnaCaloTrackCorrBaseClass* ana,
//                                              TString calorimeter , Int_t  year,
//                                              TString col         , Bool_t simulation,
//                                              Bool_t printSettings, Int_t  debug);


/// Global name to be composed of the settings, used to set the AOD branch name
TString kAnaGammaHadronCorr = "";

///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle
/// \param simulation : A bool identifying the data as simulation
/// \param year: The year the data was taken, used to configure some histograms
/// \param col: A string with the colliding system
/// \param rejectEMCTrig : An int to reject EMCal triggered events with bad trigger: 0 no rejection, 1 old runs L1 bit, 2 newer runs L1 bit
/// \param clustersArray : A string with the array of clusters not being the default (default is empty string)
/// \param tender : A bool indicating if the tender was running before this analysis
/// \param nonLinOn : A bool to set the use of the non linearity correction
/// \param shshMax : A float setting the maximum value of the shower shape of the clusters for the correlation analysis
/// \param isoCone : A float setting the isolation cone size
/// \param isoPtTh : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param isoMethod : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
/// \param isoContent : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
/// \param leading : An int setting the type of leading particle selection: 0, select all;l 1: absolute  leading of charged; 2: absolute  leading of charged and neutral; 3: near side leading absolute of charged; 4: near side leading absolute of charged and neutral
/// \param tm : track matching options: 0- no matching; 1-fixed residual cuts; 2-pT track dependent cut
/// \param minCen : An int to select the minimum centrality, -1 means no selection
/// \param maxCen : An int to select the maximum centrality, -1 means no selection
/// \param mixOn : A bool to switch the correlation mixing analysis
/// \param qaAn : A bool to switch the calorimeter QA analysis
/// \param chargedAn : A bool to switch the selected tracks QA analysis
/// \param outputfile : A string to change the name of the histograms output file, default is AnalysisResults.root
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param trigSuffix :  A string with the trigger class, abbreviated, defined in ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C
///
AliAnalysisTaskCaloTrackCorrelation * AddTaskGammaHadronCorrelation
(
 TString  calorimeter   = "EMCAL",
 Bool_t   simulation    = kFALSE,
 Int_t    year          = 2011,
 TString  col           = "pp",
 Int_t    rejectEMCTrig = 0,
 TString  clustersArray = "",
 Bool_t   tender        = kFALSE,
 Bool_t   nonLinOn      = kFALSE,
 Float_t  shshMax       = 0.27,
 Float_t  isoCone       = 0.4,
 Float_t  isoPtTh       = 0.5,
 Int_t    isoMethod     = AliIsolationCut::kPtThresIC,
 Int_t    isoContent    = AliIsolationCut::kNeutralAndCharged,
 Int_t    leading       = 0,
 Int_t    tm            = 1,
 Int_t    minCen        = -1,
 Int_t    maxCen        = -1,
 Bool_t   mixOn         = kTRUE,
 Bool_t   qaAn          = kFALSE,
 Bool_t   chargedAn     = kFALSE,
 TString  outputfile    = "",
 Bool_t   printSettings = kFALSE,
 Int_t    debug         = 0,  
 const char *trigSuffix = "EMC7"
)
{
  TString trigger = trigSuffix;
  
  printf("Passed settings:\n calorimeter <%s>, simulation <%d>, year <%d>,\n col <%s>, trigger <%s>, reject EMC <%d>, clustersArray <%s>, tender <%d>, non linearity <%d>\n shshMax <%2.2f>, isoCone <%2.2f>, isoPtTh <%2.2f>, isoMethod <%d>,isoContent <%d>,\n leading <%d>, tm <%d>, minCen <%d>, maxCen <%d>, mixOn <%d>,\n qaAn <%d>, chargedAn <%d>, outputfile <%s>, printSettings <%d>, debug <%d>\n",
         calorimeter.Data(),simulation,year,col.Data(),trigger.Data(), rejectEMCTrig, clustersArray.Data(),tender, nonLinOn, shshMax,
         isoCone,isoPtTh,isoMethod,isoContent,leading,tm,
         minCen,maxCen,mixOn,qaAn,chargedAn,outputfile.Data(),printSettings,debug);
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }
  
  // Make sure the B field is enabled for track selection, some cuts need it
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);
  
  // Name for containers
  
  //kAnaGammaHadronCorr = Form("%s_Trig%s_Cl%s_TM%d_R%1.1f_Pt%1.1f",calorimeter.Data(), trigger.Data(),kClusterArray.Data(),tm,cone,pth);
  kAnaGammaHadronCorr = Form("GammaHadron_%s_Trig%s_Col_%s_Year%d_Cl%s_Ten%d_TM%d_M02_%1.2f_IsoParam_C%1.2fPt%1.2fM%dPa%d_Lead%d_Mix%d",
                             calorimeter.Data(),trigger.Data(),col.Data(),year,clustersArray.Data(),tender,
                             tm, shshMax,isoCone,isoPtTh,isoMethod,isoContent,leading,mixOn);
  
  if(col=="PbPb" && maxCen>=0) kAnaGammaHadronCorr+=Form("Cen%d_%d",minCen,maxCen);
  
  printf("<<<< NAME: %s >>>>>\n",kAnaGammaHadronCorr.Data());
  
  // #### Configure analysis ####
  
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  
  // General frame setting and configuration
  maker->SetReader   ( ConfigureReader   (col,simulation,clustersArray,tender,calorimeter,nonLinOn,trigger,rejectEMCTrig,minCen,maxCen,printSettings,debug) );
  maker->SetCaloUtils( ConfigureCaloUtils(col,simulation,tender,nonLinOn,year,printSettings,debug) );
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important
  
  //
  // Photon analysis
  //
  maker->AddAnalysis(ConfigurePhotonAnalysis(col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Photon cluster selection
  
  
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0"        , AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Pi0 event by event selection, invariant mass and photon tagging from decay
  
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Eta"        , AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Eta event by event selection, invariant mass and photon tagging from decay
  
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0SideBand", AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Pi0 out of peak event by event selection, and photon tagging from decay
  
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("EtaSideBand", AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Eta out of peak event by event selection, and photon tagging from decay
  
  
  maker->AddAnalysis(ConfigureIsolationAnalysis("Photon", leading, isoContent,isoMethod,isoCone,isoPtTh, kFALSE,
                                                col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Photon isolation
  
  
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon", leading, kFALSE, shshMax, isoContent,isoMethod,isoCone,isoPtTh, mixOn,
                                                        col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Gamma-hadron correlation
  
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon", leading, kTRUE,  shshMax, isoContent,isoMethod,isoCone,isoPtTh, mixOn,
                                                        col,simulation,calorimeter,year,tm,printSettings,debug) , n++); // Isolated gamma hadron correlation
  
  //
  // Merged pi0 analysis
  //
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kSSCalo,kTRUE,kTRUE,
                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Pi0 event by event selection, cluster splitting
  
  maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SS", leading, isoContent,isoMethod,isoCone,isoPtTh, kFALSE,
                                                col,simulation,calorimeter,year,tm,printSettings,debug), n++);          // Pi0 isolation, cluster splits
  
  
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS", leading, kFALSE, shshMax, isoContent,isoMethod,isoCone,isoPtTh, mixOn,
                                                        col,simulation,calorimeter,year,tm,printSettings,debug), n++);  // Pi0-hadron correlation
  
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS", leading, kTRUE,  shshMax, isoContent,isoMethod,isoCone,isoPtTh, mixOn,
                                                        col,simulation,calorimeter,year,tm,printSettings,debug) , n++); // Isolated pi0-hadron correlation
  
  // Check the generated kinematics
  if(simulation)  maker->AddAnalysis(ConfigureGenKineAnalysis(isoMethod,isoCone,isoPtTh,
                                                              col,simulation,calorimeter,year,printSettings,debug), n++);
  
  // Charged analysis
  if(chargedAn)   maker->AddAnalysis(ConfigureChargedAnalysis(simulation,printSettings,debug), n++); // track selection checks
  
  // Calo QA
  if(qaAn)        maker->AddAnalysis(ConfigureQAAnalysis(col,simulation,calorimeter,year,printSettings,debug) , n++);
  
  maker->SetAnaDebug(debug)  ;
  
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  
  if( simulation || !trigger.Contains("EMC") ) maker->SwitchOffDataControlHistograms();
  
  if(simulation)
  {
    // Calculate the cross section weights, apply them to all histograms 
    // and fill xsec and trial histo. Sumw2 must be activated.
    //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionCalculation(); 
    //maker->SwitchOnSumw2Histograms();
    
    // For recent productions where the cross sections and trials are not stored in separate file
    //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionFromEventHeader() ;
    
    // Just fill cross section and trials histograms.
    maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionHistoFill(); 
    
    // Add control histogram with pT hard to control aplication of weights 
    maker->SwitchOnPtHardHistogram();
  }
  
  if(printSettings) maker->Print("");
  
  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, calorimeter.Data());
  
  //
  // Create task, pass the maker and add it to the manager
  //
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("%s",kAnaGammaHadronCorr.Data()));
  
  task->SetDebugLevel(debug);
  
  //task->SetBranches("ESD:AliESDRun.,AliESDHeader");
  //task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  
  task->SetAnalysisMaker(maker);
  
  mgr->AddTask(task);

  //
  // Select events trigger depending on trigger
  //
  maker->GetReader()->SwitchOnEventTriggerAtSE(); // on is default case
  if(!simulation)
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C");
    TString caloTriggerString = "";
    UInt_t mask = ConfigureAndGetEventTriggerMaskAndCaloTriggerString(trigger, year, caloTriggerString);

    maker->GetReader()->SetFiredTriggerClassName(caloTriggerString);

    // For mixing with AliAnaParticleHadronCorrelation switch it off
    if(mixOn)
    {
      maker->GetReader()->SwitchOffEventTriggerAtSE();
      maker->GetReader()->SetEventTriggerMask(mask); 
      // what to do with caloTriggerString?
      
      // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
      //reader->SetMixEventTriggerMask(AliVEvent::kMB); 
      maker->GetReader()->SetMixEventTriggerMask(AliVEvent::kINT7); 
      
      printf("---Trigger selection done in AliCaloTrackReader!!!\n");
    }
    else
    {
      task ->SelectCollisionCandidates( mask );
    }
  }
  
  //
  // Create containers
  //
  if(outputfile.Length()==0) outputfile = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(kAnaGammaHadronCorr, TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             Form("%s",outputfile.Data()));
  
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",kAnaGammaHadronCorr.Data()), TList::Class(),
                                                             AliAnalysisManager::kParamContainer,
                                                             "AnalysisParameters.root");
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  //if(!kData.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  
  return task;
}

///
/// Configure the class handling the events and cluster/tracks filtering.
///
AliCaloTrackReader * ConfigureReader(TString col,           Bool_t simulation,
                                     TString clustersArray, Bool_t tender,
                                     TString calorimeter,   Bool_t nonLinOn,
                                     TString trigger,       Bool_t rejectEMCTrig,
                                     Int_t   minCen,        Int_t  maxCen,
                                     Bool_t  printSettings, Int_t  debug         )
{
  // Get the data type ESD or AOD
  AliAnalysisManager * mgr = AliAnalysisManager::GetAnalysisManager();
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();
  
  AliCaloTrackReader * reader = 0;
  if     (inputDataType == "AOD") reader = new AliCaloTrackAODReader();
  else if(inputDataType == "ESD") reader = new AliCaloTrackESDReader();
  else printf("AliCaloTrackReader::ConfigureReader() - Data not known InputData=%s\n",inputDataType.Data());
  
  reader->SetDebug(debug);//10 for lots of messages
  
  //
  // MC settings
  //
  // Check if kine stack is available, independent of request of simulation
//  Bool_t useKinematics = kFALSE;
//  useKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
//  
//  if(simulation)
//  {
//    if (!useKinematics && inputDataType=="AOD") useKinematics = kTRUE; //AOD primary should be available ...
//  }
    
  // In case of Pythia pt Hard bin simulations (jet-jet, gamma-jet)
  // reject some special events that bother the cross section
  if(simulation)
  {
    // Event rejection cuts for jet-jet simulations, do not use in other
    reader->SetPtHardAndJetPtComparison(kTRUE);
    reader->SetPtHardAndJetPtFactor(4);
    
    reader->SetPtHardAndClusterPtComparison(kTRUE);
    reader->SetPtHardAndClusterPtFactor(1.5);
  }
  
  //------------------------
  // Detector input filling
  //------------------------
  
  //Min cluster/track E
  reader->SetEMCALEMin(0.3);
  reader->SetEMCALEMax(1000);
  reader->SetPHOSEMin(0.3);
  reader->SetPHOSEMax(1000);
  reader->SetCTSPtMin(0.2);
  reader->SetCTSPtMax(1000);
  
  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();
  
  // Shower shape smearing
  // Set it in the train configuration page not here for the moment
  //  if(simulation)
  //  {
  //    reader->SwitchOffShowerShapeSmearing(); // Active only on MC, off by default
  //    reader->SetShowerShapeSmearWidth(0.005);
  //  }
  
  //
  // Tracks
  //
  reader->SwitchOnCTS();
  
  reader->SwitchOffUseTrackTimeCut();
  reader->SetTrackTimeCut(0,50);
  
  reader->SwitchOnFiducialCut();
  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;
  
  reader->SwitchOffUseTrackDCACut();
  //reader->SetTrackDCACut(0,0.0105);
  //reader->SetTrackDCACut(1,0.035);
  //reader->SetTrackDCACut(2,1.1);
  
  if(inputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
    
    //AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
    //reader->SetTrackCuts(esdTrackCuts);
    
    AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
    reader->SetTrackCuts(esdTrackCuts);
    AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
    reader->SetTrackComplementaryCuts(esdTrackCuts2);
    
    reader->SwitchOnConstrainTrackToVertex();
  }
  else if(inputDataType=="AOD")
  {
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SwitchOnAODTrackSharedClusterSelection();
    reader->SetTrackStatus(AliVTrack::kITSrefit);
    
    //reader->SwitchOnAODPrimaryTrackSelection(); // Used in preliminary results of QM from Nicolas and Xiangrong?
    //reader->SwitchOnTrackHitSPDSelection();     // Check that the track has at least a hit on the SPD, not much sense to use for hybrid or TPC only tracks
    //reader->SetTrackFilterMask(128);            // Filter bit, not mask, use if off hybrid, TPC only
  }
  
  //
  // Calorimeter
  //
  if(clustersArray == "" && !tender)
  {
    printf("**************** Standard EMCAL clusters branch analysis **************** \n");
    reader->SwitchOnClusterRecalculation();
    // Check in ConfigureCaloUtils that the recalibration and bad map are ON
  }
  else
  {
    printf("**************** Input for analysis is Clusterizer %s **************** \n", clustersArray.Data());
    reader->SetEMCALClusterListName(clustersArray);
    reader->SwitchOffClusterRecalculation();
  }
  
  // Time cuts
  reader->SwitchOffUseParametrizedTimeCut();
  if(simulation)
  {
    reader->SwitchOffUseEMCALTimeCut();
    reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
  }
  else
  {
    reader->SwitchOnUseEMCALTimeCut();
    reader->SetEMCALTimeCut(-25,20);
  }
  
  // CAREFUL
  if(nonLinOn) reader->SwitchOnClusterELinearityCorrection();
  else         reader->SwitchOffClusterELinearityCorrection();
  
  if(calorimeter == "EMCAL")
  {
    reader->SwitchOnEMCALCells();
    reader->SwitchOnEMCAL();
  }
  
  if(calorimeter == "PHOS")
  { // Should be on if QA is activated with correlation on
    reader->SwitchOnPHOSCells();
    reader->SwitchOnPHOS();
  }
  
  //-----------------
  // Event selection
  //-----------------
  
  //if(!simulation) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
  
  // Event triggered by EMCal selection settings
  reader->SwitchOffTriggerPatchMatching();
  reader->SwitchOffBadTriggerEventsRemoval();
  
  if( rejectEMCTrig > 0 && !simulation && (trigger.Contains("EMC") || trigger.Contains("L")))
  {
    printf("=== Remove bad triggers === \n");
    reader->SwitchOnTriggerPatchMatching();
    reader->SwitchOnBadTriggerEventsRemoval();
    
    //    reader->SetTriggerPatchTimeWindow(8,9); // default values
    //    if     (kRunNumber < 146861) reader->SetEventTriggerL0Threshold(3.);
    //    else if(kRunNumber < 154000) reader->SetEventTriggerL0Threshold(4.);
    //    else if(kRunNumber < 165000) reader->SetEventTriggerL0Threshold(5.5);
    //    //redefine for other periods, triggers
    //
    //    if(kRunNumber < 172000)
    //    {
    //      reader->SetEventTriggerL1Bit(4,5); // current LHC11 data
    //      printf("\t Old L1 Trigger data format!\n");
    //    }
    //    else
    //    {
    //      reader->SetEventTriggerL1Bit(6,8); // LHC12-13 data
    //      printf("\t Current L1 Trigger data format!\n");
    //    }
    
    if(clustersArray != "" || tender)
    {
      printf("Trigger cluster calibration OFF\n");
      reader->SwitchOffTriggerClusterTimeRecal() ;
    }
    
  }
  
  //reader->RejectFastClusterEvents() ;
    
  reader->SetZvertexCut(10.);               // Open cut
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
  reader->SwitchOnRejectNoTrackEvents();
  
  reader->SwitchOffV0ANDSelection() ;       // and besides v0 AND
  reader->SwitchOffPileUpEventRejection();  // remove pileup by default off, apply it only for MB not for trigger
  
  if(col=="PbPb")
  {
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
    reader->SetCentralityBin(minCen,maxCen); // Accept all events, if not select range
    //reader->SwitchOnAcceptOnlyHIJINGLabels();
    // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
    reader->SetEventPlaneMethod("V0");
  }
  
  if(printSettings) reader->Print("");
  
  return reader;
}

///
/// Configure the class handling the calorimeter clusters specific methods
///
AliCalorimeterUtils* ConfigureCaloUtils(TString col,    Bool_t simulation,
                                        Bool_t  tender, Bool_t nonLinOn,      
                                        Int_t   year,   Bool_t printSettings, 
                                        Int_t   debug)
{
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  
  cu->SetDebug(debug);
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder (2);
  
  cu->SetNumberOfSuperModulesUsed(10);
  
  if     (year == 2010) cu->SetNumberOfSuperModulesUsed(4);
  else if(year <= 2013) cu->SetNumberOfSuperModulesUsed(10);
  else if(year >  2013) cu->SetNumberOfSuperModulesUsed(20);
  else                  cu->SetNumberOfSuperModulesUsed(10);
  
  printf("xxx Number of SM set to <%d> xxx\n",cu->GetNumberOfSuperModulesUsed());
  
  // Search of local maxima in cluster
  if(col=="pp")
  {
    cu->SetLocalMaximaCutE(0.1);
    cu->SetLocalMaximaCutEDiff(0.03);
  }
  else
  {
    cu->SetLocalMaximaCutE(0.2);
    cu->SetLocalMaximaCutEDiff(0.03);
  }
  
  cu->SwitchOffRecalculateClusterTrackMatching();
  
  cu->SwitchOnBadChannelsRemoval() ;
  
  // EMCAL settings
  
  if(!simulation)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  
  // calibrations
  Bool_t calibEner = kFALSE;
  Bool_t calibTime = kFALSE;
  cu->SwitchOffRecalibration();
  cu->SwitchOffRunDepCorrection();
  
  if( !tender )
  {
    cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
    cu->SwitchOnRunDepCorrection();
    
    calibEner = kTRUE;
    calibTime = kTRUE;
  }
  
  if( simulation )
  {
    calibEner = kFALSE;
    calibTime = kFALSE;
    
    cu->SwitchOffRecalibration(); // Check the reader if it is taken into account during filtering
    cu->SwitchOffRunDepCorrection();
  }
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          simulation,
                          kTRUE,      // exotic
                          nonLinOn,   // Non linearity
                          calibEner,  // E calib
                          kTRUE,      // bad map
                          calibTime); // time calib
  
  if( calibTime ) recou->SetExoticCellDiffTimeCut(50);
  
  if( nonLinOn )  cu->SwitchOnCorrectClusterLinearity();
  
  printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
  printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
  
  // PHOS
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
  
  if(printSettings) cu->Print("");
  
  return cu;
}

///
/// Configure the task doing the first photon cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection ...
///
AliAnaPhoton* ConfigurePhotonAnalysis(TString col,           Bool_t simulation,
                                      TString calorimeter,   Int_t year, Int_t tm,
                                      Bool_t  printSettings, Int_t   debug)
{
  AliAnaPhoton *ana = new AliAnaPhoton();
  
  // cluster selection cuts
  
  ana->SwitchOnRealCaloAcceptance();
  
  ana->SwitchOffFiducialCut();
  
  ana->SetCalorimeter(calorimeter);
  
  ana->SetFirstSMCoveredByTRD(6);
  
  ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection
  
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  if(tm) ana->SwitchOnTrackMatchRejection() ;
  else   ana->SwitchOffTrackMatchRejection() ;
  
  ana->SwitchOnTMHistoFill() ;
  
  if(calorimeter == "PHOS")
  {
    ana->SetNCellCut(2);// At least 3 cells
    ana->SetMinPt(0.3);
    ana->SetMinDistanceToBadChannel(2, 4, 5);
    ana->SetTimeCut(-1e10,1e10); // open cut
  }
  else
  {//EMCAL
    ana->SetNCellCut(1);// At least 2 cells
    ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
    ana->SetMaxEnergy(100);
    ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off
    // restrict to less than 100 ns when time calibration is on
    ana->SetMinDistanceToBadChannel(2, 4, 6);
    
    // NLM cut, used in all, exclude clusters with more than 2 maxima
    // Not needed if M02 cut is already strong or clusterizer V2
    ana->SetNLMCut(1, 2) ;
  }
  
  //PID cuts (shower shape)
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
  //Not used in bayesian
  
  // EMCAL
  
  //caloPID->SetEMCALLambda0CutMax(0.27);
  caloPID->SetEMCALLambda0CutMax(10); // open, full shower shape needed for isolation studies
  caloPID->SetEMCALLambda0CutMin(0.10);
  
  // Track matching
  // tm = 1, fixed cuts
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
  
  // pT track dependent cuts
  if(tm > 1) caloPID->SwitchOnEMCTrackPtDepResMatching();
  
  // PHOS
  caloPID->SetPHOSDispersionCut(2.5);
  caloPID->SetPHOSRCut(2.);
  //if(kInputData=="AOD") caloPID->SetPHOSRCut(2000.); // Open cut since dX, dZ not stored
  
  // Branch AOD settings
  ana->SetOutputAODName(Form("PhotonTrigger_%s",kAnaGammaHadronCorr.Data()));
  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  
  //Set Histograms name tag, bins and ranges
  ana->AddToHistogramsName(Form("AnaPhoton_TM%d_",tm));
  
  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug) ; // see method below
  
  if(ana->GetFirstSMCoveredByTRD() > 0)
    printf(">>> Set first SM covered by TRD, SM=%d <<< year %d \n", ana->GetFirstSMCoveredByTRD(),year);
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms (17); // 18 max
  ana->FillNPrimaryHistograms(6); // 6 max
  
  return ana;
}

///
/// Configure the task doing the pi0 even by event selection (invariant mass or split)
/// and the cluster tagging as decay in different mass windows.
///
AliAnaPi0EbE* ConfigurePi0EbEAnalysis(TString particle,      Int_t  analysis,
                                      Bool_t useSSIso,       Bool_t useAsy,
                                      TString col,           Bool_t simulation,
                                      TString calorimeter,   Int_t  year,  Int_t tm,
                                      Bool_t  printSettings, Int_t debug           )
{
  // Configuration of pi0 event by event selection
  
  AliAnaPi0EbE *ana = new AliAnaPi0EbE();
  
  ana->SetDebug(debug);
  
  ana->SetAnalysisType((AliAnaPi0EbE::anaTypes)analysis);
  TString opt = "";
  if(analysis==AliAnaPi0EbE::kIMCaloTracks) opt = "Conv";
  if(analysis==AliAnaPi0EbE::kSSCalo)       opt = "SS";
  
  ana->SwitchOffAllNLMHistoFill();
  ana->SetFirstSMCoveredByTRD(6);
  ana->SwitchOffSelectedClusterHistoFill();
  
  ana->SwitchOffFillWeightHistograms();
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  //if(!kTime && !simulation) ana->SwitchOnFillEMCALBCHistograms();
  
  if(tm) ana->SwitchOnTrackMatchRejection() ;
  else   ana->SwitchOffTrackMatchRejection() ;
  ana->SwitchOffTMHistoFill() ;
  
  ana->SetCalorimeter(calorimeter);
  
  // Branch AOD settings
  ana->SetOutputAODName(Form("%s%sTrigger_%s",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
  printf("***Out branch %s***\n",ana->GetOutputAODName().Data());
  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  
  if(analysis == AliAnaPi0EbE::kIMCaloTracks) ana->SetInputAODGammaConvName("PhotonsCTS");
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("Ana%s%sEbE_TM%d_",particle.Data(),opt.Data(),tm));
  
  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  ///////////////////////////////////
  if(analysis!=AliAnaPi0EbE::kSSCalo)
  {
    ana->SetInputAODName(Form("PhotonTrigger_%s",kAnaGammaHadronCorr.Data()));
    
    ana->SetM02CutForInvMass(0.1,0.35); // Loose SS cut
    
    ana->SwitchOnSelectPairInIsolationCone();
    ana->SetR(0.4);
    ana->SetIsolationCandidateMinPt(5);
    
    if(useSSIso)
    {
      ana->SwitchOnSelectIsolatedDecay();
      ana->AddToHistogramsName(Form("Ana%s%sEbEIsoDecay_TM%d_",particle.Data(),opt.Data(),tm));
      ana->SetOutputAODName(Form("%s%sIsoDecayTrigger_%s",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
    }
    
    if(calorimeter=="EMCAL" && !simulation) ana->SetPairTimeCut(100);
    
    AliNeutralMesonSelection *nms = ana->GetNeutralMesonSelection();
    nms->SetParticle(particle);
    
    //****
    nms->SetInvMassCutMaxParameters(0,0,0); // Overrule the setting in SetParticle for Pi0 option
    //****
    
    // Tighten a bit mass cut with respect to default window
    if(particle=="Pi0") nms->SetInvMassCutRange(0.110,0.160);
    if(particle=="Eta") nms->SetInvMassCutRange(0.520,0.580);
    
    //if(!particle.Contains("SideBand")) nms->SwitchOnAngleSelection();
    //else nms->SwitchOnAngleSelection();
    
    nms->SwitchOffAngleSelection();
    
    if(particle.Contains("Pi0SideBand")) // For pi0, do not consider left band
      nms->SetSideBandCutRanges(-1,0,0.190,0.240);
    
    if(particle.Contains("EtaSideBand")) // For pi0, do not consider left band
      nms->SetSideBandCutRanges(0.410,0.470,0.620,0.680);
    
    nms->KeepNeutralMesonSelectionHistos(kTRUE);
    //nms->SetAngleMaxParam(2,0.2);
    nms->SetHistoERangeAndNBins(0, 20, 80) ;
    //nms->SetHistoIMRangeAndNBins(0, 1, 400);
  }
  else  ///////////////////////////////////
  {
    // cluster splitting settings
    ana->SetMinEnergy(6);
    ana->SetMaxEnergy(100.);
    
    ana->SetNLMMinEnergy(0, 10);
    ana->SetNLMMinEnergy(1, 6);
    ana->SetNLMMinEnergy(2, 6);
    
    // NLM cut, used in all, exclude clusters with more than 2 maxima
    ana->SetNLMCut(1, 2) ;
    
    //
    ana->SetMinDistanceToBadChannel(2, 4, 6);
    ana->SwitchOnSplitClusterDistToBad();
    ana->SetTimeCut(-1e10,1e10); // Open time cut
    
    AliCaloPID* caloPID = ana->GetCaloPID();
    
    // Track matching
    // tm = 1, fixed cuts
    caloPID->SetEMCALDEtaCut(0.025);
    caloPID->SetEMCALDPhiCut(0.030);
    
    // pT track dependent cuts
    if(tm > 1) caloPID->SwitchOnEMCTrackPtDepResMatching();

    caloPID->SetSplitWidthSigma(3); // cut at 3 sigma of the mean pi0 peak.
    
    if(!useSSIso)
    {
      printf("Do not apply SS cut on merged pi0 analysis \n");
      caloPID->SwitchOffSplitShowerShapeCut() ;
      ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_TM%d_",particle.Data(),opt.Data(),tm));
      ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenSS",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
      caloPID->SetClusterSplittingM02Cut(0.1,10);
    }
    else
    {
      caloPID->SetClusterSplittingM02Cut(0.3,4); // Do the selection in the analysis class and not in the PID method to fill SS histograms
      caloPID->SwitchOnSplitShowerShapeCut() ;
    }
    
    if(useAsy)
    {
      caloPID->SwitchOnSplitAsymmetryCut() ;
      ana->GetCaloPID()->SetSubClusterEnergyMinimum(0,2);
      ana->GetCaloPID()->SetSubClusterEnergyMinimum(1,0.5);
      ana->GetCaloPID()->SetSubClusterEnergyMinimum(2,0.5);
    }
    else
    {
      caloPID->SwitchOffSplitAsymmetryCut() ;
      if(!useSSIso)
      {
        ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_OpenAsy_TM%d_",particle.Data(),opt.Data(),tm));
        ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenSS_OpenAsy",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
      }
      else
      {
        ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenAsy_TM%d_",particle.Data(),opt.Data(),tm));
        ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenAsy",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
      }
    }
    
    // For Pi0 only if  SwitchOnSimpleSplitMassCut()
    caloPID->SetPi0MassRange(0.10, 0.18);
    caloPID->SetEtaMassRange(0.50, 0.60);
    caloPID->SetPhotonMassRange(0.00, 0.08);
    
    caloPID->SetClusterSplittingMinNCells(6);
    
    //caloPID->SetSplitEnergyFractionMinimum(0, 0.95);
    //caloPID->SetSplitEnergyFractionMinimum(1, 0.95);
    //caloPID->SetSplitEnergyFractionMinimum(2, 0.8);
    
    if(col=="PbPb" || kAnaGammaHadronCorr.Contains("150"))
    {
      caloPID->SetClusterSplittingMinNCells(4);
      //caloPID->SetPi0MassShiftHighECell(0.005);
    }
  }
  ///////////////////////////////////
  
  return  ana;
}

///
/// Configure the task doing the trigger particle isolation
///
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle,      Int_t  leading,
                                                    Int_t   partInCone,    Int_t   thresType,
                                                    Float_t cone,          Float_t pth,        Bool_t multi,
                                                    TString col,           Bool_t  simulation,
                                                    TString calorimeter,   Int_t   year,       Int_t tm,
                                                    Bool_t  printSettings, Int_t   debug                      )
{
  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
  
  ana->SetDebug(debug);
  
  ana->SetMinPt(5);
  ana->SetCalorimeter(calorimeter);
  
  ana->SwitchOffUEBandSubtractionHistoFill();
  ana->SwitchOffCellHistoFill() ;
  
  ana->SwitchOffLeadingOnly();
  ana->SwitchOffCheckNeutralClustersForLeading();
  if( leading > 0 )   ana->SwitchOnLeadingOnly();
  if( leading == 2 ||
     leading == 4)   ana->SwitchOnCheckNeutralClustersForLeading();
  
  // MC
  ana->SwitchOnPrimariesInConeSelection();
  ana->SwitchOnPrimariesPi0DecayStudy() ;
  
  if(particle.Contains("Photon"))
  {
    ana->SwitchOnDecayTaggedHistoFill() ;
    ana->SetNDecayBits(5);
    ana->SwitchOnSSHistoFill();
  }
  else
  {
    ana->SwitchOffSSHistoFill();
  }
  
  ana->SwitchOnPtTrigBinHistoFill();
  ana->SetNPtTrigBins(6);
  //ana->SetPtTrigLimits(0,8); ana->SetPtTrigLimits(1,12); ana->SetPtTrigLimits(2,16); ana->SetPtTrigLimits(3,25);
  
  ana->SwitchOnBackgroundBinHistoFill();
  ana->SetNBackgroundBins(11);
  //ana->SetBackgroundLimits(0,0); ana->SetBackgroundLimits(1,0.2); ana->SetBackgroundLimits(2,3); ana->SetBackgroundLimits(3,0.4);
  
  ana->SetFirstSMCoveredByTRD(6);
  
  if(!tm)  ana->SwitchOnTMHistoFill();
  else      ana->SwitchOffTMHistoFill();
  
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  ana->SwitchOnRealCaloAcceptance();
  ana->SwitchOnFiducialCut();
  
  if(calorimeter=="EMCAL")
  {
    // Avoid borders of EMCal
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.60, 86, 174) ;
    
  }
  
  // Same Eta as EMCal, cut in phi if EMCAL was triggering
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    //if(trigger.Contains("EMC"))
    //  ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
    //else
    ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;
  }
  
  // Branch AOD settings
  
  ana->SetInputAODName(Form("%sTrigger_%s",particle.Data(),kAnaGammaHadronCorr.Data()));
  ana->SetAODObjArrayName(Form("IC%sTrigger_%s_R%1.1f_ThMin%1.1f",particle.Data(),kAnaGammaHadronCorr.Data(),cone,pth));
  
  //
  // Do settings for main isolation cut class
  //
  AliIsolationCut * ic =  ana->GetIsolationCut();
  ic->SetDebug(debug);
  ic->SetParticleTypeInCone(partInCone);
  ic->SetICMethod(thresType);
  ic->SetPtFraction(0.1);
  ic->SetPtThreshold(0.5); // default, change in next lines
  ic->SetSumPtThreshold(1.0); // default, change in next lines
  
  if(cone > 0 && pth > 0)
  {
    ic->SetConeSize(cone);
    ic->SetPtThresholdMax(10000);
    
    if(thresType == AliIsolationCut::kPtThresIC)
    {
      printf("*** Iso *** PtThresMin = %1.1f GeV/c *** R = %1.1f ***\n",pth,cone);
      ic->SetPtThreshold(pth);
    }
    
    if(thresType == AliIsolationCut::kSumPtIC)
    {
      printf("*** Iso *** SumPtMin = %1.1f GeV/c *** R = %1.1f ***\n",pth,cone);
      ic->SetSumPtThreshold(pth);
    }
  }
  else
  {
    if(col=="pp")
    {
      ic->SetPtThreshold(0.5);
      ic->SetSumPtThreshold(1.0) ;
      ic->SetConeSize(0.4);
    }
    if(col=="PbPb")
    {
      ic->SetPtThreshold(3.);
      ic->SetSumPtThreshold(3.0) ;
      ic->SetConeSize(0.3);
    }
  }
  
  
  // Do or not do isolation with previously produced AODs.
  // No effect if use of SwitchOnSeveralIsolation()
  ana->SwitchOffReIsolation();
  
  // Multiple IC
  if(multi)
  {
    ic->SetConeSize(1.);    // Take all for first iteration
    ic->SetPtThreshold(100);// Take all for first iteration
    ana->SwitchOnSeveralIsolation() ;
    ana->SetAODObjArrayName(Form("MultiIC%sTM%d",particle.Data(),tm));
    
    ana->SetNCones(3);
    ana->SetNPtThresFrac(2);
    ana->SetConeSizes(0,0.3);       ana->SetConeSizes(1,0.4);       ana->SetConeSizes(2,0.5);
    ana->SetPtThresholds(0, 0.5);   ana->SetPtThresholds(1, 1);     ana->SetPtThresholds(2, 1.5);  ana->SetPtThresholds(3, 2);
    ana->SetPtFractions (0, 0.05) ; ana->SetPtFractions (1, 0.1);   ana->SetPtFractions (2, 0.2) ;  ana->SetPtFractions (3, 0.3) ;
    ana->SetSumPtThresholds(0, 0.5) ; ana->SetSumPtThresholds(1, 1) ; ana->SetSumPtThresholds(2, 1.5);  ana->SetSumPtThresholds(3, 2)  ;
    //ana->SetPtThresholds(0, 0.5);
    
    ana->SwitchOffTMHistoFill();
    ana->SwitchOffSSHistoFill();
  }
  else
    ana->SwitchOffSeveralIsolation() ;
  
  AliCaloPID* caloPID = ana->GetCaloPID();
 
  // Track matching 
  // tm = 1, fixed cuts
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
  
  // pT track dependent cuts
  if(tm > 1) caloPID->SwitchOnEMCTrackPtDepResMatching();

  //Set Histograms name tag, bins and ranges
  
  if(!multi) ana->AddToHistogramsName(Form("AnaIsol%s_TM%d_",particle.Data(),tm));
  else       ana->AddToHistogramsName(Form("AnaMultiIsol%s_TM%d_",particle.Data(),tm));
  
  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  if(printSettings) ic ->Print("");
  
  return ana;
}


///
/// Configure the task doing the trigger particle hadron correlation
///
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle,      Int_t   leading,
                                                                    Bool_t  bIsolated,     Float_t shshMax,
                                                                    Int_t partInCone,      Int_t   thresType,
                                                                    Float_t cone,          Float_t pth,      Bool_t mixOn,
                                                                    TString col,           Bool_t  simulation,
                                                                    TString calorimeter,   Int_t   year,     Int_t tm,
                                                                    Bool_t  printSettings, Int_t   debug                      )
{
  AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
  
  ana->SetTriggerPtRange(5,100);
  ana->SetAssociatedPtRange(0.2,100);
  ana->SetDeltaPhiCutRange  (TMath::DegToRad()*120.,TMath::DegToRad()*240.);
  
  // Underlying event
  ana->SetUeDeltaPhiCutRange(TMath::DegToRad()*60. ,TMath::DegToRad()*120.);
  ana->SwitchOnSeveralUECalculation();
  
  ana->SwitchOffAbsoluteLeading();  // Select trigger leading particle of all the selected tracks
  ana->SwitchOffNearSideLeading();  // Select trigger leading particle of all the particles at +-90 degrees, default
  ana->SwitchOffCheckNeutralClustersForLeading();
  
  if(leading >  0 && leading <  3 ) ana->SwitchOnAbsoluteLeading();
  if(leading >  2 )                 ana->SwitchOnNearSideLeading();
  if(leading == 2 || leading == 4 ) ana->SwitchOnCheckNeutralClustersForLeading();
  
  ana->SwitchOffFillPtImbalancePerPtABinHistograms();
  ana->SwitchOffCorrelationVzBin() ;
  ana->SwitchOffFillEtaGapHistograms();
  
  ana->SwitchOffFillHighMultiplicityHistograms();
  
  ana->SwitchOffPi0TriggerDecayCorr();
  if(particle.Contains("Photon"))
  {
    ana->SwitchOnDecayTriggerDecayCorr();
    ana->SetNDecayBits(5);
    printf("**** SET M02 limits in correlation task *** \n");
    ana->SetM02Cut(0.10,shshMax);
    ana->SwitchOnInvariantMassHistograms();
    ana->SwitchOnBackgroundBinsPtInConeHistograms();
  }
  
  ana->SetMCGenType(0,7);
  
  ana->SwitchOffLeadHadronSelection(); // Open cuts, just fill histograms
  ana->SwitchOnFillLeadHadronHistograms();
  ana->SwitchOnBackgroundBinsPtInConeHistograms();
  ana->SwitchOnBackgroundBinsTaggedDecayPtInConeHistograms();
  ana->SetLeadHadronPhiCut(TMath::DegToRad()*130, TMath::DegToRad()*230.);
  ana->SetLeadHadronPtCut(0.5, 1000);
  
  // if triggering on PHOS and EMCAL is on
  ana->SwitchOffNeutralCorr(); // Do only correlation with TPC
  //ana->SetPi0AODBranchName("Pi0EMCAL_TrigEMC7_Cl_TM1");
  
  ana->SwitchOffHMPIDCorrelation();
  
  ana->SwitchOffFillBradHistograms();
  
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  ana->SetNAssocPtBins(8);
  ana->SetAssocPtBinLimit(0, 1) ;
  ana->SetAssocPtBinLimit(1, 2) ;
  ana->SetAssocPtBinLimit(2, 3) ;
  ana->SetAssocPtBinLimit(3, 4) ;
  ana->SetAssocPtBinLimit(4, 5) ;
  ana->SetAssocPtBinLimit(5, 8) ;
  ana->SetAssocPtBinLimit(6, 10) ;
  ana->SetAssocPtBinLimit(7, 100);
  
  ana->SelectIsolated(bIsolated); // do correlation with isolated photons
  
  // Mixing with own pool
  if(mixOn)
  {
    ana->SwitchOnOwnMix();
    ana->SwitchOnFillNeutralInMixedEvent();
    
    if(bIsolated)
    {
      //Do settings for main isolation cut class
      AliIsolationCut * ic =  ana->GetIsolationCut();
      ic->SetDebug(debug);
      
      if(cone >0 && pth > 0)
      {
        printf("*** Correl *** PtThres = %1.1f GeV/c *** R = %1.1f ***\n",pth,cone);
        ic->SetPtThreshold(pth);
        ic->SetConeSize(cone);
      }
      else
      {
        if(col=="pp")
        {
          ic->SetPtThreshold(0.5);
          ic->SetConeSize(0.4);
        }
        if(col=="PbPb")
        {
          ic->SetPtThreshold(3.);
          //ic->SetPtThreshold(1.);
          ic->SetConeSize(0.3);
        }
      }
      
      ic->SetPtFraction(0.1);
      ic->SetSumPtThreshold(1.0) ;
      ic->SetParticleTypeInCone(partInCone);
      ic->SetICMethod(thresType);
    }
  }
  else
    ana->SwitchOffOwnMix();
  
  ana->SetNZvertBin(20);
  
  if(col=="pp")
  {
    ana->SetNMaxEvMix(100);
    ana->SwitchOnTrackMultBins();
    ana->SetNTrackMultBin(10);
    ana->SetNRPBin(1);
  }
  else
  {
    ana->SetNMaxEvMix(10);
    ana->SwitchOffTrackMultBins(); // centrality bins
    ana->SetNCentrBin(12);
    ana->SetNRPBin(3);
    if(kAnaGammaHadronCorr.Contains("60_90"))
    {
      printf("*** Set mixing for peripheral\n");
      ana->SetNMaxEvMix(50);
      ana->SetNCentrBin(2);
    }
  }
  
  ana->SwitchOnFiducialCut();
  
  if(calorimeter=="EMCAL")
  {
    // Avoid borders of EMCal, same as for isolation
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;
    
  }
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%sTrigger_%s",particle.Data(),kAnaGammaHadronCorr.Data()));
  ana->SetAODObjArrayName(Form("%sHadronCorrIso%dTrigger_%s",particle.Data(),bIsolated,kAnaGammaHadronCorr.Data()));
  //ana->SetAODNamepTInConeHisto(Form("IC%s_%s_R%1.1f_ThMin%1.1f"           ,particle.Data(),kAnaGammaHadronCorr.Data(),cone,pth));
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("Ana%sHadronCorr_Iso%d_TM%d_",particle.Data(),bIsolated,tm));
  
  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  return ana;
}

///
/// Configure the task doing the selected tracks checking
///
AliAnaChargedParticles* ConfigureChargedAnalysis( Bool_t simulation, Bool_t printSettings, Int_t   debug )
{
  AliAnaChargedParticles *ana = new AliAnaChargedParticles();
  
  ana->SetDebug(debug);
  
  // selection cuts
  
  ana->SetMinPt(0.2);
  ana->SwitchOnFiducialCut();
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ; //more restrictive cut in reader and after in isolation
  
  ana->SwitchOffFillVertexBC0Histograms();
  
  ana->SwitchOffFillPileUpHistograms();
  ana->SwitchOffFillTrackBCHistograms();
  
  // Branch AOD settings
  
  ana->SetOutputAODName(Form("HadronTrigger_%s",kAnaGammaHadronCorr.Data()));
  ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaHadrons_");
  
  SetAnalysisCommonParameters(ana,"CTS",2012,"pp",simulation,printSettings,debug); // see method below
  
  return ana;
}

///
/// Configure the task doing standard calorimeter QA
///
AliAnaCalorimeterQA* ConfigureQAAnalysis(TString col,           Bool_t  simulation,
                                         TString calorimeter,   Int_t   year,
                                         Bool_t  printSettings, Int_t   debug      )
{
  AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
  
  ana->SetCalorimeter(calorimeter);
  
  ana->SetTimeCut(-1e10,1e10); // Open time cut
  
  ana->SwitchOnCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
  
  ana->SwitchOnRealCaloAcceptance();
  
  ana->SwitchOffFiducialCut();
  ana->SwitchOffFillAllTH3Histogram();
  ana->SwitchOffFillAllPositionHistogram();
  ana->SwitchOffFillAllPositionHistogram2();
  ana->SwitchOffStudyBadClusters() ;
  ana->SwitchOffFillAllCellTimeHisto() ;
  
  ana->SwitchOnFillAllTrackMatchingHistogram();
  
  ana->AddToHistogramsName("QA_"); // Begining of histograms name
  
  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  return ana;
}

///
/// Configure the task filling generated particle kinematics histograms
///
AliAnaGeneratorKine* ConfigureGenKineAnalysis(Int_t   thresType,     Float_t cone,
                                              Float_t pth,
                                              TString col,           Bool_t  simulation,
                                              TString calorimeter,   Int_t   year,
                                              Bool_t  printSettings, Int_t   debug       )
{
  AliAnaGeneratorKine *ana = new AliAnaGeneratorKine();
  
  // Trigger detector, acceptance and pT cut
  ana->SetTriggerDetector(calorimeter);
  ana->SetMinPt(2); // Trigger photon, pi0 minimum pT
  ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.4, 90, 170);
  
  // Particles associated to trigger or isolation cone acceptance and pT cut
  ana->SetCalorimeter(calorimeter);
  ana->SetMinChargedPt(0.2);
  ana->SetMinNeutralPt(0.3);
  ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.65, 81, 179);
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360);
  
  // Isolation paramters
  AliIsolationCut * ic = ana->GetIsolationCut();
  ic->SetDebug(debug);
  ic->SetPtThreshold(pth);
  ic->SetConeSize(cone);
  ic->SetSumPtThreshold(1.0) ;
  ic->SetICMethod(thresType);
  
  ana->AddToHistogramsName("AnaGenKine_");
  
  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
  return ana;
}


///
/// Set common histograms binning
/// plus other analysis common settings like TRD covered super modules
/// the activation of the MC dedicated histograms and the activation of
/// the debug mode
///
void SetAnalysisCommonParameters(AliAnaCaloTrackCorrBaseClass* ana,
                                 TString calorimeter,   Int_t  year,
                                 TString col,           Bool_t simulation,
                                 Bool_t  printSettings, Int_t  debug)
{
  //
  // Histograms ranges
  //
  AliHistogramRanges* histoRanges = ana->GetHistogramRanges();
  
  histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
  
  if(calorimeter=="EMCAL")
  {
    if ( year == 2010 )
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
      histoRanges->SetHistoXRangeAndNBins(-230,90,120); // QA
      histoRanges->SetHistoYRangeAndNBins(370,450,40);  // QA
    }
    else if ( year < 2014 )
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
      histoRanges->SetHistoXRangeAndNBins(-460,90,200); // QA
      histoRanges->SetHistoYRangeAndNBins(100,450,100); // QA    
    }
    else // Run2
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 329*TMath::DegToRad(), 250) ;
      histoRanges->SetHistoXRangeAndNBins(-460,460,230); // QA
      histoRanges->SetHistoYRangeAndNBins(-450,450,225); // QA
    }
    
    histoRanges->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  }
  else if(calorimeter=="PHOS")
  {
    histoRanges->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 320*TMath::DegToRad(), 70) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
  }
  else if(calorimeter=="CTS")
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 4.9, 500);
  
  // Invariant mass histo
  histoRanges->SetHistoMassRangeAndNBins(0., 1., 200) ;
  histoRanges->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  
  // check if time calibration is on
  //histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  histoRanges->SetHistoTimeRangeAndNBins(-400.,400,400);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistodRRangeAndNBins(0.,0.15,150);//QA
  
  // QA, electron, charged
  histoRanges->SetHistoPOverERangeAndNBins(0,2.,200);
  histoRanges->SetHistodEdxRangeAndNBins(0.,200.,200);
  
  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  histoRanges->SetHistoZRangeAndNBins(-350,350,175);
  histoRanges->SetHistoRRangeAndNBins(430,460,30);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  
  // QA, correlation
  if(col=="PbPb")
  {
    histoRanges->SetHistoNClusterCellRangeAndNBins(0,100,100);
    histoRanges->SetHistoNClustersRangeAndNBins(0,500,50);
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,2000,200);
  }
  else
  {
    histoRanges->SetHistoNClusterCellRangeAndNBins(0,50,50);
    histoRanges->SetHistoNClustersRangeAndNBins(0,50,50);
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,200,200);
  }
  
  // xE, zT
  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,200);
  histoRanges->SetHistoHBPRangeAndNBins  (0.,10.,200);
  
  // Isolation
  histoRanges->SetHistoPtInConeRangeAndNBins(0, 50 , 250);
  histoRanges->SetHistoPtSumRangeAndNBins   (0, 100, 250);
  
  //
  // TRD SM
  //
  if     (year == 2011) ana->SetFirstSMCoveredByTRD( 6);
  else if(year == 2012 ||
          year == 2013) ana->SetFirstSMCoveredByTRD( 4);
  else                  ana->SetFirstSMCoveredByTRD(-1);
  
  //
  // MC histograms?
  //
  if(simulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else           ana->SwitchOffDataMC() ;
  
  //Set here generator name, default pythia
  //ana->GetMCAnalysisUtils()->SetMCGenerator("");
  
  //
  // Debug
  //
  if(printSettings) ana->Print("");
  
  ana->SetDebug(debug); // 10 for lots of messages
}
