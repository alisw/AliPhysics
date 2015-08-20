///// \file AddTaskGammaHadronCorrelation.C
///// \brief Configuration of gamma-hadron and pi0-hadron + isolation, correlation analysis
/////
/////
/////
///// Configuration macro for analysis of gamma-hadron and pi0-hadron correlation analysis
///// both pi0 and gamma isolated or not.
///// It does: 
/////   * Photon selection in calorimeter with AliAnaPhoton: Track matching, mild shower shape, NLM ... cuts
/////   * Tagging of photon candidate clusters as decay from pi0 or eta or from their side bands with AliAnaPi0EbE (4 calls, one pi0, one eta, one pi0 side band and one eta side band)
/////   * Tagging of photon candidate as isolated with AliAnaParticleIsolation
/////   * Correlation of the photon and charged tracks, twice, one without isolation condition, and other with isolation condition
/////   * Pi0 selection with the identification of merged clusters with AliAnaPi0EbE
/////   * Tagging of pi0 as isolated with AliAnaParticleIsolation
/////   * Correlation of the pi0 and charged tracks, twice, one without isolation condition, and other with isolation condition
/////   * Optionally, the QA tasks AliAnaCalorimeterQA and AliAnaChargedParticle are executed
/////
///// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
//
//// Set includes for compilation
//
////#include <TString.h>
////#include <TROOT.h>
////
////#include "AliLog.h"
////#include "AliAnalysisTaskCaloTrackCorrelation.h"
////#include "AliCaloTrackESDReader.h"
////#include "AliCaloTrackAODReader.h"
////#include "AliCalorimeterUtils.h"
////#include "AliAnaPhoton.h"
////#include "AliAnaPi0EbE.h"
////#include "AliHistogramRanges.h"
////#include "AliAnaParticleIsolation.h"
////#include "AliAnaParticleHadronCorrelation.h"
////#include "AliAnaChargedParticles.h"
////#include "AliAnaCalorimeterQA.h"
////#include "AliAnaGeneratorKine.h"
////#include "AliAnalysisTaskCaloTrackCorrelation.h"
////#include "AliAnaCaloTrackCorrMaker.h"
////#include "AliAnalysisManager.h"
////#include "AliInputEventHandler.h"
////#include "AliVTrack.h"
////
////// Declare methods for compilation
////
////AliCaloTrackReader  * ConfigureReader        (TString col,           Bool_t simulation, 
////                                              TString clustersArray, Bool_t tender, 
////                                              TString calorimeter,   Bool_t nonLinOn,
////                                              TString trigger,       Bool_t rejectEMCTrig, 
////                                              Int_t   minCen,        Int_t  maxCen,
////                                              Bool_t  mixOn,         Bool_t printSettings, 
////                                              Int_t   debug                                 );
////AliCalorimeterUtils * ConfigureCaloUtils     (TString col,           Bool_t simulation,
////                                              TString clustersArray, Bool_t tender, 
////                                              Bool_t  nonLinOn,      Int_t  year, 
////                                              Bool_t  printSettings, Int_t debug            );
////AliAnaPhoton        * ConfigurePhotonAnalysis(TString col,           Bool_t simulation, 
////                                              TString calorimeter,   Int_t  year,  Int_t tm,
////                                              Bool_t  printSettings, Int_t debug            );
////AliAnaPi0EbE        * ConfigurePi0EbEAnalysis(TString particle,      Int_t  analysis, 
////                                              Bool_t useSSIso,       Bool_t useAsy,
////                                              TString col,           Bool_t simulation, 
////                                              TString calorimeter,   Int_t  year, Int_t tm,
////                                              Bool_t  printSettings, Int_t debug            );
////AliAnaParticleIsolation* ConfigureIsolationAnalysis
////                                             (TString particle,      Int_t  leading,
////                                              Int_t partInCone, Int_t thresType,
////                                              Float_t cone,          Float_t pth,      Bool_t multi,
////                                              TString col,           Bool_t simulation, 
////                                              TString calorimeter,   Int_t  year,      Int_t tm,
////                                              Bool_t  printSettings, Int_t debug                    );
////AliAnaParticleHadronCorrelation * ConfigureHadronCorrelationAnalysis
////                                             (TString particle,      Int_t leading,
////                                              Bool_t  bIsolated,     Float_t shshMax,
////                                              Int_t partInCone,      Int_t   thresType,
////                                              Float_t cone,          Float_t pth,      Bool_t mixOn,
////                                              TString col,           Bool_t  simulation, 
////                                              TString calorimeter,   Int_t   year,     Int_t tm,
////                                              Bool_t  printSettings, Int_t   debug                     );
////AliAnaChargedParticles* ConfigureChargedAnalysis
////                                             (Bool_t simulation,     Bool_t  printSettings, Int_t   debug);
////AliAnaCalorimeterQA * ConfigureQAAnalysis    (TString col,           Bool_t  simulation, 
////                                              TString calorimeter,   Int_t   year,     
////                                              Bool_t  printSettings, Int_t   debug                     );
////AliAnaGeneratorKine* ConfigureGenKineAnalysis(Int_t   thresType,     Float_t cone,      Float_t pth,  
////                                              TString col,           Bool_t  simulation, 
////                                              TString calorimeter,   Int_t   year,     
////                                              Bool_t  printSettings, Int_t   debug                     );
////
////void SetAnalysisCommonParameters     (AliAnaCaloTrackCorrBaseClass* ana, 
////                                              TString calorimeter, Int_t year, 
////                                              TString col, Bool_t simulation, 
////                                              Bool_t printSettings, Int_t debug);
////UInt_t SetTriggerMaskFromName                (TString trigger);
//
///// Global name to be composed of the settings, used to set the AOD branch name
//TString kAnaGammaHadronCorr = "";
//
/////
///// Main method calling all the configuration 
///// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
/////
///// The options that can be passed to the macro are:
///// \param calorimeter : A string with he calorimeter used to measure the trigger particle
///// \param simulation : A bool identifying the data as simulation
///// \param year: The year the data was taken, used to configure some histograms
///// \param col: A string with the colliding system
///// \param trigger : A string with the trigger class, abbreviated, defined in method belowSetTriggerMaskFromName()
///// \param rejectEMCTrig : An int to reject EMCal triggered events with bad trigger: 0 no rejection, 1 old runs L1 bit, 2 newer runs L1 bit
///// \param clustersArray : A string with the array of clusters not being the default (default is empty string)
///// \param tender : A bool indicating if the tender was running before this analysis
///// \param nonLinOn : A bool to set the use of the non linearity correction
///// \param shshMax : A float setting the maximum value of the shower shape of the clusters for the correlation analysis
///// \param isoCone : A float setting the isolation cone size
///// \param isoPtTh : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
///// \param isoMethod : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
///// \param isoContent : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
///// \param leading : An int setting the type of leading particle selection: 0, select all;l 1: absolute  leading of charged; 2: absolute  leading of charged and neutral; 3: near side leading absolute of charged; 4: near side leading absolute of charged and neutral
///// \param tm : A bool to select neutral clusters as triggers
///// \param minCen : An int to select the minimum centrality, -1 means no selection
///// \param maxCen : An int to select the maximum centrality, -1 means no selection
///// \param mixOn : A bool to switch the correlation mixing analysis
///// \param qaAn : A bool to switch the calorimeter QA analysis
///// \param chargedAn : A bool to switch the selected tracks QA analysis
///// \param outputfile : A string to change the name of the histograms output file, default is AnalysisResults.root
///// \param printSettings : A bool to enable the print of the settings per task
///// \param debug : An int to define the debug level of all the tasks
/////
//AliAnalysisTaskCaloTrackCorrelation * AddTaskGammaHadronCorrelation
//(
// TString  calorimeter   = "EMCAL", 
// Bool_t   simulation    = kFALSE,
// Int_t    year          = 2011,
// TString  col           = "pp", 
// TString  trigger       = "EMC7", 
// Int_t    rejectEMCTrig = 0, 
// TString  clustersArray = "",
// Bool_t   tender        = kFALSE,
// Bool_t   nonLinOn      = kFALSE,
// Float_t  shshMax       = 0.27,
// Float_t  isoCone       = 0.4,
// Float_t  isoPtTh       = 0.5,
// Int_t    isoMethod     = AliIsolationCut::kPtThresIC,
// Int_t    isoContent    = AliIsolationCut::kNeutralAndCharged,
// Int_t    leading       = 0,
// Bool_t   tm            = kTRUE,
// Int_t    minCen        = -1,
// Int_t    maxCen        = -1,
// Bool_t   mixOn         = kTRUE,
// Bool_t   qaAn          = kFALSE,
// Bool_t   chargedAn     = kFALSE,
// TString  outputfile    = "",
// Bool_t   printSettings = kFALSE, 
// Int_t    debug         = 0  // Debug level
// )
//
//{  
//  // Get the pointer to the existing analysis manager via the static access method.
//  
//  printf("Passed settings:\n calorimeter <%s>, simulation <%d>, year <%d>,\n col <%s>, trigger <%s>, reject EMC <%d>, clustersArray <%s>, tender <%d>, non linearity <%d>\n shshMax <%2.2f>, isoCone <%2.2f>, isoPtTh <%2.2f>, isoMethod <%d>,isoContent <%d>,\n leading <%d>, tm <%d>, minCen <%d>, maxCen <%d>, mixOn <%d>,\n qaAn <%d>, chargedAn <%d>, outputfile <%s>, printSettings <%d>, debug <%d>\n", 
//         calorimeter.Data(),simulation,year,col.Data(),trigger.Data(), rejectEMCTrig, clustersArray.Data(),tender, nonLinOn, shshMax,
//         isoCone,isoPtTh,isoMethod,isoContent,leading,tm,
//         minCen,maxCen,mixOn,qaAn,chargedAn,outputfile.Data(),printSettings,debug);
//  
//  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
//  if (!mgr) 
//  {
//    ::Error("AddTask", "No analysis manager to connect to.");
//    return NULL;
//  }  
//  
//  // Check the analysis type using the event handlers connected to the analysis manager.
//  
//  if (!mgr->GetInputEventHandler()) 
//  {
//    ::Error("AddTask", "This task requires an input event handler");
//    return NULL;
//  }
//  
//  // Make sure the B field is enabled for track selection, some cuts need it
//  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);
//  
//  // Name for containers
//  
//  //kAnaGammaHadronCorr = Form("%s_Trig%s_Cl%s_TM%d_R%1.1f_Pt%1.1f",calorimeter.Data(), trigger.Data(),kClusterArray.Data(),tm,cone,pth);
//  kAnaGammaHadronCorr = Form("GammaHadron_%s_Trig%s_Col_%s_Year%d_Cl%s_Ten%d_TM%d_M02_%1.2f_IsoParam_C%1.2fPt%1.2fM%dPa%d_Lead%d_Mix%d",
//                             calorimeter.Data(),trigger.Data(),col.Data(),year,clustersArray.Data(),tender,
//                             tm, shshMax,isoCone,isoPtTh,isoMethod,isoContent,leading,mixOn);
//
//  if(col=="PbPb" && maxCen>=0) kAnaGammaHadronCorr+=Form("Cen%d_%d",minCen,maxCen);
//    
//  printf("<<<< NAME: %s >>>>>\n",kAnaGammaHadronCorr.Data());
//    
//  // #### Configure analysis ####
//    
//  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
// 
//  // General frame setting and configuration
//  maker->SetReader   ( ConfigureReader   (col,simulation,clustersArray,tender,calorimeter,nonLinOn,trigger,rejectEMCTrig,minCen,maxCen,mixOn,printSettings,debug) ); 
//  maker->SetCaloUtils( ConfigureCaloUtils(col,simulation,clustersArray,tender,nonLinOn,year,                                                 printSettings,debug) );
//                       
//  // Analysis tasks setting and configuration
//  Int_t n = 0;//Analysis number, order is important
//  
//  //
//  // Photon analysis
//  //
//  maker->AddAnalysis(ConfigurePhotonAnalysis(col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Photon cluster selection
//  
//  
//  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0"        , AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
//                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Pi0 event by event selection, invariant mass and photon tagging from decay
//  
//  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Eta"        , AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
//                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Eta event by event selection, invariant mass and photon tagging from decay
//  
//  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0SideBand", AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
//                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Pi0 out of peak event by event selection, and photon tagging from decay
//  
//  maker->AddAnalysis(ConfigurePi0EbEAnalysis("EtaSideBand", AliAnaPi0EbE::kIMCalo,kFALSE,kFALSE,
//                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Eta out of peak event by event selection, and photon tagging from decay
//
//  
//  maker->AddAnalysis(ConfigureIsolationAnalysis("Photon", leading, isoContent,isoMethod,isoCone,isoPtTh, kFALSE,
//                                                col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Photon isolation
//  
//  
//  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon", leading, kFALSE, shshMax, isoContent,isoMethod,isoCone,isoPtTh, mixOn,
//                                                        col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Gamma-hadron correlation
//  
//  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon", leading, kTRUE,  shshMax, isoContent,isoMethod,isoCone,isoPtTh, mixOn,
//                                                        col,simulation,calorimeter,year,tm,printSettings,debug) , n++); // Isolated gamma hadron correlation
//  
//  //
//  // Merged pi0 analysis
//  //
//  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kSSCalo,kTRUE,kTRUE,
//                                             col,simulation,calorimeter,year,tm,printSettings,debug), n++); // Pi0 event by event selection, cluster splitting
//  
//  maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SS", leading, isoContent,isoMethod,isoCone,isoPtTh, kFALSE,
//                                                col,simulation,calorimeter,year,tm,printSettings,debug), n++);          // Pi0 isolation, cluster splits
//  
//  
//  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS", leading, kFALSE, shshMax, isoContent,isoMethod,isoCone,isoPtTh, mixOn,
//                                                        col,simulation,calorimeter,year,tm,printSettings,debug), n++);  // Pi0-hadron correlation
//  
//  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS", leading, kTRUE,  shshMax, isoContent,isoMethod,isoCone,isoPtTh, mixOn,
//                                                        col,simulation,calorimeter,year,tm,printSettings,debug) , n++); // Isolated pi0-hadron correlation
//  
//  // Check the generated kinematics
//  if(simulation)  maker->AddAnalysis(ConfigureGenKineAnalysis(isoMethod,isoCone,isoPtTh,
//                                                              col,simulation,calorimeter,year,printSettings,debug), n++);
//  
//  // Charged analysis
//  if(chargedAn)   maker->AddAnalysis(ConfigureChargedAnalysis(simulation,printSettings,debug), n++); // track selection checks
//  
//  // Calo QA
//  if(qaAn)        maker->AddAnalysis(ConfigureQAAnalysis(col,simulation,calorimeter,year,printSettings,debug) , n++); 
//  
//  maker->SetAnaDebug(debug)  ;
//  
//  maker->SwitchOnHistogramsMaker()  ;
//  maker->SwitchOnAODsMaker()  ;
//  
//  if( simulation || !trigger.Contains("EMC") ) maker->SwitchOffDataControlHistograms();
//  
//  if(printSettings) maker->Print("");
//  
//  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, calorimeter.Data());
//
//  // Create task
//  
//  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("%s",kAnaGammaHadronCorr.Data()));
//  
//  task->SetDebugLevel(debug);
//
//  //task->SetBranches("ESD:AliESDRun.,AliESDHeader");
//  //task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
//
//  task->SetAnalysisMaker(maker);
//  
//  mgr->AddTask(task);
//  
//  //Create containers
//  
//  if(outputfile.Length()==0) outputfile = AliAnalysisManager::GetCommonFileName(); 
//  
//  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(kAnaGammaHadronCorr, TList::Class(), 
//                                                             AliAnalysisManager::kOutputContainer, 
//                                                             Form("%s",outputfile.Data()));
//	
//  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",kAnaGammaHadronCorr.Data()), TList::Class(), 
//                                                             AliAnalysisManager::kParamContainer, 
//                                                             "AnalysisParameters.root");
//  
//  // Create ONLY the output containers for the data produced by the task.
//  // Get and connect other common input/output containers via the manager as below
//  //==============================================================================
//  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
//  //if(!kData.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
//  mgr->ConnectOutput (task, 1, cout_pc);
//  mgr->ConnectOutput (task, 2, cout_cuts);
//    
//  if(!mixOn)
//  {    
//    UInt_t mask =  SetTriggerMaskFromName(trigger);
//    task->SelectCollisionCandidates(mask);
//  } 
//  
//  return task;
//}
//
/////
///// Configure the class handling the events and cluster/tracks filtering.
/////
//AliCaloTrackReader * ConfigureReader(TString col,           Bool_t simulation, 
//                                     TString clustersArray, Bool_t tender, 
//                                     TString calorimeter,   Bool_t nonLinOn,
//                                     TString trigger,       Bool_t rejectEMCTrig, 
//                                     Int_t   minCen,        Int_t  maxCen,
//                                     Bool_t  mixOn,         Bool_t printSettings, 
//                                     Int_t   debug                                )
//{
//  // Get the data type ESD or AOD
//  AliAnalysisManager * mgr = AliAnalysisManager::GetAnalysisManager();
//  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); 
//  
//  AliCaloTrackReader * reader = 0;
//  if     (inputDataType == "AOD") reader = new AliCaloTrackAODReader();
//  else if(inputDataType == "ESD") reader = new AliCaloTrackESDReader();
//  else printf("AliCaloTrackReader::ConfigureReader() - Data not known InputData=%s\n",inputDataType.Data());
//  
//  reader->SetDebug(debug);//10 for lots of messages
//  
//  //
//  // MC settings
//  //
//  // Check if kine stack is available, independent of request of simulation
//  Bool_t useKinematics = kFALSE;
//  useKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
//  
//  if(simulation)
//  {
//    if (!useKinematics && inputDataType=="AOD") useKinematics = kTRUE; //AOD primary should be available ...
//  }
//
//  if(useKinematics)
//  {
//    if(inputDataType == "ESD")
//    {
//      reader->SwitchOnStack();          
//      reader->SwitchOffAODMCParticles(); 
//    }
//    else if(inputDataType == "AOD")
//    {
//      reader->SwitchOffStack();          
//      reader->SwitchOnAODMCParticles(); 
//    }
//  }  
//  
//  // In case of Pythia pt Hard bin simulations (jet-jet, gamma-jet)
//  // reject some special events that bother the cross section
//  if(simulation)
//  {
//    // Event rejection cuts for jet-jet simulations, do not use in other
//    reader->SetPtHardAndJetPtComparison(kTRUE);
//    reader->SetPtHardAndJetPtFactor(4);
//    
//    reader->SetPtHardAndClusterPtComparison(kTRUE);
//    reader->SetPtHardAndClusterPtFactor(1.5);
//  }
//  
//  //------------------------
//  // Detector input filling
//  //------------------------
//
//  //Min cluster/track E
//  reader->SetEMCALEMin(0.3); 
//  reader->SetEMCALEMax(1000); 
//  reader->SetPHOSEMin(0.3);
//  reader->SetPHOSEMax(1000);
//  reader->SetCTSPtMin(0.2);
//  reader->SetCTSPtMax(1000);
//
//  reader->SwitchOffRecalculateVertexBC();
//  reader->SwitchOffVertexBCEventSelection();
//  
//  // Shower shape smearing
//  // Set it in the train configuration page not here for the moment
////  if(simulation)
////  {
////    reader->SwitchOffShowerShapeSmearing(); // Active only on MC, off by default
////    reader->SetShowerShapeSmearWidth(0.005);  
////  }
//
//  //
//  // Tracks
//  //
//  reader->SwitchOnCTS();
//
//  reader->SwitchOffUseTrackTimeCut();
//  reader->SetTrackTimeCut(0,50);
//  
//  reader->SwitchOnFiducialCut();
//  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;
//
//  reader->SwitchOffUseTrackDCACut();
//  //reader->SetTrackDCACut(0,0.0105);
//  //reader->SetTrackDCACut(1,0.035);
//  //reader->SetTrackDCACut(2,1.1);
//  
//  if(inputDataType=="ESD")
//  {
//    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
//    
//    //AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
//    //reader->SetTrackCuts(esdTrackCuts);
//    
//    AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
//    reader->SetTrackCuts(esdTrackCuts);
//    AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
//    reader->SetTrackComplementaryCuts(esdTrackCuts2);
//    
//    reader->SwitchOnConstrainTrackToVertex();
//  }
//  else if(inputDataType=="AOD")
//  {
//    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
//    reader->SwitchOnAODTrackSharedClusterSelection();
//    reader->SetTrackStatus(AliVTrack::kITSrefit);
//    
//    //reader->SwitchOnAODPrimaryTrackSelection(); // Used in preliminary results of QM from Nicolas and Xiangrong?
//    //reader->SwitchOnTrackHitSPDSelection();     // Check that the track has at least a hit on the SPD, not much sense to use for hybrid or TPC only tracks
//    //reader->SetTrackFilterMask(128);            // Filter bit, not mask, use if off hybrid, TPC only
//  }
//  
//  //
//  // Calorimeter
//  //
//  if(clustersArray == "" && !tender) 
//  {
//    printf("**************** Standard EMCAL clusters branch analysis **************** \n");
//    reader->SwitchOnClusterRecalculation();
//    // Check in ConfigureCaloUtils that the recalibration and bad map are ON 
//  }
//  else 
//  {
//    printf("**************** Input for analysis is Clusterizer %s **************** \n", clustersArray.Data());
//    reader->SetEMCALClusterListName(clustersArray);
//    reader->SwitchOffClusterRecalculation();
//  }  
//  
//  // Time cuts
//  reader->SwitchOffUseParametrizedTimeCut();
//  if(simulation) 
//  {
//    reader->SwitchOffUseEMCALTimeCut();
//    reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
//  }
//  else
//  {
//    reader->SwitchOnUseEMCALTimeCut();
//    reader->SetEMCALTimeCut(-25,20);
//  }
//
//  // CAREFUL
//  if(nonLinOn) reader->SwitchOnClusterELinearityCorrection();
//  else         reader->SwitchOffClusterELinearityCorrection();
//  
//  if(calorimeter == "EMCAL") 
//  {
//    reader->SwitchOnEMCALCells();  
//    reader->SwitchOnEMCAL();
//  }
//  
//  if(calorimeter == "PHOS") 
//  { // Should be on if QA is activated with correlation on
//    reader->SwitchOnPHOSCells();
//    reader->SwitchOnPHOS();
//  }
//  
//  //-----------------
//  // Event selection
//  //-----------------
//  
//  //if(!simulation) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
//  
//  // Event triggered by EMCal selection settings
//  reader->SwitchOffTriggerPatchMatching();
//  reader->SwitchOffBadTriggerEventsRemoval();
//  
//  if( rejectEMCTrig > 0 && !simulation && (trigger.Contains("EMC") || trigger.Contains("L")))
//  {
//    printf("=== Remove bad triggers === \n");
//    reader->SwitchOnTriggerPatchMatching();
//    reader->SwitchOnBadTriggerEventsRemoval();
//    
////    reader->SetTriggerPatchTimeWindow(8,9); // default values
////    if     (kRunNumber < 146861) reader->SetEventTriggerL0Threshold(3.);
////    else if(kRunNumber < 154000) reader->SetEventTriggerL0Threshold(4.);
////    else if(kRunNumber < 165000) reader->SetEventTriggerL0Threshold(5.5);
////    //redefine for other periods, triggers
////    
////    if(kRunNumber < 172000)
////    {
////      reader->SetEventTriggerL1Bit(4,5); // current LHC11 data
////      printf("\t Old L1 Trigger data format!\n");
////    }
////    else
////    {
////      reader->SetEventTriggerL1Bit(6,8); // LHC12-13 data
////      printf("\t Current L1 Trigger data format!\n");
////    }
//    
//    if(clustersArray != "" || tender)
//    {
//      printf("Trigger cluster calibration OFF\n");
//      reader->SwitchOffTriggerClusterTimeRecal() ;
//    }
//    
//  }
//
//  //reader->RejectFastClusterEvents() ;
//  
//  // For mixing with AliAnaParticleHadronCorrelation switch it off
//  reader->SwitchOnEventTriggerAtSE(); // on is default case
//  if(mixOn)
//  {
//    reader->SwitchOffEventTriggerAtSE();
//    UInt_t mask = SetTriggerMaskFromName(trigger);
//    reader->SetEventTriggerMask(mask); // Only for mixing and SwitchOffEventTriggerAtSE();
//    //reader->SetMixEventTriggerMask(AliVEvent::kMB); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
//    reader->SetMixEventTriggerMask(AliVEvent::kINT7); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
//
//    printf("---Trigger selection done in AliCaloTrackReader!!!\n");
//  }
//  
//  reader->SetZvertexCut(10.);               // Open cut
//  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
//  reader->SwitchOnRejectNoTrackEvents();
//
//  reader->SwitchOffV0ANDSelection() ;       // and besides v0 AND
//  reader->SwitchOffPileUpEventRejection();  // remove pileup by default off, apply it only for MB not for trigger
//  
//  if(col=="PbPb") 
//  {
//    // Centrality
//    reader->SetCentralityClass("V0M");
//    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
//    reader->SetCentralityBin(minCen,maxCen); // Accept all events, if not select range
//    reader->SwitchOnAcceptOnlyHIJINGLabels();
//    // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
//    reader->SetEventPlaneMethod("V0");
//  }
//  
//  if(printSettings) reader->Print("");
//  
//  return reader;
//}
//
/////
///// Configure the class handling the calorimeter clusters specific methods
/////
//AliCalorimeterUtils* ConfigureCaloUtils(TString col,           Bool_t simulation,
//                                        TString clustersArray, Bool_t tender, 
//                                        Bool_t  nonLinOn,      Int_t year, 
//                                        Bool_t  printSettings, Int_t   debug)
//{
//  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
//  
//  cu->SetDebug(debug);
//
//  // Remove clusters close to borders, at least max energy cell is 1 cell away 
//  cu->SetNumberOfCellsFromEMCALBorder(1);
//  cu->SetNumberOfCellsFromPHOSBorder (2);
//  
//  cu->SetNumberOfSuperModulesUsed(10);
//  
//  if     (year == 2010) cu->SetNumberOfSuperModulesUsed(4);
//  else if(year <= 2013) cu->SetNumberOfSuperModulesUsed(10);
//  else if(year >  2013) cu->SetNumberOfSuperModulesUsed(20);
//  else                  cu->SetNumberOfSuperModulesUsed(10);
//  
//  printf("xxx Number of SM set to <%d> xxx\n",cu->GetNumberOfSuperModulesUsed());
//  
//  // Search of local maxima in cluster
//  if(col=="pp")
//  {
//    cu->SetLocalMaximaCutE(0.1);
//    cu->SetLocalMaximaCutEDiff(0.03);
//  }
//  else 
//  {
//    cu->SetLocalMaximaCutE(0.2);
//    cu->SetLocalMaximaCutEDiff(0.03);
//  }
//  
//  cu->SwitchOffRecalculateClusterTrackMatching();
//  
//  cu->SwitchOnBadChannelsRemoval() ;
//  
//  // EMCAL settings
//
//  if(!simulation)
//    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
//  
//  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
//  
//  // calibrations
//  Bool_t calibEner = kFALSE;
//  Bool_t calibTime = kFALSE;
//  cu->SwitchOffRecalibration(); 
//  cu->SwitchOffRunDepCorrection();
//  
//  if( !tender )
//  {
//    cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
//    cu->SwitchOnRunDepCorrection();
//    
//    calibEner = kTRUE;
//    calibTime = kTRUE;
//  }
//  
//  if( simulation )
//  {
//    calibEner = kFALSE;
//    calibTime = kFALSE;
//
//    cu->SwitchOffRecalibration(); // Check the reader if it is taken into account during filtering
//    cu->SwitchOffRunDepCorrection();
//  }
//  
//  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
//  ConfigureEMCALRecoUtils(recou,
//                          simulation,                             
//                          kTRUE,      // exotic
//                          nonLinOn,   // Non linearity
//                          calibEner,  // E calib
//                          kTRUE,      // bad map
//                          calibTime); // time calib   
//  
//  if( calibTime ) recou->SetExoticCellDiffTimeCut(50);
//  
//  if( nonLinOn )  cu->SwitchOnCorrectClusterLinearity();
//    
//  printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
//  printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
//    
//  // PHOS 
//  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
//    
//  if(printSettings) cu->Print("");
//  
//  return cu;
//}
//
/////
///// Configure the task doing the first photon cluster selections
///// Basically the track matching, minor shower shape cut, NLM selection ...
/////
//AliAnaPhoton* ConfigurePhotonAnalysis(TString col,           Bool_t simulation, 
//                                      TString calorimeter,   Int_t year, Int_t tm,
//                                      Bool_t  printSettings, Int_t   debug)
//{
//  AliAnaPhoton *ana = new AliAnaPhoton();
//
//  // cluster selection cuts
//  
//  ana->SwitchOnRealCaloAcceptance();
//
//  ana->SwitchOffFiducialCut();
//
//  ana->SetCalorimeter(calorimeter);
//    
//  ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection
//
// //if(!simulation) ana->SwitchOnFillPileUpHistograms();
//  
//  if(tm) ana->SwitchOnTrackMatchRejection() ;
//  else   ana->SwitchOffTrackMatchRejection() ;
//  
//  ana->SwitchOnTMHistoFill() ;
//  
//  if(calorimeter == "PHOS")
//  {
//    ana->SetNCellCut(2);// At least 3 cells
//    ana->SetMinPt(0.3);
//    ana->SetMinDistanceToBadChannel(2, 4, 5);
//    ana->SetTimeCut(-1e10,1e10); // open cut
//  }
//  else 
//  {//EMCAL
//    ana->SetNCellCut(1);// At least 2 cells
//    ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
//    ana->SetMaxEnergy(100);
//    ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off
//    // restrict to less than 100 ns when time calibration is on
//    ana->SetMinDistanceToBadChannel(2, 4, 6);
//    
//    // NLM cut, used in all, exclude clusters with more than 2 maxima
//    // Not needed if M02 cut is already strong or clusterizer V2
//    ana->SetNLMCut(1, 2) ;    
//  }
//  
//  //PID cuts (shower shape)
//  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
//  AliCaloPID* caloPID = ana->GetCaloPID();
//  //Not used in bayesian
//  
//  // EMCAL
//  
//  //caloPID->SetEMCALLambda0CutMax(0.27);
//  caloPID->SetEMCALLambda0CutMax(10); // open, full shower shape needed for isolation studies 
//  caloPID->SetEMCALLambda0CutMin(0.10);
//  
//  // Track matching
//  caloPID->SetEMCALDEtaCut(0.025);
//  caloPID->SetEMCALDPhiCut(0.030);
//    
//  // PHOS
//  caloPID->SetPHOSDispersionCut(2.5);
//  caloPID->SetPHOSRCut(2.);
//  //if(kInputData=="AOD") caloPID->SetPHOSRCut(2000.); // Open cut since dX, dZ not stored
//
//  // Branch AOD settings
//  ana->SetOutputAODName(Form("PhotonTrigger_%s",kAnaGammaHadronCorr.Data()));
//  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
//  
//  //Set Histograms name tag, bins and ranges
//  ana->AddToHistogramsName(Form("AnaPhoton_TM%d_",tm));
//  
//  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug) ; // see method below
//  
//  if(ana->GetFirstSMCoveredByTRD() > 0)
//    printf(">>> Set first SM covered by TRD, SM=%d <<< year %d \n", ana->GetFirstSMCoveredByTRD(),year);
//  
//  // Number of particle type MC histograms
//  ana->FillNOriginHistograms (14); // 14 max
//  ana->FillNPrimaryHistograms(6); // 6 max
//      
//  return ana;
//}
//
/////
///// Configure the task doing the pi0 even by event selection (invariant mass or split)
///// and the cluster tagging as decay in different mass windows.
/////
//AliAnaPi0EbE* ConfigurePi0EbEAnalysis(TString particle,      Int_t  analysis, 
//                                      Bool_t useSSIso,       Bool_t useAsy,
//                                      TString col,           Bool_t simulation, 
//                                      TString calorimeter,   Int_t  year,  Int_t tm,
//                                      Bool_t  printSettings, Int_t debug           )
//{
//  // Configuration of pi0 event by event selection
//  
//  AliAnaPi0EbE *ana = new AliAnaPi0EbE();
//  
//  ana->SetAnalysisType((AliAnaPi0EbE::anaTypes)analysis);
//  TString opt = "";
//  if(analysis==AliAnaPi0EbE::kIMCaloTracks) opt = "Conv";
//  if(analysis==AliAnaPi0EbE::kSSCalo)       opt = "SS";
//
//  ana->SwitchOffAllNLMHistoFill();
//  ana->SwitchOffSelectedClusterHistoFill();
//
//  ana->SwitchOffFillWeightHistograms();
//  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
//  //if(!kTime && !simulation) ana->SwitchOnFillEMCALBCHistograms();
//
//  if(tm) ana->SwitchOnTrackMatchRejection() ;
//  else   ana->SwitchOffTrackMatchRejection() ;
//  ana->SwitchOffTMHistoFill() ;
//
//  ana->SetCalorimeter(calorimeter);
//  
//  // Branch AOD settings
//  ana->SetOutputAODName(Form("%s%sTrigger_%s",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
//  printf("***Out branch %s***\n",ana->GetOutputAODName().Data());
//  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
// 
//  if(analysis == AliAnaPi0EbE::kIMCaloTracks) ana->SetInputAODGammaConvName("PhotonsCTS");
//  
//  //Set Histograms name tag, bins and ranges
//  
//  ana->AddToHistogramsName(Form("Ana%s%sEbE_TM%d_",particle.Data(),opt.Data(),tm));
//
//  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
//
//  ///////////////////////////////////
//  if(analysis!=AliAnaPi0EbE::kSSCalo)
//  {
//    ana->SetInputAODName(Form("PhotonTrigger_%s",kAnaGammaHadronCorr.Data()));
//    
//    ana->SetM02CutForInvMass(0.1,0.35); // Loose SS cut
//
//    ana->SwitchOnSelectPairInIsolationCone();
//    ana->SetR(0.4);
//    ana->SetIsolationCandidateMinPt(5);
//
//    if(useSSIso)
//    {
//      ana->SwitchOnSelectIsolatedDecay();
//      ana->AddToHistogramsName(Form("Ana%s%sEbEIsoDecay_TM%d_",particle.Data(),opt.Data(),tm));
//      ana->SetOutputAODName(Form("%s%sIsoDecayTrigger_%s",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
//    }
//    
//    if(calorimeter=="EMCAL" && !simulation) ana->SetPairTimeCut(100);
//
//    AliNeutralMesonSelection *nms = ana->GetNeutralMesonSelection();
//    nms->SetParticle(particle);
//    
//    //****
//    nms->SetInvMassCutMaxParameters(0,0,0); // Overrule the setting in SetParticle for Pi0 option
//    //****
//    
//    // Tighten a bit mass cut with respect to default window
//    if(particle=="Pi0") nms->SetInvMassCutRange(0.110,0.160);
//    if(particle=="Eta") nms->SetInvMassCutRange(0.520,0.580);
//    
//    //if(!particle.Contains("SideBand")) nms->SwitchOnAngleSelection();
//    //else nms->SwitchOnAngleSelection();
//    
//    nms->SwitchOffAngleSelection();
//    
//    if(particle.Contains("Pi0SideBand")) // For pi0, do not consider left band
//      nms->SetSideBandCutRanges(-1,0,0.190,0.240);
//    
//    if(particle.Contains("EtaSideBand")) // For pi0, do not consider left band
//      nms->SetSideBandCutRanges(0.410,0.470,0.620,0.680);
//    
//    nms->KeepNeutralMesonSelectionHistos(kTRUE);
//    //nms->SetAngleMaxParam(2,0.2);
//    nms->SetHistoERangeAndNBins(0, 20, 80) ;
//    //nms->SetHistoIMRangeAndNBins(0, 1, 400);
//  }
//  else  ///////////////////////////////////
//  {
//    // cluster splitting settings
//    ana->SetMinEnergy(6);
//    ana->SetMaxEnergy(100.);
//    
//    ana->SetNLMMinEnergy(0, 10);
//    ana->SetNLMMinEnergy(1, 6);
//    ana->SetNLMMinEnergy(2, 6);
//    
//    // NLM cut, used in all, exclude clusters with more than 2 maxima
//    ana->SetNLMCut(1, 2) ;
//
//    //
//    ana->SetMinDistanceToBadChannel(2, 4, 6);
//    ana->SwitchOnSplitClusterDistToBad();
//    ana->SetTimeCut(-1e10,1e10); // Open time cut
//    
//    AliCaloPID* caloPID = ana->GetCaloPID();
//    
//    caloPID->SetSplitWidthSigma(3); // cut at 3 sigma of the mean pi0 peak.
//    
//    if(!useSSIso)
//    {
//      printf("Do not apply SS cut on merged pi0 analysis \n");
//      caloPID->SwitchOffSplitShowerShapeCut() ;
//      ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_TM%d_",particle.Data(),opt.Data(),tm));
//      ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenSS",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
//      caloPID->SetClusterSplittingM02Cut(0.1,10);
//    }
//    else
//    {
//      caloPID->SetClusterSplittingM02Cut(0.3,4); // Do the selection in the analysis class and not in the PID method to fill SS histograms
//      caloPID->SwitchOnSplitShowerShapeCut() ;
//    }
//    
//    if(useAsy)
//    {
//      caloPID->SwitchOnSplitAsymmetryCut() ;
//      ana->GetCaloPID()->SetSubClusterEnergyMinimum(0,2);
//      ana->GetCaloPID()->SetSubClusterEnergyMinimum(1,0.5);
//      ana->GetCaloPID()->SetSubClusterEnergyMinimum(2,0.5);
//    }
//    else
//    {
//      caloPID->SwitchOffSplitAsymmetryCut() ;
//      if(!useSSIso)
//      {
//        ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_OpenAsy_TM%d_",particle.Data(),opt.Data(),tm));
//        ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenSS_OpenAsy",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
//      }
//      else
//      {
//        ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenAsy_TM%d_",particle.Data(),opt.Data(),tm));
//        ana->SetOutputAODName(Form("%s%sTrigger_%s_OpenAsy",particle.Data(), opt.Data(), kAnaGammaHadronCorr.Data()));
//      }
//    }
//    
//    // For Pi0 only if  SwitchOnSimpleSplitMassCut()
//    caloPID->SetPi0MassRange(0.10, 0.18);
//    caloPID->SetEtaMassRange(0.50, 0.60);
//    caloPID->SetPhotonMassRange(0.00, 0.08);
//    
//    caloPID->SetClusterSplittingMinNCells(6);
//    
//    //caloPID->SetSplitEnergyFractionMinimum(0, 0.95);
//    //caloPID->SetSplitEnergyFractionMinimum(1, 0.95);
//    //caloPID->SetSplitEnergyFractionMinimum(2, 0.8);
//    
//    if(col=="PbPb" || kAnaGammaHadronCorr.Contains("150"))
//    {
//      caloPID->SetClusterSplittingMinNCells(4);
//      //caloPID->SetPi0MassShiftHighECell(0.005);
//    }
//  }
//  ///////////////////////////////////
//
//  return  ana;
//}
//
/////
///// Configure the task doing the trigger particle isolation
/////
//AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle,      Int_t  leading, 
//                                                    Int_t   partInCone,    Int_t   thresType,
//                                                    Float_t cone,          Float_t pth,        Bool_t multi,
//                                                    TString col,           Bool_t  simulation, 
//                                                    TString calorimeter,   Int_t   year,       Int_t tm,
//                                                    Bool_t  printSettings, Int_t   debug                      )
//{
//  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
//
//  ana->SetMinPt(5);
//  ana->SetCalorimeter(calorimeter);
//
//  ana->SwitchOffUEBandSubtractionHistoFill();
//  ana->SwitchOffCellHistoFill() ;
//
//  ana->SwitchOffLeadingOnly();
//  ana->SwitchOffCheckNeutralClustersForLeading();
//  if( leading > 0 )   ana->SwitchOnLeadingOnly();
//  if( leading == 2 || 
//      leading == 4)   ana->SwitchOnCheckNeutralClustersForLeading();
//    
//  // MC
//  ana->SwitchOnPrimariesInConeSelection();
//  ana->SwitchOnPrimariesPi0DecayStudy() ;
//
//  if(particle.Contains("Photon"))
//  {
//    ana->SwitchOnDecayTaggedHistoFill() ;
//    ana->SetNDecayBits(5);
//    ana->SwitchOnSSHistoFill();
//  }
//  else
//  {
//    ana->SwitchOffSSHistoFill();
//  }
//
//  ana->SwitchOnPtTrigBinHistoFill();  
//  ana->SetNPtTrigBins(6);
//  //ana->SetPtTrigLimits(0,8); ana->SetPtTrigLimits(1,12); ana->SetPtTrigLimits(2,16); ana->SetPtTrigLimits(3,25);
//  
//  ana->SwitchOnBackgroundBinHistoFill();
//  ana->SetNBackgroundBins(11);
//  //ana->SetBackgroundLimits(0,0); ana->SetBackgroundLimits(1,0.2); ana->SetBackgroundLimits(2,3); ana->SetBackgroundLimits(3,0.4);
//    
//  if(!tm)  ana->SwitchOnTMHistoFill();
//  else      ana->SwitchOffTMHistoFill();
//  
//  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
//
//  ana->SwitchOnRealCaloAcceptance();
//  ana->SwitchOnFiducialCut();
//
//  if(calorimeter=="EMCAL")
//  {
//    // Avoid borders of EMCal
//    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.60, 86, 174) ;
//
//  }
//  
//  // Same Eta as EMCal, cut in phi if EMCAL was triggering
//  if(particle=="Hadron"  || particle.Contains("CTS"))
//  {
//    //if(trigger.Contains("EMC"))
//    //  ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
//    //else
//    ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;
//  }
//  
//  // Branch AOD settings
//
//  ana->SetInputAODName(Form("%sTrigger_%s",particle.Data(),kAnaGammaHadronCorr.Data()));
//  ana->SetAODObjArrayName(Form("IC%sTrigger_%s_R%1.1f_ThMin%1.1f",particle.Data(),kAnaGammaHadronCorr.Data(),cone,pth));
//  
//  //
//  // Do settings for main isolation cut class
//  //
//  AliIsolationCut * ic =  ana->GetIsolationCut();
//  ic->SetDebug(debug);
//  ic->SetParticleTypeInCone(partInCone);
//  ic->SetICMethod(thresType);
//  ic->SetPtFraction(0.1);
//  ic->SetPtThreshold(0.5); // default, change in next lines
//  ic->SetSumPtThreshold(1.0); // default, change in next lines
//
//  if(cone > 0 && pth > 0)
//  {
//    ic->SetConeSize(cone);
//    ic->SetPtThresholdMax(10000);
//
//    if(thresType == AliIsolationCut::kPtThresIC)
//    {
//      printf("*** Iso *** PtThresMin = %1.1f GeV/c *** R = %1.1f ***\n",pth,cone);
//      ic->SetPtThreshold(pth);
//    }
//    
//    if(thresType == AliIsolationCut::kSumPtIC)
//    {
//      printf("*** Iso *** SumPtMin = %1.1f GeV/c *** R = %1.1f ***\n",pth,cone);
//      ic->SetSumPtThreshold(pth);
//    }
//  }
//  else
//  {
//    if(col=="pp")
//    {
//      ic->SetPtThreshold(0.5);
//      ic->SetSumPtThreshold(1.0) ;
//      ic->SetConeSize(0.4);
//    }
//    if(col=="PbPb")
//    {
//      ic->SetPtThreshold(3.);
//      ic->SetSumPtThreshold(3.0) ;
//      ic->SetConeSize(0.3);
//    }
//  }
//  
// 
//  // Do or not do isolation with previously produced AODs.
//  // No effect if use of SwitchOnSeveralIsolation()
//  ana->SwitchOffReIsolation();
//  
//  // Multiple IC
//  if(multi)
//  {
//    ic->SetConeSize(1.);    // Take all for first iteration
//    ic->SetPtThreshold(100);// Take all for first iteration
//    ana->SwitchOnSeveralIsolation() ;
//    ana->SetAODObjArrayName(Form("MultiIC%sTM%d",particle.Data(),tm));
//    
//    ana->SetNCones(3);
//    ana->SetNPtThresFrac(2);
//    ana->SetConeSizes(0,0.3);       ana->SetConeSizes(1,0.4);       ana->SetConeSizes(2,0.5);
//    ana->SetPtThresholds(0, 0.5);   ana->SetPtThresholds(1, 1);     ana->SetPtThresholds(2, 1.5);  ana->SetPtThresholds(3, 2);
//    ana->SetPtFractions (0, 0.05) ; ana->SetPtFractions (1, 0.1);   ana->SetPtFractions (2, 0.2) ;  ana->SetPtFractions (3, 0.3) ;
//    ana->SetSumPtThresholds(0, 0.5) ; ana->SetSumPtThresholds(1, 1) ; ana->SetSumPtThresholds(2, 1.5);  ana->SetSumPtThresholds(3, 2)  ;
//    //ana->SetPtThresholds(0, 0.5);
//    
//    ana->SwitchOffTMHistoFill();
//    ana->SwitchOffSSHistoFill();
//  }
//  else
//    ana->SwitchOffSeveralIsolation() ;
//  
//  AliCaloPID* caloPID = ana->GetCaloPID();
//  caloPID->SetEMCALDEtaCut(0.025);
//  caloPID->SetEMCALDPhiCut(0.030);
//  
//  //Set Histograms name tag, bins and ranges
//  
//  if(!multi) ana->AddToHistogramsName(Form("AnaIsol%s_TM%d_",particle.Data(),tm));
//  else       ana->AddToHistogramsName(Form("AnaMultiIsol%s_TM%d_",particle.Data(),tm));
//   
//  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
//  
//  if(particle=="Hadron"  || particle.Contains("CTS"))
//  {
//    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
//    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
//  }
//  
//  if(printSettings) ic ->Print("");
//  
//  return ana;
//}
//
//
/////
///// Configure the task doing the trigger particle hadron correlation
/////
//AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle,      Int_t   leading, 
//                                                                    Bool_t  bIsolated,     Float_t shshMax,
//                                                                    Int_t partInCone,      Int_t   thresType,
//                                                                    Float_t cone,          Float_t pth,      Bool_t mixOn,
//                                                                    TString col,           Bool_t  simulation, 
//                                                                    TString calorimeter,   Int_t   year,     Int_t tm,
//                                                                    Bool_t  printSettings, Int_t   debug                      )
//{
//  AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
//  
//  ana->SetTriggerPtRange(5,100);
//  ana->SetAssociatedPtRange(0.2,100);
//  ana->SetDeltaPhiCutRange  (TMath::DegToRad()*120.,TMath::DegToRad()*240.);
//  
//  // Underlying event
//  ana->SetUeDeltaPhiCutRange(TMath::DegToRad()*60. ,TMath::DegToRad()*120.);
//  ana->SwitchOnSeveralUECalculation();
//  
//  ana->SwitchOffAbsoluteLeading();  // Select trigger leading particle of all the selected tracks
//  ana->SwitchOffNearSideLeading();  // Select trigger leading particle of all the particles at +-90 degrees, default
//  ana->SwitchOffCheckNeutralClustersForLeading();
//  
//  if(leading >  0 && leading <  3 ) ana->SwitchOnAbsoluteLeading();
//  if(leading >  2 )                 ana->SwitchOnNearSideLeading();
//  if(leading == 2 || leading == 4 ) ana->SwitchOnCheckNeutralClustersForLeading();
//  
//  ana->SwitchOffFillPtImbalancePerPtABinHistograms();
//  ana->SwitchOffCorrelationVzBin() ;
//  ana->SwitchOffFillEtaGapHistograms();
//
//  ana->SwitchOffFillHighMultiplicityHistograms();
//  
//  ana->SwitchOffPi0TriggerDecayCorr();
//  if(particle.Contains("Photon"))
//  {
//    ana->SwitchOnDecayTriggerDecayCorr();
//    ana->SetNDecayBits(5);
//    printf("**** SET M02 limits in correlation task *** \n");
//    ana->SetM02Cut(0.10,shshMax);
//    ana->SwitchOnInvariantMassHistograms();
//    ana->SwitchOnBackgroundBinsPtInConeHistograms();
//  }
//  
//  ana->SetMCGenType(0,7);
//  
//  ana->SwitchOffLeadHadronSelection(); // Open cuts, just fill histograms
//  ana->SwitchOnFillLeadHadronHistograms();
//  ana->SwitchOnBackgroundBinsTaggedDecayPtInConeHistograms();
//  ana->SetLeadHadronPhiCut(TMath::DegToRad()*130, TMath::DegToRad()*230.);
//  ana->SetLeadHadronPtCut(0.5, 1000);
//  
//  // if triggering on PHOS and EMCAL is on
//  ana->SwitchOffNeutralCorr(); // Do only correlation with TPC
//  //ana->SetPi0AODBranchName("Pi0EMCAL_TrigEMC7_Cl_TM1");
//  
//  ana->SwitchOffHMPIDCorrelation();
//  
//  ana->SwitchOffFillBradHistograms();
//
//  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
//  
//  ana->SetNAssocPtBins(8);
//  ana->SetAssocPtBinLimit(0, 1) ;
//  ana->SetAssocPtBinLimit(1, 2) ;
//  ana->SetAssocPtBinLimit(2, 3) ;
//  ana->SetAssocPtBinLimit(3, 4) ;
//  ana->SetAssocPtBinLimit(4, 5) ;
//  ana->SetAssocPtBinLimit(5, 8) ;
//  ana->SetAssocPtBinLimit(6, 10) ;
//  ana->SetAssocPtBinLimit(7, 100);
//  
//  ana->SelectIsolated(bIsolated); // do correlation with isolated photons
//  
//  // Mixing with own pool
//  if(mixOn)
//  {
//    ana->SwitchOnOwnMix();
//    ana->SwitchOnFillNeutralInMixedEvent();
//    
//    if(bIsolated)
//    {
//      //Do settings for main isolation cut class
//      AliIsolationCut * ic =  ana->GetIsolationCut();
//      ic->SetDebug(debug);
//      
//      if(cone >0 && pth > 0)
//      {
//        printf("*** Correl *** PtThres = %1.1f GeV/c *** R = %1.1f ***\n",pth,cone);
//        ic->SetPtThreshold(pth);
//        ic->SetConeSize(cone);
//      }
//      else
//      {
//        if(col=="pp")
//        {
//          ic->SetPtThreshold(0.5);
//          ic->SetConeSize(0.4);
//        }
//        if(col=="PbPb")
//        {
//          ic->SetPtThreshold(3.);
//          //ic->SetPtThreshold(1.);
//          ic->SetConeSize(0.3);
//        }
//      }
//      
//      ic->SetPtFraction(0.1);
//      ic->SetSumPtThreshold(1.0) ;
//      ic->SetParticleTypeInCone(partInCone);
//      ic->SetICMethod(thresType);
//    }
//  }
//  else
//    ana->SwitchOffOwnMix();
//  
//  ana->SetNZvertBin(20);
//  
//  if(col=="pp")
//  {
//    ana->SetNMaxEvMix(100);
//    ana->SwitchOnTrackMultBins();
//    ana->SetNTrackMultBin(10);
//    ana->SetNRPBin(1);
//  }
//  else
//  {
//    ana->SetNMaxEvMix(10);
//    ana->SwitchOffTrackMultBins(); // centrality bins
//    ana->SetNCentrBin(12);
//    ana->SetNRPBin(3);
//    if(kAnaGammaHadronCorr.Contains("60_90"))
//    {
//      printf("*** Set mixing for peripheral\n");
//      ana->SetNMaxEvMix(50);
//      ana->SetNCentrBin(2);
//    }
//  }
//  
//  ana->SwitchOnFiducialCut();
//  
//  if(calorimeter=="EMCAL")
//  {  
//    // Avoid borders of EMCal, same as for isolation
//    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;
//    
//  }
//    
//  // Input / output delta AOD settings
//  
//  ana->SetInputAODName(Form("%sTrigger_%s",particle.Data(),kAnaGammaHadronCorr.Data()));
//  ana->SetAODObjArrayName(Form("%sHadronCorrIso%dTrigger_%s",particle.Data(),bIsolated,kAnaGammaHadronCorr.Data()));
//  ana->SetAODNamepTInConeHisto(Form("IC%s_%s_R%1.1f_ThMin%1.1f"           ,particle.Data(),kAnaGammaHadronCorr.Data(),cone,pth));
//  
//  //Set Histograms name tag, bins and ranges
//  
//  ana->AddToHistogramsName(Form("Ana%sHadronCorr_Iso%d_TM%d_",particle.Data(),bIsolated,tm));
//
//  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
//  
//  if(particle=="Hadron"  || particle.Contains("CTS"))
//  {
//    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
//    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
//  }
//  
//  return ana;
//}
//
/////
///// Configure the task doing the selected tracks checking
/////
//AliAnaChargedParticles* ConfigureChargedAnalysis( Bool_t simulation, Bool_t printSettings, Int_t   debug )
//{
//  AliAnaChargedParticles *ana = new AliAnaChargedParticles();
//  
//  // selection cuts
//  
//  ana->SetMinPt(0.2);
//  ana->SwitchOnFiducialCut();
//  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ; //more restrictive cut in reader and after in isolation
//  
//  ana->SwitchOffFillVertexBC0Histograms();
//
//  ana->SwitchOffFillPileUpHistograms();
//  ana->SwitchOffFillTrackBCHistograms();
//  
//  // Branch AOD settings
//  
//  ana->SetOutputAODName(Form("HadronTrigger_%s",kAnaGammaHadronCorr.Data()));
//  ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
//  
//  //Set Histograms name tag, bins and ranges
//  
//  ana->AddToHistogramsName("AnaHadrons_");
//
//  SetAnalysisCommonParameters(ana,"CTS",2012,"pp",simulation,printSettings,debug); // see method below
//  
//  return ana;
//}
//
/////
///// Configure the task doing standard calorimeter QA
/////
//AliAnaCalorimeterQA* ConfigureQAAnalysis(TString col,           Bool_t  simulation, 
//                                         TString calorimeter,   Int_t   year,     
//                                         Bool_t  printSettings, Int_t   debug      )
//{
//  AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
//
//  ana->SetCalorimeter(calorimeter);
//  
//  ana->SetTimeCut(-1e10,1e10); // Open time cut
//  
//  ana->SwitchOnCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
//   
//  ana->SwitchOnRealCaloAcceptance();
//  
//  ana->SwitchOffFiducialCut();
//  ana->SwitchOffFillAllTH3Histogram();
//  ana->SwitchOffFillAllPositionHistogram();
//  ana->SwitchOffFillAllPositionHistogram2();
//  ana->SwitchOffStudyBadClusters() ;
//  ana->SwitchOffStudyClustersAsymmetry();
//  ana->SwitchOffStudyWeight();
//  ana->SwitchOffFillAllCellTimeHisto() ;
//  
//  ana->SwitchOnFillAllTrackMatchingHistogram();
//
//  ana->AddToHistogramsName("QA_"); // Begining of histograms name
//
//  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
//  
//  return ana;
//}
//
/////
///// Configure the task filling generated particle kinematics histograms
/////
//AliAnaGeneratorKine* ConfigureGenKineAnalysis(Int_t   thresType,     Float_t cone,      
//                                              Float_t pth,  
//                                              TString col,           Bool_t  simulation, 
//                                              TString calorimeter,   Int_t   year,     
//                                              Bool_t  printSettings, Int_t   debug       )
//{  
//  AliAnaGeneratorKine *ana = new AliAnaGeneratorKine();
//  
//  // Trigger detector, acceptance and pT cut
//  ana->SetTriggerDetector(calorimeter);
//  ana->SetMinPt(2); // Trigger photon, pi0 minimum pT
//  ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.4, 90, 170);
//  
//  // Particles associated to trigger or isolation cone acceptance and pT cut
//  ana->SetCalorimeter(calorimeter);
//  ana->SetMinChargedPt(0.2);
//  ana->SetMinNeutralPt(0.3);
//  ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.65, 81, 179);
//  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360);
//  
//  // Isolation paramters
//  AliIsolationCut * ic = ana->GetIsolationCut();
//  ic->SetDebug(debug);
//  ic->SetPtThreshold(pth);
//  ic->SetConeSize(cone);
//  ic->SetSumPtThreshold(1.0) ;
//  ic->SetICMethod(thresType); 
//  
//  ana->AddToHistogramsName("AnaGenKine_");
//
//  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
//  
//  return ana;  
//}
//
//
/////
///// Set common histograms binning 
///// plus other analysis common settings like TRD covered super modules
///// the activation of the MC dedicated histograms and the activation of 
///// the debug mode
/////
//void SetAnalysisCommonParameters(AliAnaCaloTrackCorrBaseClass* ana, 
//                                 TString calorimeter,   Int_t  year, 
//                                 TString col,           Bool_t simulation,
//                                 Bool_t  printSettings, Int_t  debug)
//{
//  //
//  // Histograms ranges
//  //
//  AliHistogramRanges* histoRanges = ana->GetHistogramRanges();
//  
//  histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
//  
//  if(calorimeter=="EMCAL")
//  {
//    if ( year == 2010 )
//    {
//      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
//      histoRanges->SetHistoXRangeAndNBins(-230,90,120); // QA
//      histoRanges->SetHistoYRangeAndNBins(370,450,40);  // QA
//    }
//    else if ( year < 2014 )
//    {           
//      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
//      histoRanges->SetHistoXRangeAndNBins(-460,90,200); // QA
//      histoRanges->SetHistoYRangeAndNBins(100,450,100); // QA    }
//      else // Run2
//      {
//        histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 329*TMath::DegToRad(), 250) ;
//        histoRanges->SetHistoXRangeAndNBins(-460,460,230); // QA
//        histoRanges->SetHistoYRangeAndNBins(-450,450,225); // QA
//      }
//      
//      histoRanges->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
//    }
//    else if(calorimeter=="PHOS") 
//    {
//      histoRanges->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 320*TMath::DegToRad(), 70) ;
//      histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
//    }
//    else if(calorimeter=="CTS")
//    {
//      ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
//      ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
//    }
//  }
//  
//  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 4.9, 500);
//  
//  // Invariant mass histo
//  histoRanges->SetHistoMassRangeAndNBins(0., 1., 200) ;
//  histoRanges->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
//  
//  // check if time calibration is on
//  //histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
//  histoRanges->SetHistoTimeRangeAndNBins(-400.,400,400);
//  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
//  
//  // track-cluster residuals
//  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.15,0.15,300);
//  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.15,0.15,300);
//  histoRanges->SetHistodRRangeAndNBins(0.,0.15,150);//QA
//
//  // QA, electron, charged
//  histoRanges->SetHistoPOverERangeAndNBins(0,2.,200);
//  histoRanges->SetHistodEdxRangeAndNBins(0.,200.,200);
//  
//  // QA
//  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
//  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
//  histoRanges->SetHistoZRangeAndNBins(-350,350,175);
//  histoRanges->SetHistoRRangeAndNBins(430,460,30);
//  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
//  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
//  
//  // QA, correlation
//  if(col=="PbPb")
//  {
//    histoRanges->SetHistoNClusterCellRangeAndNBins(0,100,100);
//    histoRanges->SetHistoNClustersRangeAndNBins(0,500,50);
//    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,2000,200);
//  }
//  else
//  {
//    histoRanges->SetHistoNClusterCellRangeAndNBins(0,50,50);
//    histoRanges->SetHistoNClustersRangeAndNBins(0,50,50);
//    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,200,200);
//  }
//  
//  // xE, zT
//  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,200);
//  histoRanges->SetHistoHBPRangeAndNBins  (0.,10.,200);
//  
//  // Isolation
//  histoRanges->SetHistoPtInConeRangeAndNBins(0, 50 , 250);
//  histoRanges->SetHistoPtSumRangeAndNBins   (0, 100, 250);
//  
//  //
//  // TRD SM
//  //
//  if     (year == 2011) ana->SetFirstSMCoveredByTRD( 6);
//  else if(year == 2012 || 
//          year == 2013) ana->SetFirstSMCoveredByTRD( 4);
//  else                  ana->SetFirstSMCoveredByTRD(-1);
//
//  //
//  // MC histograms?
//  //
//  if(simulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
//  else           ana->SwitchOffDataMC() ;
//  
//  //Set here generator name, default pythia
//  //ana->GetMCAnalysisUtils()->SetMCGenerator("");
//
//  //
//  // Debug
//  //
//  if(printSettings) ana->Print("");
//  
//  ana->SetDebug(debug); // 10 for lots of messages
//}
//
/////
///// Set the trigger requested for the analysis, 
///// depending on a string given
/////
//UInt_t SetTriggerMaskFromName(TString trigger)
//{
//  if(trigger=="EMC7")
//  {
//    printf("CaloTrackCorr trigger EMC7\n");
//    return AliVEvent::kEMC7;
//  }
//  else if (trigger=="INT7")
//  {
//    printf("CaloTrackCorr trigger INT7\n");
//    return AliVEvent::kINT7;
//  }
//  else if(trigger=="EMC1")
//  {
//    printf("CaloTrackCorr trigger EMC1\n");
//    return AliVEvent::kEMC1;
//  }
//  else if(trigger=="MB")
//  {
//    printf("CaloTrackCorr trigger MB\n");
//    return AliVEvent::kMB;
//  }  
//  else if(trigger=="PHOS")
//  {
//    printf("CaloTrackCorr trigger PHOS\n");
//    return AliVEvent::kPHI7;
//  }  
//  else if(trigger=="PHOSPb")
//  {
//    printf("CaloTrackCorr trigger PHOSPb\n");
//    return AliVEvent::kPHOSPb;
//  }
//  else if(trigger=="AnyINT")
//  {
//    printf("CaloTrackCorr trigger AnyINT\n");
//    return AliVEvent::kAnyINT;
//  }  
//  else if(trigger=="INT")
//  {
//    printf("CaloTrackCorr trigger AnyINT\n");
//    return AliVEvent::kAny;
//  }
//  else if(trigger=="EMCEGA")
//  {
//    printf("CaloTrackCorr trigger EMC Gamma\n");
//    return AliVEvent::kEMCEGA;
//  } 
//  else if(trigger=="EMCEJE")
//  {
//    printf("CaloTrackCorr trigger EMC Jet\n");
//    return AliVEvent::kEMCEJE;
//  }
//  else if(trigger=="Central")
//  {
//    printf("CaloTrackCorr trigger Central\n");
//    return AliVEvent::kCentral;
//  }
//  else if(trigger=="CentralEGA")
//  {
//    printf("CaloTrackCorr trigger Central+EMCEGA\n");
//    return (AliVEvent::kCentral | AliVEvent::kEMCEGA);
//  }
//  else if(trigger=="SemiCentral")
//  {
//    printf("CaloTrackCorr trigger SemiCentral\n");
//    return AliVEvent::kSemiCentral;
//  }
//  else if(trigger=="SemiOrCentral")
//  {
//    printf("CaloTrackCorr trigger SemiCentral Or Central\n");
//    return (AliVEvent::kSemiCentral | AliVEvent::kCentral);
//  }  
//  else return AliVEvent::kAny;
//}
//





Bool_t  kPrint         = kFALSE;
Bool_t  kSimulation    = kFALSE;
Bool_t  kUseKinematics = kFALSE;
Bool_t  kOutputAOD     = kFALSE;
Bool_t  kEventSelection= kFALSE;
Bool_t  kExotic        = kTRUE;
Bool_t  kNonLinearity  = kFALSE;
Int_t   kYears         = 2011;
TString kCollisions    = "pp";
TString kTrig          = "EMC7" ;
TString kClusterArray  = "";
TString kData          = ""; // MC or deltaAOD
TString kInputDataType = "ESD";
TString kCalorimeter   = "EMCAL";
Bool_t  kTM            = kTRUE;
Bool_t  kRecalTM       = kTRUE;
Int_t   kMinCen        = -1;
Int_t   kMaxCen        = -1;
TString kName          = "";
Int_t   kDebug         = -1;
Bool_t  kQA            = kFALSE;
Bool_t  kHadronAN      = kFALSE;
Bool_t  kCalibE        = kTRUE;
Bool_t  kCalibT        = kTRUE;
Bool_t  kBadMap        = kTRUE;
Bool_t  kTender        = kFALSE;
Bool_t  kMix           = kFALSE;
Int_t   kRunNumber     = -1;

AliAnalysisTaskCaloTrackCorrelation *AddTaskGammaHadronCorrelation(const TString  data          = "",
                                                          const TString  calorimeter   = "EMCAL",
                                                          const Bool_t   simulation    = kFALSE,
                                                          const Bool_t   eventsel      = kFALSE,
                                                          const Bool_t   exotic        = kTRUE,
                                                          const Bool_t   nonlin        = kFALSE,
                                                          TString        outputfile    = "",
                                                          const Int_t    year          = 2010,
                                                          const TString  col           = "pp",
                                                          const TString  trigger       = "MB",
                                                          const TString  clustersArray = "V1",
                                                          const Bool_t   mix           = kTRUE,
                                                          const Bool_t   recaltm       = kTRUE,
                                                          const Bool_t   tm            = kTRUE,
                                                          const Int_t    minCen        = -1,
                                                          const Int_t    maxCen        = -1,
                                                          const Bool_t   qaan          = kFALSE,
                                                          const Bool_t   hadronan      = kFALSE,
                                                          const Bool_t   calibE        = kTRUE,
                                                          const Bool_t   badmap        = kTRUE,
                                                          const Bool_t   calibT        = kTRUE,
                                                          const Bool_t   tender        = kFALSE,
                                                          const Bool_t   outputAOD     = kFALSE,
                                                          const Bool_t   printSettings = kFALSE,
                                                          const Double_t scaleFactor   = -1,
                                                          const Int_t    runNumber     = -1
                                                          )
{
    // Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
    
    kPrint         = printSettings;
    kSimulation    = simulation;
    kYears         = year;
    kCollisions    = col;
    kExotic        = exotic;
    kNonLinearity  = nonlin;
    kTrig          = trigger;
    kClusterArray  = clustersArray;
    kData          = data;
    kCalorimeter   = calorimeter;
    kOutputAOD     = outputAOD;
    kTM            = tm;
    kRecalTM       = recaltm;
    kMinCen        = minCen;
    kMaxCen        = maxCen;
    kEventSelection= eventsel;
    kQA            = qaan;
    kHadronAN      = hadronan;
    kCalibE        = calibE;
    kCalibT        = calibT;
    kBadMap        = badmap;
    kTender        = tender;
    kMix           = mix;
    kRunNumber     = runNumber;
    
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
    
    ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);
    
    kInputDataType = "AOD";
    if(!kData.Contains("delta"))
        kInputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    if(kSimulation)
    {
        kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
        if (!kUseKinematics && data=="AOD" && kInputDataType != "ESD") kUseKinematics = kTRUE; //AOD primary should be available ...
    }
    
    cout<<"********* ACCESS KINE? "<<kUseKinematics<<endl;
    
    // Name for containers
    
    //kName = Form("%s_Trig%s_Cl%s_TM%d_R%1.1f_Pt%1.1f",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data(),kTM,cone,pth);
    kName = Form("%s_Trig%s_Cl%s_TM%d",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data(),kTM);
    
    
    if(kCollisions=="PbPb" && kMaxCen>=0) kName+=Form("Cen%d_%d",kMinCen,kMaxCen);
    
    printf("<<<< NAME: %s >>>>>\n",kName.Data());
    
    // #### Configure analysis ####
    
    AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
    //printf("SCALE FACTOR %e\n",scaleFactor);
    maker->SetScaleFactor(scaleFactor); // for MC, negative (not scaled) by default
    
    // General frame setting and configuration
    maker->SetReader   (ConfigureReader()   );
    maker->SetCaloUtils(ConfigureCaloUtils());
    
    // Analysis tasks setting and configuration
    Int_t n = 0;//Analysis number, order is important
    
    // Isolation settings
    Int_t partInCone = AliIsolationCut::kNeutralAndCharged; // kOnlyCharged;
    Int_t thresType  = AliIsolationCut::kPtThresIC;
    //Int_t thresType  = AliIsolationCut::kSumPtIC;
    Float_t cone   = 0.4;
    Float_t pth    = 0.5;
    Float_t pthMax = 1000;
    
    //  maker->AddAnalysis(ConfigureQAAnalysis(), n++);
    
    // Photon analysis
    maker->AddAnalysis(ConfigurePhotonAnalysis(), n++); // Photon cluster selection
    maker->AddAnalysis(ConfigurePi0Analysis(), n++);
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0"        , AliAnaPi0EbE::kIMCalo,kFALSE), n++); // Pi0 event by event selection, invariant mass and photon tagging from decay
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Eta"        , AliAnaPi0EbE::kIMCalo,kFALSE), n++); // Eta event by event selection, invariant mass and photon tagging from decay
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0SideBand", AliAnaPi0EbE::kIMCalo,kFALSE), n++); // Pi0 out of peak event by event selection, and photon tagging from decay
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("EtaSideBand", AliAnaPi0EbE::kIMCalo,kFALSE), n++); // Eta out of peak event by event selection, and photon tagging from decay
    printf("\n\nBEFORE ISOLATION PHOTON\n\n");
    maker->AddAnalysis(ConfigureIsolationAnalysis("Photon", partInCone,thresType, cone, pth), n++); // Photon isolation
    printf("\n\nBEFORE CORRELATION PHOTON NON ISO\n\n");
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kFALSE), n++); // Gamma hadron correlation
    printf("\n\nBEFORE CORRELATION PHOTON ISO\n\n");
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kTRUE, partInCone,thresType, cone, pth) , n++); // Isolated gamma hadron correlation
    printf("");
    
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0"        , AliAnaPi0EbE::kIMCalo,kTRUE), n++); // Pi0 event by event+isolation selection, invariant mass and photon tagging from decay
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Eta"        , AliAnaPi0EbE::kIMCalo,kTRUE), n++); // Eta event by event+isolation selection, invariant mass and photon tagging from decay
    //
    ////  // Merged pi0 analysis
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kSSCalo,kTRUE,kTRUE), n++); // Pi0 event by event selection, cluster splitting
    printf("\n\nBEFORE ISOLATION PI0\n\n");
    maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SS", partInCone,thresType, cone, pth), n++);       // Pi0 isolation, cluster splits
    printf("\n\nBEFORE CORRELATION PI0 NON ISO\n\n");
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS" ,kFALSE), n++); // Pi0 hadron correlation
    printf("\n\nBEFORE CORRELATION PHOTON ISO\n\n");
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS" ,kTRUE, partInCone,thresType, cone, pth) , n++); // Isolated pi0 hadron correlation
    
    
    // Charged analysis
    maker->AddAnalysis(ConfigureChargedAnalysis(), n++); // track selection
    
    maker->SetAnaDebug(kDebug)  ;
    maker->SwitchOnHistogramsMaker()  ;
    if(kData.Contains("delta")) maker->SwitchOffAODsMaker() ;
    else                        maker->SwitchOnAODsMaker()  ;
    
    //if(kSimulation || !kTrig.Contains("EMC"))
    maker->SwitchOffDataControlHistograms();
    
    if(kPrint) maker->Print("");
    
    printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, kCalorimeter.Data());
    
    // Create task
    
    AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("CaloTrackCorr%s",kName.Data()));
    task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
    task->SetDebugLevel(kDebug);
    //task->SetDebugLevel(100);
    //task->SetBranches("ESD:AliESDRun.,AliESDHeader");
    task->SetAnalysisMaker(maker);
    mgr->AddTask(task);
    
    //Create containers
    
    if(outputfile.Length()==0) outputfile = AliAnalysisManager::GetCommonFileName();
    
    AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(kName, TList::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               Form("%s",outputfile.Data()));
    
    AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",kName.Data()), TList::Class(),
                                                               AliAnalysisManager::kParamContainer,
                                                               "AnalysisParameters.root");
    
    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //==============================================================================
    mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
    // AOD output slot will be used in a different way in future
    if(!kData.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
    mgr->ConnectOutput (task, 1, cout_pc);
    mgr->ConnectOutput (task, 2, cout_cuts);
    
    if(!kMix)
    {
        UInt_t mask =  SetTriggerMaskFromName();
        task->SelectCollisionCandidates(mask);
    }
    
    AliAnaCaloTrackCorrMaker* maker = task->GetAnalysisMaker();
    AliAnaParticleHadronCorrelation* corrIso0 = (AliAnaParticleHadronCorrelation*)maker->GetListOfAnalysisContainers().At(6);
    AliAnaParticleHadronCorrelation* corrIso1 = (AliAnaParticleHadronCorrelation*)maker->GetListOfAnalysisContainers().At(7);
    corrIso0->SwitchOnBackgroundBinsPtInConeHistograms();
    corrIso1->SwitchOnBackgroundBinsPtInConeHistograms();
    
    return task;
}

//____________________________________
AliCaloTrackReader * ConfigureReader()
{
    
    AliCaloTrackReader * reader = 0;
    if     (kInputDataType == "ESD"&& kData=="MC" )
        reader = new AliCaloTrackMCReader();
    else if(kInputDataType=="AOD" || kData.Contains("AOD"))
        reader = new AliCaloTrackAODReader();
    else if(kInputDataType=="ESD")
        reader = new AliCaloTrackESDReader();
    else
        printf("AliCaloTrackReader::ConfigureReader() - Data combination not known kData=%s, kInputData=%s\n",kData.Data(),kInputDataType.Data());
    
    reader->SetDebug(kDebug);//10 for lots of messages
    //reader->SetDebug(100);
    
    // Event triggered selection settings
    if(!kSimulation && (kTrig.Contains("EMC") || kTrig.Contains("L")))
    {
        printf("=== Remove bad triggers === \n");
        reader->SwitchOnTriggerPatchMatching();
        reader->SwitchOnBadTriggerEventsRemoval();
        
        reader->SetTriggerPatchTimeWindow(8,9);
        if     (kRunNumber < 146861) reader->SetEventTriggerL0Threshold(3.);
        else if(kRunNumber < 154000) reader->SetEventTriggerL0Threshold(4.);
        else if(kRunNumber < 165000) reader->SetEventTriggerL0Threshold(5.5);
        //redefine for other periods, triggers
        
        if(kRunNumber < 172000)
        {
            reader->SetEventTriggerL1Bit(4,5); // current LHC11 data
            printf("\t Old L1 Trigger data format!\n");
        }
        else
        {
            reader->SetEventTriggerL1Bit(6,8); // LHC12-13 data
            printf("\t Current L1 Trigger data format!\n");
        }
        
        if(kClusterArray != "" || kTender)
        {
            printf("Trigger cluster calibration OFF\n");
            reader->SwitchOffTriggerClusterTimeRecal() ;
        }
        
    }
    else
    {
        reader->SwitchOffTriggerPatchMatching();
        reader->SwitchOffBadTriggerEventsRemoval();
    }
    // In case of Pythia pt Hard bin simulations (jet-jet, gamma-jet)
    // reject some special events that bother the cross section
    if(kSimulation)
    {
        // Event rejection cuts for jet-jet simulations, do not use in other
        reader->SetPtHardAndJetPtComparison(kTRUE);
        reader->SetPtHardAndJetPtFactor(4);
        
        reader->SetPtHardAndClusterPtComparison(kTRUE);
        reader->SetPtHardAndClusterPtFactor(1.5);
    }
    
    //Delta AOD?
    //reader->SetDeltaAODFileName("");
    if(kOutputAOD) reader->SwitchOnWriteDeltaAOD()  ;
    
    // MC settings
    if(kUseKinematics){
        if(kInputDataType == "ESD"){
            reader->SwitchOnStack();
            reader->SwitchOffAODMCParticles();
        }
        else if(kInputDataType == "AOD"){
            reader->SwitchOffStack();
            reader->SwitchOnAODMCParticles();
        }
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
    if(kSimulation)
    {
        reader->SwitchOffShowerShapeSmearing(); // Active only on MC
        //reader->SetShowerShapeSmearWidth(0.002);
        //reader->SetShowerShapeSmearWidth(0.003);
        reader->SetShowerShapeSmearWidth(0.005);
    }
    
    // Time cuts
    if(kSimulation)
    {
        reader->SwitchOffUseTrackTimeCut();
        reader->SwitchOffUseParametrizedTimeCut();
        reader->SwitchOffUseEMCALTimeCut();
        reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
    }
    else
    {
        if(kCalibT)
        {
            printf("Set time cut parameters for run %d\n",kRunNumber);
            reader->SwitchOnUseEMCALTimeCut();
            reader->SwitchOnUseParametrizedTimeCut();
            
            //Absolute window
            reader->SetEMCALTimeCut(-25,20);
            
            //Parametrization
            if     (kRunNumber >= 151636 && kRunNumber <= 155384 )
            {
                printf("Set time parameters for LHC11c\n");
                reader->SetEMCALParametrizedMinTimeCut(0,-5  ); reader->SetEMCALParametrizedMinTimeCut(1,-1 ); reader->SetEMCALParametrizedMinTimeCut(2, 1.87); reader->SetEMCALParametrizedMinTimeCut(3, 0.4);
                reader->SetEMCALParametrizedMaxTimeCut(0, 3.5); reader->SetEMCALParametrizedMaxTimeCut(1, 50); reader->SetEMCALParametrizedMaxTimeCut(2, 0.15); reader->SetEMCALParametrizedMaxTimeCut(3, 1.6);
            }
            else if(kRunNumber >= 156447 && kRunNumber <= 159635 )
            {
                printf("Set time parameters for LHC11d\n");
                reader->SetEMCALParametrizedMinTimeCut(0,-5);  reader->SetEMCALParametrizedMinTimeCut(1,-1 );  reader->SetEMCALParametrizedMinTimeCut(2, 3.5 ); reader->SetEMCALParametrizedMinTimeCut(3, 1.  );
                reader->SetEMCALParametrizedMaxTimeCut(0, 5);  reader->SetEMCALParametrizedMaxTimeCut(1, 50);  reader->SetEMCALParametrizedMaxTimeCut(2, 0.45); reader->SetEMCALParametrizedMaxTimeCut(3, 1.25);
            }
            else
            {
                reader->SwitchOffUseParametrizedTimeCut();
            }
        }
        else
        {
            reader->SwitchOffUseParametrizedTimeCut();
            reader->SwitchOffUseEMCALTimeCut();
            reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
        }
    }
    
    reader->SwitchOffUseTrackTimeCut();
    reader->SetTrackTimeCut(0,50);
    
    reader->SwitchOnFiducialCut();
    reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;
    
    reader->SwitchOffUseTrackDCACut();
    //reader->SetTrackDCACut(0,0.0105);
    //reader->SetTrackDCACut(1,0.035);
    //reader->SetTrackDCACut(2,1.1);
    
    // Tracks
    reader->SwitchOnCTS();
    
    if(kInputDataType=="ESD")
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
    else if(kInputDataType=="AOD")
    {
        reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
        //reader->SwitchOffAODTrackSharedClusterSelection();
        //reader->SetTrackStatus(AliVTrack::kITSrefit);
        //reader->SwitchOnTrackHitSPDSelection();    // Check that the track has at least a hit on the SPD, not much sense to use for hybrid or TPC only tracks
        //reader->SetTrackFilterMask(128);           // Filter bit, not mask, use if off hybrid, TPC only
    }
    
    // Calorimeter
    
    reader->SetEMCALClusterListName(kClusterArray);
    if(kClusterArray == "" && !kTender)
    {
        printf("**************** Standard EMCAL clusters branch analysis **************** \n");
        reader->SwitchOnClusterRecalculation();
        // Check in ConfigureCaloUtils that the recalibration and bad map are ON
    }
    else
    {
        printf("**************** Input for analysis is Clusterizer %s **************** \n", kClusterArray.Data());
        reader->SwitchOffClusterRecalculation();
    }
    
    // CAREFUL
    if(!kNonLinearity)reader->SwitchOffClusterELinearityCorrection();
    
    //  if(kQA) reader->SwitchOnClusterRecalculation();
    //  else    reader->SwitchOffClusterRecalculation();
    
    //if(kCalorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();
    reader->SwitchOnEMCAL();
    //}
    //if(kCalorimeter == "PHOS") { // Should be on if QA is activated with correlation on
    reader->SwitchOffPHOSCells();
    reader->SwitchOffPHOS();
    //}
    
    // for case data="deltaAOD", no need to fill the EMCAL/PHOS cluster lists
    if(kData.Contains("delta"))
    {
        reader->SwitchOffEMCAL();
        reader->SwitchOffPHOS();
        reader->SwitchOffEMCALCells();
        reader->SwitchOffPHOSCells();
    }
    
    //-----------------
    // Event selection
    //-----------------
    
    //if(!kUseKinematics) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
    
    reader->RejectFastClusterEvents() ;
    
    // For mixing with AliAnaParticleHadronCorrelation switch it off
    if(kMix)
    {
        reader->SwitchOffEventTriggerAtSE();
        UInt_t mask =  SetTriggerMaskFromName();
        reader->SetEventTriggerMask(mask); // Only for mixing and SwitchOffEventTriggerAtSE();
        //reader->SetMixEventTriggerMask(AliVEvent::kMB); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
        reader->SetMixEventTriggerMask(AliVEvent::kINT7); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
        
        printf("---Trigger selection done in AliCaloTrackReader!!!\n");
    }
    else
        reader->SwitchOnEventTriggerAtSE();
    
    reader->SetZvertexCut(10.);                // Open cut
    reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
    reader->SwitchOnRejectNoTrackEvents();
    
    reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
    reader->SwitchOffPileUpEventRejection();         // remove pileup by default
    
    if(kCollisions=="PbPb")
    {
        // Centrality
        reader->SetCentralityClass("V0M");
        reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
        reader->SetCentralityBin(kMinCen,kMaxCen); // Accept all events, if not select range
        
        // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
        reader->SetEventPlaneMethod("V0");
    }
    
    if(kPrint) reader->Print("");
    
    return reader;
    
}

//_______________________________________
AliCalorimeterUtils* ConfigureCaloUtils()
{
    
    AliCalorimeterUtils *cu = new AliCalorimeterUtils;
    cu->SetDebug(kDebug);
    //cu->SetDebug(100);
    
    // Remove clusters close to borders, at least max energy cell is 1 cell away
    cu->SetNumberOfCellsFromEMCALBorder(1);
    cu->SetNumberOfCellsFromPHOSBorder(2);
    
    cu->SetNumberOfSuperModulesUsed(10);
    
    // Search of local maxima in cluster
    if(kCollisions=="pp")
    {
        cu->SetLocalMaximaCutE(0.1);
        cu->SetLocalMaximaCutEDiff(0.03);
        
        if(kName.Contains("150"))
        {
            printf("Reclusterize with 150 threshold, set PbPb settings\n");
            cu->SetLocalMaximaCutE(0.2);
            cu->SetLocalMaximaCutEDiff(0.03);
        }
    }
    else
    {
        cu->SetLocalMaximaCutE(0.2);
        cu->SetLocalMaximaCutEDiff(0.03);
    }
    
    cu->SwitchOffClusterPlot();
    
    if(kRecalTM) cu->SwitchOnRecalculateClusterTrackMatching(); // Done in clusterization
    else         cu->SwitchOffRecalculateClusterTrackMatching();
    
    cu->SwitchOnBadChannelsRemoval() ;
    
    //EMCAL settings
    
    if(!kSimulation)
        cu->SwitchOnLoadOwnEMCALGeometryMatrices();
    
    AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
    
    if(kClusterArray == "" && !kTender)
    {
        cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
        cu->SwitchOnRunDepCorrection();
    }
    
    if(!kSimulation)
    {
        cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
        cu->SwitchOnRunDepCorrection();
    }
    else
    {
        cu->SwitchOffRecalibration(); // Check the reader if it is taken into account during filtering
        cu->SwitchOffRunDepCorrection();
    }
    
    //gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigureEMCALRecoUtils.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
    ConfigureEMCALRecoUtils(recou,
                            kSimulation,
                            kExotic,
                            kNonLinearity,
                            kCalibE,
                            kBadMap,
                            kCalibT);
    recou->SetExoticCellDiffTimeCut(50);
    
    if( kNonLinearity )
    {
        //    printf("ConfigureCaloUtils() - Apply non linearity to EMCAL\n");
        //    //CAREFUL only for the latest simulation
        //    recou->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrected);
        //    recou->SetNonLinearityParam(0,9.81039e-01);
        //    recou->SetNonLinearityParam(1,1.13508e-01);
        //    recou->SetNonLinearityParam(2,1.00173e+00);
        //    recou->SetNonLinearityParam(3,9.67998e-02);
        //    recou->SetNonLinearityParam(4,2.19381e+02);
        //    recou->SetNonLinearityParam(5,6.31604e+01);
        //    recou->SetNonLinearityParam(6,1);
        if(kRunNumber <= 150000 && kSimulation) {// CAREFUL only for old productions
            printf("\t ------>>>>> Set NonLin kPi0MC\n");
            recou->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MC);
        }
        cu->SwitchOnCorrectClusterLinearity();
    }
    
    
    //  if(!kQA)
    //  {
    //    //cu->SwitchOffEMCALOADB() ;
    //    cu->SwitchOffRecalibration();
    //    cu->SwitchOffRunDepCorrection();
    //    cu->SwitchOffTimeRecalibration();
    //  }
    
    printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
    printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
    
    
    // PHOS
    cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
    if(kPrint) cu->Print("");
    
    return cu;
    
}

//_____________________________________
AliAnaPhoton* ConfigurePhotonAnalysis()
{
    
    AliAnaPhoton *ana = new AliAnaPhoton();
    ana->SetDebug(kDebug); //10 for lots of messages
    
    // cluster selection cuts
    
    ana->SwitchOnRealCaloAcceptance();
    
    ana->SwitchOffFiducialCut();
    
    ana->SetCalorimeter(kCalorimeter);
    
    ana->SetFirstSMCoveredByTRD(-1);
    
    ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection
    
    //if(!kSimulation) ana->SwitchOnFillPileUpHistograms();
    
    if(kTM) ana->SwitchOnTrackMatchRejection() ;
    else    ana->SwitchOffTrackMatchRejection() ;
    ana->SwitchOffTMHistoFill() ;
    
    if(kCalorimeter == "PHOS")
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
    
    //EMCAL
    //caloPID->SetEMCALLambda0CutMax(0.27);
    caloPID->SetEMCALLambda0CutMax(10);
    caloPID->SetEMCALLambda0CutMin(0.10);
    
    caloPID->SetEMCALDEtaCut(0.025);
    caloPID->SetEMCALDPhiCut(0.030);
    
    //PHOS
    caloPID->SetPHOSDispersionCut(2.5);
    caloPID->SetPHOSRCut(2.);
    if(kInputData=="AOD") caloPID->SetPHOSRCut(2000.); // Open cut since dX, dZ not stored
    
    // Input / output delta AOD settings
    
    if(!kData.Contains("delta"))
    {
        ana->SetOutputAODName(Form("Photon%s",kName.Data()));
        //ana->SetOutputAODClassName("AliAODCaloTrackParticleCorrelation");
        ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
        //ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
    }
    else ana->SetInputAODName(Form("Photon%s",kName.Data()));
    
    //Set Histograms name tag, bins and ranges
    
    ana->AddToHistogramsName(Form("AnaPhoton_TM%d_",kTM));
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    // Number of particle type MC histograms
    ana->FillNOriginHistograms (4); // 14 max
    ana->FillNPrimaryHistograms(4); // 6 max
    
    ConfigureMC(ana);
    
    if(kPrint) ana->Print("");
    
    return ana;
    
}

//_______________________________
AliAnaPi0* ConfigurePi0Analysis()
{
    
    AliAnaPi0 *ana = new AliAnaPi0();
    
    ana->SetDebug(kDebug);//10 for lots of messages
    
    // Input delta AOD settings
    ana->SetInputAODName(Form("Photon%s",kName.Data()));
    
    ana->SwitchOnRealCaloAcceptance();
    
    //ana->SwitchOnFiducialCut();
    //ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.65, 81, 179);
    
    // Calorimeter settings
    ana->SetCalorimeter(kCalorimeter);
    
    //settings for pp collision mixing
    ana->SwitchOffOwnMix(); //Off when mixing done with general mixing frame
    
    // Cuts
    if(kCalorimeter=="EMCAL") ana->SetPairTimeCut(70);
    
    ana->SetNAsymCuts(1); // no asymmetry cut, previous studies showed small effect.
    //In EMCAL assymetry cut prevents combination of assymetric decays which is the main source of pi0 at high E.
    
    if     (kCollisions=="pp"  )
    {
        ana->SetNCentrBin(1);
        ana->SetNZvertBin(10);
        ana->SetNRPBin(1);
        ana->SetNMaxEvMix(100);
        ana->SwitchOnSMCombinations();
    }
    else if(kCollisions=="PbPb")
    {
        ana->SetNCentrBin(10);
        ana->SetNZvertBin(10);
        ana->SetNRPBin(4);
        ana->SetNMaxEvMix(10);
        ana->SwitchOffSMCombinations();
    }
    
    ana->SwitchOffMultipleCutAnalysis();
    
    //Set Histograms name tag, bins and ranges
    
    ana->AddToHistogramsName(Form("AnaPi0_TM%d_",kTM));
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    ana->SwitchOnFillOriginHisto();
    ConfigureMC(ana);
    
    if(kPrint) ana->Print("");
    
    return ana;
    
}

//_____________________________________________________
AliAnaPi0EbE* ConfigurePi0EbEAnalysis(TString particle,
                                      Int_t analysis, Bool_t useSSIso = kTRUE, Bool_t useAsy = kTRUE)
{
    
    AliAnaPi0EbE *ana = new AliAnaPi0EbE();
    ana->SetDebug(kDebug);//10 for lots of messages
    
    ana->SetAnalysisType(analysis);
    TString opt = "";
    if(analysis==AliAnaPi0EbE::kIMCaloTracks) opt = "Conv";
    if(analysis==AliAnaPi0EbE::kSSCalo)       opt = "SS";
    
    ana->SwitchOffAllNLMHistoFill();
    ana->SetFirstSMCoveredByTRD(-1);
    ana->SwitchOffSelectedClusterHistoFill();
    
    //  if(opt.Contains("SS"))  ana->SwitchOnSelectedClusterHistoFill();
    //  else                    ana->SwitchOffSelectedClusterHistoFill();
    
    ana->SwitchOffFillWeightHistograms();
    //if(!kSimulation) ana->SwitchOnFillPileUpHistograms();
    //if(!kTime && !kSimulation) ana->SwitchOnFillEMCALBCHistograms();
    
    if(kTM) ana->SwitchOnTrackMatchRejection() ;
    else    ana->SwitchOffTrackMatchRejection() ;
    ana->SwitchOffTMHistoFill() ;
    
    ana->SetCalorimeter(kCalorimeter);
    
    // Output delta AOD settings
    
    if(!kInputDataType.Contains("delta"))
    {
        ana->SetOutputAODName(Form("%s%s%s",particle.Data(), opt.Data(), kName.Data()));
        printf("***Out branch %s***\n",ana->GetOutputAODName().Data());
        //ana->SetOutputAODClassName("AliAODCaloTrackParticleCorrelation");
        ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    }
    else
        ana->SetInputAODName(Form("%s%s%s",particle.Data(),opt.Data(),kName.Data()));
    
    if(analysis == AliAnaPi0EbE::kIMCaloTracks) ana->SetInputAODGammaConvName("PhotonsCTS");
    
    //Set Histograms name tag, bins and ranges
    
    ana->AddToHistogramsName(Form("Ana%s%sEbE_TM%d_",particle.Data(),opt.Data(),kTM));
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    if(kPrint) ana->Print("");
    
    ConfigureMC(ana);
    
    ///////////////////////////////////
    if(analysis!=AliAnaPi0EbE::kSSCalo)
    {
        ana->SetInputAODName(Form("Photon%s",kName.Data()));
        
        ana->SwitchOnSelectPairInIsolationCone();
        ana->SetR(0.4);
        ana->SetIsolationCandidateMinPt(5);
        
        if(useSSIso)
        {
            ana->SwitchOnSelectIsolatedDecay();
            ana->AddToHistogramsName(Form("Ana%s%sEbEIsoDecay_TM%d_",particle.Data(),opt.Data(),kTM));
            ana->SetOutputAODName(Form("%s%sIsoDecay%s",particle.Data(), opt.Data(), kName.Data()));
        }
        
        if(kCalorimeter=="EMCAL" && !kSimulation) ana->SetPairTimeCut(100);
        
        AliNeutralMesonSelection *nms = ana->GetNeutralMesonSelection();
        nms->SetParticle(particle);
        
        // Tighten a bit mass cut with respect to default window
        if(particle=="Pi0") nms->SetInvMassCutRange(0.110,0.160);
        if(particle=="Eta") nms->SetInvMassCutRange(0.520,0.580);
        
        //if(!particle.Contains("SideBand")) nms->SwitchOnAngleSelection();
        //else nms->SwitchOnAngleSelection();
        
        nms->SwitchOffAngleSelection();
        if(particle.Contains("Pi0SideBand")) // For pi0, do not consider left band
            nms->SetSideBandCutRanges(-1,0,0.170,0.210);
        
        if(particle.Contains("EtaSideBand")) // For pi0, do not consider left band
            nms->SetSideBandCutRanges(0.450,0.510,0.590,0.650);
        
        nms->KeepNeutralMesonSelectionHistos(kTRUE);
        //nms->SetAngleMaxParam(2,0.2);
        nms->SetHistoERangeAndNBins(0, 20, 80) ;
        //nms->SetHistoIMRangeAndNBins(0, 1, 400);
    }
    else  ///////////////////////////////////
    {
        // cluster splitting settings
        ana->SetMinEnergy(6);
        ana->SetMaxEnergy(200.);
        
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
        
        caloPID->SetSplitWidthSigma(3); // cut at 3 sigma of the mean pi0 peak.
        
        if(!useSSIso)
        {
            printf("Do not apply SS cut on merged pi0 analysis \n");
            caloPID->SwitchOffSplitShowerShapeCut() ;
            ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_TM%d_",particle.Data(),opt.Data(),kTM));
            ana->SetOutputAODName(Form("%s%s%s_OpenSS",particle.Data(), opt.Data(), kName.Data()));
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
                ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_OpenAsy_TM%d_",particle.Data(),opt.Data(),kTM));
                ana->SetOutputAODName(Form("%s%s%s_OpenSS_OpenAsy",particle.Data(), opt.Data(), kName.Data()));
            }
            else
            {
                ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenAsy_TM%d_",particle.Data(),opt.Data(),kTM));
                ana->SetOutputAODName(Form("%s%s%s_OpenAsy",particle.Data(), opt.Data(), kName.Data()));
            }
        }
        
        //For Pi0 only if  SwitchOnSimpleSplitMassCut()
        caloPID->SetPi0MassRange(0.10, 0.18);
        caloPID->SetEtaMassRange(0.40, 0.60);
        caloPID->SetPhotonMassRange(0.00, 0.08);
        
        caloPID->SetClusterSplittingMinNCells(6);
        
        //caloPID->SetSplitEnergyFractionMinimum(0, 0.95);
        //caloPID->SetSplitEnergyFractionMinimum(1, 0.95);
        //caloPID->SetSplitEnergyFractionMinimum(2, 0.8);
        
        if(kCollisions=="PbPb" || kName.Contains("150"))
        {
            caloPID->SetClusterSplittingMinNCells(4);
            //caloPID->SetPi0MassShiftHighECell(0.005);
        }
    }
    ///////////////////////////////////
    
    return  ana;
    
}

//____________________________________________________________________________________________________
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle="Photon",
                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                    Float_t cone = 0.4,
                                                    Float_t pth  = 0.5, Float_t pthMax = -1,
                                                    Bool_t multi = kFALSE)
{
    
    AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
    ana->SetDebug(kDebug);
    
    ana->SetMinPt(5);
    ana->SetCalorimeter(kCalorimeter);
    
    ana->SwitchOffUEBandSubtractionHistoFill();
    ana->SwitchOffCellHistoFill() ;
    
    ana->SwitchOffLeadingOnly();
    ana->SwitchOffCheckNeutralClustersForLeading();
    
    // MC
    ana->SwitchOnPrimariesInConeSelection();
    ana->SwitchOnPrimariesPi0DecayStudy() ;
    
    if(particle.Contains("Photon"))
    {
        ana->SwitchOnDecayTaggedHistoFill() ;
        ana->SetNDecayBits(4);
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
    
    ana->SetFirstSMCoveredByTRD(-1);
    
    if(!kTM)  ana->SwitchOnTMHistoFill();
    else      ana->SwitchOffTMHistoFill();
    
    //if(!kSimulation) ana->SwitchOnFillPileUpHistograms();
    
    ana->SwitchOnRealCaloAcceptance();
    ana->SwitchOnFiducialCut();
    
    //Avoid borders of EMCal
    if(kCalorimeter=="EMCAL")
        ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.60, 86, 174) ;
    
    // Same Eta as EMCal, cut in phi if EMCAL was triggering
    if(particle=="Hadron"  || particle.Contains("CTS"))
    {
        //if(kTrig.Contains("EMC"))
        //  ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
        //else
        ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;
    }
    
    // Input / output delta AOD settings
    
    ana->SetInputAODName(Form("%s%s",particle.Data(),kName.Data()));
    if(pthMax > pth) ana->SetAODObjArrayName(Form("IC%s_%s_R%1.1f_ThMin%1.1f_ThMax%1.1f",particle.Data(),kName.Data(),cone,pth,pthMax));
    else             ana->SetAODObjArrayName(Form("IC%s_%s_R%1.1f_ThMin%1.1f"           ,particle.Data(),kName.Data(),cone,pth));
    
    
    //Do settings for main isolation cut class
    AliIsolationCut * ic =  ana->GetIsolationCut();
    ic->SetDebug(kDebug);
    ic->SetParticleTypeInCone(partInCone);
    ic->SetICMethod(thresType);
    ic->SetPtFraction(0.1);
    ic->SetPtThreshold(0.5); // default, change in next lines
    ic->SetSumPtThreshold(1.0); // default, change in next lines
    
    if(cone >0 && pth > 0)
    {
        ic->SetConeSize(cone);
        
        if(thresType == AliIsolationCut::kPtThresIC)
        {
            printf("*** Iso *** PtThresMin = %1.1f GeV/c *** PtThresMax = %1.1f GeV/c *** R = %1.1f ***\n",pth,pthMax,cone);
            ic->SetPtThreshold(pth);
            if(pthMax > pth)ic->SetPtThresholdMax(pthMax);
            else            ic->SetPtThresholdMax(100000);
        }
        
        if(thresType == AliIsolationCut::kSumPtIC)
        {
            printf("*** Iso *** SumPtMin = %1.1f GeV/c *** SumPtMax = %1.1f GeV/c *** R = %1.1f ***\n",pth,pthMax,cone);
            ic->SetSumPtThreshol(pth);
            if(pthMax > pth)ic->SetSumPtThresholdMax(pthMax);
            else            ic->SetSumPtThresholdMax(100000);
        }
    }
    else
    {
        if(kCollisions=="pp")
        {
            ic->SetPtThreshold(0.5);
            ic->SetSumPtThreshold(1.0) ;
            ic->SetConeSize(0.4);
        }
        if(kCollisions=="PbPb")
        {
            ic->SetPtThreshold(3.);
            ic->SetSumPtThreshold(3.0) ;
            ic->SetConeSize(0.3);
        }
    }
    
    
    //Do or not do isolation with previously produced AODs.
    //No effect if use of SwitchOnSeveralIsolation()
    ana->SwitchOffReIsolation();
    
    //Multiple IC
    if(multi)
    {
        ic->SetConeSize(1.);    // Take all for first iteration
        ic->SetPtThreshold(100);// Take all for first iteration
        ana->SwitchOnSeveralIsolation() ;
        ana->SetAODObjArrayName(Form("MultiIC%sTM%d",particle.Data(),kTM));
        
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
    caloPID->SetEMCALDEtaCut(0.025);
    caloPID->SetEMCALDPhiCut(0.030);
    
    //Set Histograms name tag, bins and ranges
    
    if(!multi) ana->AddToHistogramsName(Form("AnaIsol%s_TM%d_",particle.Data(),kTM));
    else       ana->AddToHistogramsName(Form("AnaMultiIsol%s_TM%d_",particle.Data(),kTM));
    
    //  if(pthMax > pth) ana->AddToHistogramsName(Form("AnaIsol%s_R%1.1f_ThMin%1.1f_ThMax%1.1f_",particle.Data(),cone,pth,pthMax));
    //  else             ana->AddToHistogramsName(Form("AnaIsol%s_R%1.1f_ThMin%1.1f_"           ,particle.Data(),cone,pth));
    
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    if(particle=="Hadron"  || particle.Contains("CTS"))
    {
        ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
        ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
    }
    
    ConfigureMC(ana);
    
    if(kPrint) ic ->Print("");
    if(kPrint) ana->Print("");
    
    return ana;
    
}

//___________________________________________________________________________________
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle,
                                                                    Bool_t bIsolated,
                                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                                    Float_t cone = 0.3,
                                                                    Float_t pth  = 0.3)
{
    
    printf("1111111111/n");
    AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
    ana->SetDebug(kDebug);
    
    ana->SetTriggerPtRange(5,100);
    ana->SetAssociatedPtRange(0.2,100);
    //ana->SetDeltaPhiCutRange( TMath::Pi()/2,3*TMath::Pi()/2 ); //[90 deg, 270 deg]
    ana->SetDeltaPhiCutRange  (TMath::DegToRad()*120.,TMath::DegToRad()*240.);
    // Underlying event
    ana->SetUeDeltaPhiCutRange(TMath::DegToRad()*60. ,TMath::DegToRad()*120.);
    ana->SwitchOnSeveralUECalculation();
    ana->SwitchOffAbsoluteLeading();  // Select trigger leading particle of all the selected tracks
    ana->SwitchOffNearSideLeading(); // Select trigger leading particle of all the particles at +-90 degrees, default
    ana->SwitchOffCheckNeutralClustersForLeading();
    
    ana->SwitchOffFillPtImbalancePerPtABinHistograms();
    ana->SwitchOffCorrelationVzBin() ;
    ana->SwitchOffFillEtaGapHistograms();
    
    ana->SwitchOffFillHighMultiplicityHistograms();
    ana->SwitchOffPi0TriggerDecayCorr();
    if(particle.Contains("Photon"))
    {
        ana->SwitchOnDecayTriggerDecayCorr();
        ana->SetNDecayBits(4);
        printf("**** SET M02 limits *** \n");
        ana->SetM02Cut(0.10,0.27);
        ana->SwitchOnInvariantMassHistograms();
        //ana->SwitchOnBackgroundBinsPtInConeHistograms();
    }
    
    ana->SetMCGenType(0,7);
    
    ana->SwitchOffLeadHadronSelection(); // Open cuts, just fill histograms
    ana->SwitchOnFillLeadHadronHistograms();
    ana->SwitchOnBackgroundBinsPtInConeHistograms();
    ana->SwitchOnBackgroundBinsTaggedDecayPtInConeHistograms();
    printf("*****SWITCH ON DECAY BIT PT IN CONE HISTO******\n");
    ana->SetLeadHadronPhiCut(TMath::DegToRad()*130, TMath::DegToRad()*230.);
    ana->SetLeadHadronPtCut(0.5, 1000);
    // if triggering on PHOS and EMCAL is on
    //if(kCalorimeter=="PHOS") ana->SwitchOnNeutralCorr();
    ana->SwitchOffNeutralCorr(); // Do only correlation with TPC
    ana->SetPi0AODBranchName("Pi0EMCAL_TrigEMC7_Cl_TM1");
    
    ana->SwitchOffHMPIDCorrelation();
    
    ana->SwitchOffFillBradHistograms();
    
    //if(!kSimulation) ana->SwitchOnFillPileUpHistograms();
    
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
    if(kMix)
    {
        ana->SwitchOnOwnMix();
        ana->SwitchOnFillNeutralInMixedEvent();
        
        if(bIsolated)
        {
            //Do settings for main isolation cut class
            AliIsolationCut * ic =  ana->GetIsolationCut();
            ic->SetDebug(kDebug);
            
            if(cone >0 && pth > 0)
            {
                printf("*** Correl *** PtThres = %1.1f GeV/c *** R = %1.1f ***\n",pth,cone);
                ic->SetPtThreshold(pth);
                ic->SetConeSize(cone);
            }
            else
            {
                if(kCollisions=="pp")
                {
                    ic->SetPtThreshold(0.5);
                    ic->SetConeSize(0.4);
                }
                if(kCollisions=="PbPb")
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
    
    if(kCollisions=="pp")
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
        if(kName.Contains("60_90"))
        {
            printf("*** Set mixing for peripheral\n");
            ana->SetNMaxEvMix(50);
            ana->SetNCentrBin(2);
        }
    }
    
    ana->SwitchOnFiducialCut();
    
    //Avoid borders of EMCal, same as for isolation
    if(kCalorimeter=="EMCAL")
        ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;
    
    // Same Eta as EMCal, cut in phi if EMCAL was triggering
    if(particle=="Hadron" || particle.Contains("CTS"))
    {
        //if(kTrig.Contains("EMC"))
        //  ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
        //else
        ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;
    }
    
    // Input / output delta AOD settings
    
    ana->SetInputAODName(Form("%s%s",particle.Data(),kName.Data()));
    ana->SetAODObjArrayName(Form("%sHadronCorrIso%d_%s",particle.Data(),bIsolated,kName.Data()));
    
    
    //Set Histograms name tag, bins and ranges
    
    ana->AddToHistogramsName(Form("Ana%sHadronCorr_Iso%d_TM%d_",particle.Data(),bIsolated,kTM));
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    if(particle=="Hadron"  || particle.Contains("CTS"))
    {
        ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
        ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
    }
    
    ConfigureMC(ana);
    
    if(kPrint) ana->Print("");
    printf("aaaaaaaaaaa\n");
    
    return ana;
    
}

//________________________________________________________________________________
AliAnaElectron* ConfigureElectronAnalysis()
{
    
    AliAnaElectron *ana = new AliAnaElectron();
    ana->SetDebug(kDebug); //10 for lots of messages
    
    ana->FillAODWithElectrons() ;
    //ana->FillAODWithHadrons();
    //ana->FillAODWithAny();
    
    if(kCalorimeter == "PHOS")
    {
        ana->SetNCellCut(2);// At least 2 cells
        ana->SetMinPt(0.3);
        ana->SetMinDistanceToBadChannel(2, 4, 5);
    }
    else
    {//EMCAL
        ana->SetNCellCut(1);// At least 2 cells
        ana->SetMinPt(0.5); // no effect minium EMCAL cut.
        ana->SetMaxPt(100);
        //ana->SetTimeCut(400,900);// Time window of [400-900] ns
        ana->SetMinDistanceToBadChannel(2, 4, 6);
    }
    
    //Electron selection cuts with tracks
    ana->SetEOverP(0.85, 1.2);
    
    // TO DO, find a more suitable way to set this
    if     (kRunNumber < 146861) ana->SetdEdxCut(72, 90);
    else if(kRunNumber < 154000) ana->SetdEdxCut(54, 70);
    else                         ana->SetdEdxCut(74, 90);
    
    if(kSimulation)  ana->SetdEdxCut(80, 100);
    
    ana->SetCalorimeter(kCalorimeter);
    
    ana->SwitchOnCaloPID();
    
    AliCaloPID* caloPID = ana->GetCaloPID();
    
    caloPID->SetEMCALLambda0CutMax(0.27);
    caloPID->SetEMCALLambda0CutMin(0.10);
    
    ana->SwitchOffFillShowerShapeHistograms();
    ana->SwitchOffFillWeightHistograms()  ;
    ana->SwitchOffFiducialCut();
    
    
    if(!kData.Contains("delta"))
    {
        ana->SetOutputAODName(Form("Electron%s",kName.Data()));
        ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    }
    else ana->SetInputAODName(Form("Electron%s",kName.Data()));
    
    //Set Histograms name tag, bins and ranges
    
    ana->AddToHistogramsName(Form("AnaElectron_TM%d_",kTM));
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    ConfigureMC(ana);
    
    if(kPrint) ana->Print("");
    
    return ana ;
    
}


//_______________________________________________
AliAnaChargedParticles* ConfigureChargedAnalysis()
{
    
    AliAnaChargedParticles *ana = new AliAnaChargedParticles();
    ana->SetDebug(kDebug); //10 for lots of messages
    
    // selection cuts
    
    ana->SetMinPt(0.2);
    ana->SwitchOnFiducialCut();
    ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ; //more restrictive cut in reader and after in isolation
    
    ana->SwitchOffFillVertexBC0Histograms();
    //if(!kSimulation)
    ana->SwitchOffFillPileUpHistograms();
    ana->SwitchOffFillTrackBCHistograms();
    
    // Input / output delta AOD settings
    
    if(!kData.Contains("delta"))
    {
        ana->SetOutputAODName(Form("Hadron%s",kName.Data()));
        //ana->SetOutputAODClassName("AliAODCaloTrackParticleCorrelation");
        //ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
        ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
    }
    else
        ana->SetInputAODName(Form("Hadron%s",kName.Data()));
    printf("Set Hadron%s\n",kName.Data());
    //Set Histograms name tag, bins and ranges
    
    ana->AddToHistogramsName("AnaHadrons_");
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
    
    ConfigureMC(ana);
    
    if(kPrint) ana->Print("");
    
    return ana;
    
}

//________________________________________
AliAnaCalorimeterQA* ConfigureQAAnalysis()
{
    
    AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
    ana->SetDebug(kDebug); //10 for lots of messages
    ana->SetCalorimeter(kCalorimeter);
    
    ana->SetTimeCut(-1e10,1e10); // Open time cut
    
    // Study inter detector correlation (PHOS, EMCAL, Tracks, V0)
    if(kCalorimeter=="PHOS"  && kTrig=="PHOS")
        ana->SwitchOnCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
    if(kCalorimeter=="EMCAL" && kClusterArray=="")
        ana->SwitchOnCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
    else
        ana->SwitchOffCorrelation();
    
    // Study exotic clusters PHOS and EMCAL
    if(kClusterArray=="") ana->SwitchOnStudyBadClusters() ;
    else                  ana->SwitchOffStudyBadClusters() ;
    
    ana->SwitchOnRealCaloAcceptance();
    
    ana->SwitchOffFiducialCut();
    ana->SwitchOffFillAllTH3Histogram();
    ana->SwitchOffFillAllPositionHistogram();
    ana->SwitchOffFillAllPositionHistogram2();
    if(!kExotic)ana->SwitchOnStudyBadClusters();
    else        ana->SwitchOffStudyBadClusters();
    ana->SwitchOffStudyClustersAsymmetry();
    ana->SwitchOffStudyWeight();
    ana->SwitchOnFillAllTrackMatchingHistogram();
    ana->SwitchOnFillAllCellTimeHisto() ;
    
    ana->AddToHistogramsName("QA_"); //Begining of histograms name
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    ConfigureMC(ana);
    
    if(kPrint) ana->Print("");
    
    return ana;
    
}

//________________________________________________________________
AliAnaGeneratorKine* ConfigureGenKineAnalysis()
{
    // Analysis for parton, jets correlation with photon and pi0
    
    AliAnaGeneratorKine *ana = new AliAnaGeneratorKine();
    ana->SetDebug(kDebug); //10 for lots of messages
    
    // Trigger detector, acceptance and pT cut
    ana->SetTriggerDetector("EMCAL");
    ana->SetMinPt(10); // Trigger photon, pi0 minimum pT
    ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.4, 90, 170);
    
    // Particles associated to trigger or isolation cone acceptance and pT cut
    ana->SetCalorimeter("EMCAL");
    ana->SetMinChargedPt(0.2);
    ana->SetMinNeutralPt(0.3);
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.65, 81, 179);
    ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360);
    
    // Isolation paramters
    AliIsolationCut * ic =  ana->GetIsolationCut();
    ic->SetDebug(kDebug);
    ic->SetPtThreshold(0.5);
    ic->SetConeSize(0.4);
    ic->SetSumPtThreshold(1.0) ;
    ic->SetICMethod(AliIsolationCut::kSumPtIC); // kSumPtIC
    
    ana->AddToHistogramsName("AnaGenKine_");
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    if(kPrint) ana->Print("");
    
    return ana;
    
}

//________________________________________________________________
AliAnaEMCALTriggerClusters* ConfigureEMCALTriggerClusterAnalysis()
{
    // For filling all histograms meaninfully, in the reader, time cut must be off
    // and bad triggered events not rejected, and of course analyze triggered events.
    
    AliAnaEMCALTriggerClusters *ana = new AliAnaEMCALTriggerClusters();
    ana->SetDebug(kDebug); //10 for lots of messages
    
    // cluster selection cuts
    
    ana->SwitchOffFiducialCut();
    ana->SetNCellCut(1);// At least 2 cells
    ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
    ana->SetMaxEnergy(1000);
    ana->SetM02(1, 2) ;
    ana->SwitchOnTrackMatchRejection() ;
    
    ana->AddToHistogramsName("EMCTriggerClusters_");
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    if(kPrint) ana->Print("");
    
    return ana;
    
}

//___________________________________________________
AliAnaClusterPileUp* ConfigureClusterPileUpAnalysis()
{
    // For filling all histograms meaninfully, in the reader, time cut must be off
    // and bad triggered events in different BC not rejected
    
    AliAnaClusterPileUp *ana = new AliAnaClusterPileUp();
    ana->SetDebug(kDebug); //10 for lots of messages
    
    // cluster selection cuts
    
    ana->SwitchOffFiducialCut();
    ana->SetNCellCut(1);// At least 2 cells
    ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
    ana->SetMaxEnergy(1000);
    
    ana->AddToHistogramsName("ClusterPileUp_");
    SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
    if(kPrint) ana->Print("");
    
    return ana;
    
}

//________________________________________________________
void ConfigureMC(AliAnaCaloTrackCorrBaseClass* ana)
{
    if(kSimulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
    else            ana->SwitchOffDataMC() ;
    
    //Set here generator name, default pythia
    //ana->GetMCAnalysisUtils()->SetMCGenerator("");
}

//________________________________________________________
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges)
{
    // Set common bins for all analysis and MC histograms filling
    
    histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
    
    if(kCalorimeter=="EMCAL")
    {
        if(kYears==2010)
        {
            histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
            histoRanges->SetHistoXRangeAndNBins(-230,90,120); // QA
            histoRanges->SetHistoYRangeAndNBins(370,450,40);  // QA
        }
        else if(kYears==2011)
        {
            histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
            histoRanges->SetHistoXRangeAndNBins(-600,90,200); // QA
            histoRanges->SetHistoYRangeAndNBins(100,450,100); // QA
        }
        else
        {
            histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 190*TMath::DegToRad(), 122) ;
            histoRanges->SetHistoXRangeAndNBins(-100,90,200); // QA
            histoRanges->SetHistoYRangeAndNBins(50,450,100);  // QA
        }
        
        histoRanges->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
    }
    else
    {
        histoRanges->SetHistoPhiRangeAndNBins(260*TMath::DegToRad(), 320*TMath::DegToRad(), 60) ;
        histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
    }
    
    histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 4.9, 500);
    
    // Invariant mass histoRangeslysis
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
    histoRanges->SetHistoZRangeAndNBins(-400,400,200);
    histoRanges->SetHistoRRangeAndNBins(400,450,25);
    histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
    histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
    
    // QA, correlation
    if(kCollisions=="PbPb")
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
    
}

//_____________________________
UInt_t SetTriggerMaskFromName()
{
    if(kTrig=="EMC7")
    {
        printf("CaloTrackCorr trigger EMC7\n");
        return AliVEvent::kEMC7;
    }
    else if (kTrig=="INT7")
    {
        printf("CaloTrackCorr trigger INT7\n");
        return AliVEvent::kINT7;
    }
    else if(kTrig=="EMC1")
    {
        printf("CaloTrackCorr trigger EMC1\n");
        return AliVEvent::kEMC1;
    }
    else if(kTrig=="MB")
    {
        printf("CaloTrackCorr trigger MB\n");
        return AliVEvent::kMB;
    }  
    else if(kTrig=="PHOS")
    {
        printf("CaloTrackCorr trigger PHOS\n");
        return AliVEvent::kPHI7;
    }  
    else if(kTrig=="PHOSPb")
    {
        printf("CaloTrackCorr trigger PHOSPb\n");
        return AliVEvent::kPHOSPb;
    }
    else if(kTrig=="AnyINT")
    {
        printf("CaloTrackCorr trigger AnyINT\n");
        return AliVEvent::kAnyINT;
    }  
    else if(kTrig=="INT")
    {
        printf("CaloTrackCorr trigger AnyINT\n");
        return AliVEvent::kAny;
    }
    else if(kTrig=="EMCEGA")
    {
        printf("CaloTrackCorr trigger EMC Gamma\n");
        return AliVEvent::kEMCEGA;
    } 
    else if(kTrig=="EMCEJE")
    {
        printf("CaloTrackCorr trigger EMC Jet\n");
        return AliVEvent::kEMCEJE;
    }
    else if(kTrig=="Central")
    {
        printf("CaloTrackCorr trigger Central\n");
        return AliVEvent::kCentral;
    }
    else if(kTrig=="CentralEGA")
    {
        printf("CaloTrackCorr trigger Central+EMCEGA\n");
        return (AliVEvent::kCentral | AliVEvent::kEMCEGA);
    }
    else if(kTrig=="SemiCentral")
    {
        printf("CaloTrackCorr trigger SemiCentral\n");
        return AliVEvent::kSemiCentral;
    }
    else if(kTrig=="SemiOrCentral")
    {
        printf("CaloTrackCorr trigger SemiCentral Or Central\n");
        return (AliVEvent::kSemiCentral | AliVEvent::kCentral);
    }  
}





