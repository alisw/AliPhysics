/// \file AddTaskClusterShape.C
/// \ingroup CaloTrackCorrMacros
/// \brief Configuration of EMCal cluster shape analysis
///
/// Configuration macro for analysis of EMCal cluster shape analysis. Many histograms correlating different
/// shower shape parameters and correlations are filled.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

// Set includes for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

// ROOT
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>

// ALIROOT/ALIPHYSICS
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVTrack.h"
#include "AliLog.h"

// CaloTrackCorrr
#include "AliAnalysisTaskCaloTrackCorrelation.h"
#include "AliCaloTrackESDReader.h"
#include "AliCaloTrackAODReader.h"
#include "AliCalorimeterUtils.h"
#include "AliAnaClusterShapeCorrelStudies.h"
#include "AliHistogramRanges.h"
#include "AliAnaCalorimeterQA.h"
#include "AliAnaCaloTrackCorrMaker.h"

// Macros
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#include "ConfigureEMCALRecoUtils.C"
#include "PWGGA/CaloTrackCorrelations/macros/ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C"
#include "PWGGA/CaloTrackCorrelations/macros/CheckActiveEMCalTriggerPerPeriod.C"
#include "PWGGA/CaloTrackCorrelations/macros/GetAlienGlobalProductionVariables.C"

#endif // CINT

 // Declare methods for compilation

AliCaloTrackReader  * ConfigureReader        
(TString col,           Bool_t simulation,
 TString clustersArray, Bool_t tender,
 TString calorimeter,   Bool_t nonLinOn,
 TString trigger,       Bool_t rejectEMCTrig,
 Int_t   minCen,        Int_t  maxCen,
 Bool_t  printSettings, Int_t  debug       );

AliCalorimeterUtils * ConfigureCaloUtils     
(TString col,           Bool_t simulation,
 Bool_t tender,
 Bool_t  nonLinOn,      Int_t  year,
 Bool_t  printSettings, Int_t  debug );

AliAnaClusterShapeCorrelStudies* ConfigureClusterShape
(Bool_t tmPtDep, TString col , Bool_t  simulation, 
 TString calorimeter, Int_t year, Bool_t  printSettings, Int_t   debug);

AliAnaCalorimeterQA * ConfigureQAAnalysis    
(TString col,           Bool_t  simulation,
 TString calorimeter,   Int_t   year,
 Bool_t  printSettings, Int_t   debug );

void SetAnalysisCommonParameters             
(AliAnaCaloTrackCorrBaseClass* ana,
 TString calorimeter , Int_t  year,
 TString col         , Bool_t simulation,
 Bool_t printSettings, Int_t  debug);


/// Global name to be composed of the settings, used to set the AOD branch name
TString kAnaClusterShape = "";

///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle.
/// \param simulation : A bool identifying the data as simulation.
/// \param collision: A string with the colliding system. If empty, alien environment used.
/// \param period : A string with the data period: LHC11h, LHC15n ... from it we extract the year.
/// \param rejectEMCTrig : An int to reject EMCal triggered events with bad trigger: 0 no rejection, 1 old runs L1 bit, 2 newer runs L1 bit
/// \param clustersArray : A string with the array of clusters not being the default (default is empty string)
/// \param tender : A bool indicating if the tender was running before this analysis
/// \param nonLinOn : A bool to set the use of the non linearity correction
/// \param tmDep: Track matching option
/// \param qaAn : A bool to switch the calorimeter QA analysis
/// \param outputfile : A string to change the name of the histograms output file, default is AnalysisResults.root
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param trigSuffix :  A string with the trigger class, abbreviated, defined in ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C
///
AliAnalysisTaskCaloTrackCorrelation * AddTaskClusterShape
(
 TString  calorimeter   = "EMCAL",
 Bool_t   simulation    = kFALSE,
 TString  collision     = "pp",
 TString  period        = "",
 Int_t    rejectEMCTrig = 0,
 TString  clustersArray = "",
 Bool_t   tender        = kFALSE,
 Bool_t   nonLinOn      = kFALSE,
 Bool_t   tmDep         = kTRUE,
 Bool_t   qaAn          = kFALSE,
 TString  outputfile    = "",
 Bool_t   printSettings = kFALSE,
 Int_t    debug         = 0,  // Debug level
 const char *trigSuffix = "EMC7"
)
{
  // Check the global variables, and reset the provided ones if empty.
  //
  TString trigger  = trigSuffix;
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/GetAlienGlobalProductionVariables.C");
  
  Int_t   year        = 2017;
  Bool_t  printGlobal = kTRUE;
 
  GetAlienGlobalProductionVariables(simulation,collision,period,year,printGlobal);

  // Get the pointer to the existing analysis manager via the static access method.
  //  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskClusterShape", "No analysis manager to connect to.");
    return NULL;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskClusterShape", "This task requires an input event handler");
    return NULL;
  }
  
  //
  // Create task
  //
  
  // Name for containers
  kAnaClusterShape = Form("ShapeAna%s_%s",trigger.Data(),clustersArray.Data());
  
  Int_t  minCen = -1, maxCen = -1;
  
  if(collision=="PbPb" && maxCen>=0) kAnaClusterShape+=Form("Cen%d_%d",minCen,maxCen);
  
  // Create task
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("%s",kAnaClusterShape.Data()));
  
  task->SetDebugLevel(debug);
  
  //  task->SetFirstEvent(1800);
  //  task->SetLastEvent (2000);  
  
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  
  task->SetAnalysisMaker(maker);

  mgr->AddTask(task);
  
  //
  // Create containers  
  
  if(outputfile.Length()==0) outputfile = AliAnalysisManager::GetCommonFileName();

  TString containerName = Form("Shape_%s",calorimeter.Data());
  if(clustersArray.Length()>0) containerName = Form("%s_%s",containerName.Data(),clustersArray.Data());
  
  TString subcontainerName = Form("%s",trigger.Data());
  if(clustersArray.Length()>0) subcontainerName = Form("%s_%s",subcontainerName.Data(),clustersArray.Data());
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(subcontainerName, TList::Class(),
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s:%s",outputfile.Data(),containerName.Data()));
  
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",subcontainerName.Data()), TList::Class(),
                                                             AliAnalysisManager::kParamContainer, 
                                                             Form("%s_Parameters.root",containerName.Data()));

  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  //==============================================================================
  
  // Do not configure the wagon for certain analysis combinations
  // But create the task so that the sub-wagon train can run
  //
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/CheckActiveEMCalTriggerPerPeriod.C");
  Bool_t doAnalysis = CheckActiveEMCalTriggerPerPeriod(simulation,trigger,period,year);
  if(!doAnalysis) 
  {
    maker->SwitchOffProcessEvent();
    return task;
  }

  // #### Start analysis configuration ####
  //  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Make sure the B field is enabled for track selection, some cuts need it
  //
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);

  // Print settings to check all is as expected
  //
  printf("AddTaskClusterShape() - Task NAME: %s; Passed settings\n",kAnaClusterShape.Data());

  printf("\t calorimeter <%s>, simulation <%d>, year <%d>, collision <%s>, trigger <%s>, reject EMC <%d>,\n"
         "\t clustersArray <%s>, tender <%d>, non linearity <%d>, TM dep <%d>, QA on <%d>,\n"
         "\t outputfile <%s>, printSettings <%d>, debug <%d>\n",
         calorimeter.Data(),simulation,year,collision.Data(),trigger.Data(), rejectEMCTrig, 
         clustersArray.Data(),tender, nonLinOn, tmDep, qaAn, outputfile.Data(),printSettings,debug);

  
  // General frame setting and configuration
  maker->SetReader   ( ConfigureReader   (collision,simulation,clustersArray,tender,calorimeter,nonLinOn,trigger,rejectEMCTrig,minCen,maxCen,printSettings,debug) );
  maker->SetCaloUtils( ConfigureCaloUtils(collision,simulation,tender,nonLinOn,year,printSettings,debug) );
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important
  
  // SHAPE
  maker->AddAnalysis(ConfigureClusterShape(tmDep, collision,simulation,calorimeter,year,printSettings,debug) , n++);
  
 // Detector QA
 if (qaAn) maker->AddAnalysis(ConfigureQAAnalysis(collision,simulation,calorimeter,year,printSettings,debug) , n++);

  
  maker->SetAnaDebug(debug)  ;
  
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  
  //if( simulation || !trigger.Contains("EMC") ) 
  maker->SwitchOnDataControlHistograms();
  
  if(simulation)
  {
    // Calculate the cross section weights, apply them to all histograms 
    // and fill xsec and trial histo. Sumw2 must be activated.
    //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionCalculation(); 
    //maker->SwitchOnSumw2Histograms();
        
    // Just fill cross section and trials histograms.
    maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionHistoFill(); 
    
    // For recent productions where the cross sections and trials are not stored in separate file
    TString prodType = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTYPE");
    printf("AddTaskClusterShape() - MC production name: %s\n",prodType.Data());
    if ( prodType.Contains("LHC16c") ) // add here any other affected periods, for the moment jet-jet 8 TeV
    {   
      printf("\t use the cross section from EventHeader per Event\n");
      maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionFromEventHeader() ;
    }
    
    // Add control histogram with pT hard to control aplication of weights 
    maker->SwitchOnPtHardHistogram();
  }
  
  if(printSettings) maker->Print("");
  
  printf("AddTaskClusterShape() - << End Configuration of %d analysis for calorimeter %s >>\n",n, calorimeter.Data());
  
  
  //
  // Select events trigger depending on trigger
  //
  if(!simulation)
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C");
    TString caloTriggerString = "";
    UInt_t mask = ConfigureAndGetEventTriggerMaskAndCaloTriggerString(trigger, year, caloTriggerString);
    
    maker->GetReader()->SetFiredTriggerClassName(caloTriggerString);
    
    task ->SelectCollisionCandidates( mask );
  }
  
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
                                     Bool_t printSettings,  Int_t   debug        )
{
  // Get the data type ESD or AOD
  AliAnalysisManager * mgr = AliAnalysisManager::GetAnalysisManager();
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();
  
  AliCaloTrackReader * reader = 0;
  if     (inputDataType == "AOD") reader = new AliCaloTrackAODReader();
  else if(inputDataType == "ESD") reader = new AliCaloTrackESDReader();
  else printf("AddTaskClusterShape::ConfigureReader() - Data not known InputData=%s\n",inputDataType.Data());
  
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
    reader->SetPtHardAndJetPtFactor(2);
    
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
  
  reader->SetEMCALBadChannelMinDist(5);
  reader->SetEMCALNCellsCut(2);
  
  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();
  
  // Shower shape smearing
  // Set it in the train configuration page not here for the moment
  //  if(simulation)
  //  {
  //    reader->SwitchOffShowerShapeSmearing(); // Active only on MC, off by default
  //    reader->SetShowerShapeSmearWidth(0.005);
  //  }
  reader->SwitchOffShowerShapeSmearing();
  
  reader->SwitchOffFiducialCut();
  
  //
  // No track in this analysis
  reader->SwitchOffCTS();
  
  //
  // Calorimeter
  //
  if(clustersArray == "" && !tender)
  {
    printf("AddTaskClusterShape::ConfigureReader() - **************** Standard EMCAL clusters branch analysis **************** \n");
    reader->SwitchOnClusterRecalculation();
    // Check in ConfigureCaloUtils that the recalibration and bad map are ON
  }
  else
  {
    printf("AddTaskClusterShape::ConfigureReader() - **************** Input for analysis is Clusterizer %s **************** \n", clustersArray.Data());
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
    printf("AddTaskClusterShape::ConfigureReader() === Remove bad triggers === \n");
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
      printf("AddTaskClusterShape::ConfigureReader() - Trigger cluster calibration OFF\n");
      reader->SwitchOffTriggerClusterTimeRecal() ;
    }
    
  }
  
  //reader->RejectFastClusterEvents() ;
  
  // For mixing with AliAnaParticleHadronCorrelation switch it off
  reader->SwitchOnEventTriggerAtSE(); // on is default case
  
  reader->SetZvertexCut(10.);               // Open cut
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
  
  reader->SwitchOffV0ANDSelection() ;       // and besides v0 AND
  reader->SwitchOnPileUpEventRejection();   // remove pileup by default off, apply it only for MB not for trigger
  
  if(col=="PbPb")
  {
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
    reader->SetCentralityBin(minCen,maxCen); // Accept all events, if not select range
    // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
    reader->SetEventPlaneMethod("V0");
  }
  
  if(printSettings) reader->Print("");
  
  return reader;
}

///
/// Configure the class handling the calorimeter clusters specific methods
///
AliCalorimeterUtils* ConfigureCaloUtils(TString col,           Bool_t simulation,
                                        Bool_t tender,
                                        Bool_t  nonLinOn,      Int_t year,
                                        Bool_t  printSettings, Int_t   debug)
{
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  
  cu->SetDebug(debug);
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away
  cu->SetNumberOfCellsFromEMCALBorder(3);
  cu->SetNumberOfCellsFromPHOSBorder (2);
  
  cu->SetNumberOfSuperModulesUsed(10);
  
  if     (year == 2010) cu->SetNumberOfSuperModulesUsed(4);
  else if(year <= 2013) cu->SetNumberOfSuperModulesUsed(10);
  else if(year >  2013) cu->SetNumberOfSuperModulesUsed(20);
  else                  cu->SetNumberOfSuperModulesUsed(10);
  
  printf("AddTaskClusterShape::ConfigureCaloUtils() xxx Number of SM set to <%d> xxx\n",
         cu->GetNumberOfSuperModulesUsed());
  
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
  
AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
//
//  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
//  ConfigureEMCALRecoUtils(recou,
//                          simulation,
//                          kTRUE,      // exotic
//                          nonLinOn,   // Non linearity
//                          calibEner,  // E calib
//                          kTRUE,      // bad map
//                          calibTime); // time calib
  
  cu->ConfigureEMCALRecoUtils(simulation,
                              kTRUE,      // exotic
                              nonLinOn,   // Non linearity
                              calibEner,  // E calib
                              kTRUE,      // bad map
                              calibTime); // time calib
  
  //if( calibTime ) recou->SetExoticCellDiffTimeCut(1e6);
  
  if( nonLinOn )  cu->SwitchOnCorrectClusterLinearity();
  
  printf("AddTaskClusterShape::ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
  printf("AddTaskClusterShape::ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
  
  // PHOS
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
  
  if(printSettings) cu->Print("");
  
  return cu;
}

///
/// Configure the task doing cluster shape studies
///
AliAnaClusterShapeCorrelStudies* ConfigureClusterShape
(Bool_t tmPtDep, TString col , Bool_t  simulation, 
 TString calorimeter, Int_t year, Bool_t  printSettings, Int_t   debug)
{
  AliAnaClusterShapeCorrelStudies *ana = new AliAnaClusterShapeCorrelStudies();
  
  ana->SetM02Min(-1);
  ana->SetNCellsPerClusterMin(-1);
  
  ana->SetCalorimeter(calorimeter);
  
  if(simulation) ana->SetConstantTimeShift(615);
  
  ana->SwitchOffFiducialCut();
    
  ana->SwitchOnStudyClusterShape();
  
  ana->SwitchOnStudyEMCalModuleCells();
  
  ana->SwitchOffStudyClusterShapeParam();
  
  ana->SwitchOffStudyMatchedPID() ;
  
  ana->SwitchOffStudyWeight();
  
  ana->SetNCellBinLimits(-1); // no analysis on predefined bins in nCell
  
  ana->SwitchOffStudyTCardCorrelation() ;
  ana->SwitchOffStudyExotic();
  ana->SwitchOffStudyInvariantMass();
  ana->SwitchOffStudyColRowFromCellMax() ;
  ana->SwitchOffStudyCellTime() ;
  
  // PID cuts (Track-matching)
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
    
  if(tmPtDep)
  {
    // track pT dependent cut
    caloPID->SwitchOnEMCTrackPtDepResMatching();
    
    // Begining of histograms name
    ana->AddToHistogramsName("Shape_TMDep_");
  }
  else
  {
    // Fix
    caloPID->SwitchOffEMCTrackPtDepResMatching();
    caloPID->SetEMCALDEtaCut(0.025);
    caloPID->SetEMCALDPhiCut(0.030);
    
    // Begining of histograms name
    ana->AddToHistogramsName("Shape_TMFix_"); 
  }
    
  SetAnalysisCommonParameters(ana,calorimeter,year,col,simulation,printSettings,debug); // see method below
  
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
  
  ana->SetEMCALM02Min(-1);
  ana->SetEMCALNCellsPerClusterMin(-1);
  ana->SetEMCALCellAmpMin(0.1);
  
  ana->SetCalorimeter(calorimeter);
  
  ana->SetTimeCut(-1e10,1e10); // Open time cut
  if(simulation) ana->SetConstantTimeShift(615);
  
  ana->SwitchOffCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
  
  ana->SwitchOnRealCaloAcceptance();
  
  ana->SwitchOffFiducialCut();
  ana->SwitchOffFillAllTH3Histogram();
  ana->SwitchOffFillAllPositionHistogram();
  ana->SwitchOffFillAllPositionHistogram2();
  ana->SwitchOffStudyBadClusters() ;
  //  ana->SwitchOffStudyClustersAsymmetry();
  //  ana->SwitchOffStudyWeight();
  //  ana->SwitchOnStudyTCardCorrelation() ;
  //  ana->SwitchOffStudyM02Dependence() ;
  ana->SwitchOnFillAllTrackMatchingHistogram();
  //  ana->SwitchOnStudyExotic();
  
  ana->SwitchOffFillAllCellTimeHisto() ;
  
  ana->SwitchOnFillAllCellHistogram();
    
  ana->AddToHistogramsName("QA_"); // Begining of histograms name
  
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
  
  histoRanges->SetHistoPtRangeAndNBins(0, 75, 150) ; // Energy and pt histograms
  
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
  
  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 2.9, 300);
  
  // Invariant mass histo
  histoRanges->SetHistoMassRangeAndNBins(0., 1., 200) ;
  histoRanges->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  
  // check if time calibration is on
  //histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  histoRanges->SetHistoTimeRangeAndNBins(-400.,400,400);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.06,0.06,120);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.06,0.06,120);
  histoRanges->SetHistodRRangeAndNBins(0.,0.06,60);//QA
  
  // QA, electron, charged
  histoRanges->SetHistoEOverPRangeAndNBins(0,1.5,150);
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
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,150,150);
  }
  
  // xE, zT
  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,200);
  histoRanges->SetHistoHBPRangeAndNBins  (0.,10.,200);
  
  // Isolation
  histoRanges->SetHistoPtInConeRangeAndNBins(0, 50 , 100);
  histoRanges->SetHistoPtSumRangeAndNBins   (0, 100, 100);
  
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

  //
  // Debug
  //
  if(printSettings) ana->Print("");
  
  ana->SetDebug(debug); // 10 for lots of messages
}

