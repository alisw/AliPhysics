/// \file AddTaskMultipleTrackCutIsoConeAnalysis.C
/// \ingroup CaloTrackCorrMacros
/// \brief Configuration of (isolated) gamma/pi0-hadron analysis with multiple cuts
///
/// Configuration macro for analysis of gamma-hadron and pi0-hadron correlation analysis
/// both pi0 and gamma isolated (or not) with multiple cuts (Track matching, isolation cone etc).
/// It calls 
///   * the main task intializationg the base of the analysis AddTaskCaloTrackBase.C 
///   * a configuration file for the analysis, as many times as cuts needed to vary
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

// Set includes for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

// Root
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>

// Analysis
#include "AliIsolationCut.h"
#include "AliAnalysisTaskCaloTrackCorrelation.h"

// Macros
#include "AddTaskCaloTrackCorrBase.C"
#include "ConfigureCaloTrackCorrAnalysis.C"
#include "GetAlienGlobalProductionVariables.C"

#endif

///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle: EMCAL, DCAL, PHOS
/// \param simulation : A bool identifying the data as simulation
/// \param year: The year the data was taken, used to configure some histograms
/// \param col: A string with the colliding system
/// \param period: A string with data period 
/// \param rejectEMCTrig : An int to reject EMCal triggered events with bad trigger: 0 no rejection, 1 old runs L1 bit, 2 newer runs L1 bit
/// \param clustersArray : A string with the array of clusters not being the default (default is empty string)
/// \param gloCutsString : A string with list of global cuts/parameters (activate pile-up ...)
/// \param nonLinOn : A bool to set the use of the non linearity correction
/// \param analysisString : String that contains what analysis to activate, options: Photon, DecayPi0, MergedPi0, Charged, QA, Isolation, Correlation, Generator
/// \param shshMax : A float setting the maximum value of the shower shape of the clusters for the correlation analysis
/// \param isoCone : A float setting the isolation cone size higher limit
/// \param isoPtTh : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param isoMethod : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
/// \param isoContent : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
/// \param leading : An int setting the type of leading particle selection: 0, select all;l 1: absolute  leading of charged; 2: absolute  leading of charged and neutral; 3: near side leading absolute of charged; 4: near side leading absolute of charged and neutral
/// \param minCen : An int to select the minimum centrality, -1 means no selection
/// \param maxCen : An int to select the maximum centrality, -1 means no selection
/// \param mixOn : A bool to switch the correlation mixing analysis
/// \param calibrate : Use own calibration tools, do not rely on EMCal correction framewor or clusterizer
/// \param outputfile : A string to change the name of the histograms output file, default is AnalysisResults.root
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
/// \param trigSuffix :  A string with the trigger class, abbreviated, defined in ConfigureAndGetEventTriggerMaskAndCaloTriggerString.C
///
AliAnalysisTaskCaloTrackCorrelation * AddTaskMultipleTrackCutIsoConeAnalysis
(
 TString  calorimeter   = "EMCAL", // "DCAL", "PHOS"
 Bool_t   simulation    = kFALSE,
 Int_t    year          = -1,
 TString  col           = "",
 TString  period        = "",
 Int_t    rejectEMCTrig = 0,
 TString  clustersArray = "",
 TString  gloCutsString = "",
 Bool_t   calibrate     = kFALSE,
 Bool_t   nonLinOn      = kFALSE,
 TString  analysisString= "Photon_MergedPi0_DecayPi0_Isolation_Correlation_QA_Charged",
 Float_t  shshMax       = 0.27,
 Float_t  isoCone       = 0.4,
 Float_t  isoPtTh       = 2,
 Int_t    isoMethod     = AliIsolationCut::kSumPtIC,
 Int_t    isoContent    = AliIsolationCut::kNeutralAndCharged,
 Int_t    leading       = 0,
 Int_t    minCen        = -1,
 Int_t    maxCen        = -1,
 Bool_t   mixOn         = kTRUE,
 TString  outputfile    = "",
 Bool_t   printSettings = kFALSE,
 Int_t    debug         = 0,  
 const char *trigSuffix = "EMC7"
)
{
  printf("AddTaskMultipleTrackCutIsoConeAnalysis::Start configuration\n");
  
#if defined(__CINT__)
  
  printf("AddTaskMultipleTrackCutIsoConeAnalysis::Load macros\n");
  // Load macros
  //
  // Load macros
  //
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/AddTaskCaloTrackCorrBase.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/ConfigureCaloTrackCorrAnalysis.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/GetAlienGlobalProductionVariables.C");

#endif

  // First check the ALIEN environment settings
  //
  GetAlienGlobalProductionVariables(simulation,col,period,year,kTRUE);
  
  // Init base task
  //
  AliAnalysisTaskCaloTrackCorrelation * task = AddTaskCaloTrackCorrBase
  (calorimeter, simulation, year, col, period, rejectEMCTrig, clustersArray, gloCutsString,
   calibrate, nonLinOn, minCen, maxCen, mixOn, outputfile, printSettings, debug, trigSuffix);
  
  if ( !task ) return NULL;
  
  // No need to continue configuration if event is not processed
  if ( !task->GetAnalysisMaker()->IsEventProcessed() ) return task ;
  
  TList * anaList = task->GetAnalysisMaker()->GetListOfAnalysisContainers();
  printf("TList name: %s\n",anaList->GetName());
  
  // Configure the different analysis
  //
  Float_t rMin[] = {-1,0.05,0.1};
  
  // 3 exclusion of isolation cone options, no cut and 2 cuts
  for(Int_t irmin = 0; irmin < 3; irmin++)
  {
    // Add this string to deacticate the correlation without isolation
    // and other histograms in the Photon and Pi0 selection task
    if(irmin == 1 ) analysisString+="_MultiIso";
    
    // Test 3 track matching options (no track mathing)
    for(Int_t itm = 0; itm < 3; itm++)
    {
      
      TString histoString = Form("TM%d",itm);
      if ( irmin > 0 ) histoString+=Form("_Rmin%1.2f",rMin[irmin]);
      
      
      ConfigureCaloTrackCorrAnalysis
      ( anaList, calorimeter, simulation, year, col, analysisString, histoString, 
       shshMax, isoCone, rMin[irmin], isoPtTh, isoMethod, isoContent,
       leading, itm, mixOn, printSettings, debug);
    }
  }
  
  // Execute some control task only
  ConfigureCaloTrackCorrAnalysis
  ( anaList, calorimeter, simulation, year, col, "QA_Charged", "", 
   -1, -1, -1, -1, -1, -1,-1,-1,0, printSettings, debug);
  
  printf("AddTaskMultipleTrackCutIsoConeAnalysis::End configuration\n");
  
  return task;
}


