/// \file AddTaskGenKine.C
/// \ingroup CaloTrackCorrMacros
/// \brief Analysis at generator level configuration.
///
/// Example of configuration AliAnaGeneratorKine task
/// of the package CaloTrackCorrelations. Do analysis at generator
/// level only on high-pT photon/pi0/eta and correlate with charged particles
/// jets, partons.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

// Set includes for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

// ROOT
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>

// AliPhysics
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

// CaloTrackCorrelations frame
#include "AliCaloTrackMCReader.h"
#include "AliAnalysisTaskCaloTrackCorrelation.h"
#include "AliAnaCaloTrackCorrMaker.h"
#include "AliIsolationCut.h"
#include "AliAnaGeneratorKine.h"

#endif

/// Global variables to be accessed by the different methods.
TString  kCalorimeter = "EMCAL"; /// Detector acceptance of trigger particle: EMCAL,  DCAL, PHOS, FullCalo
Int_t    kYears       = 2011;    /// Year configuration for geometry setting.
TString  kCollisName  = "pp";    /// A string with the colliding system
Int_t    kDebug       = -1;      /// Debug level.
Bool_t   kPrintConf   = kFALSE;  /// Print configuration settings.
Float_t  kConeSize    =  0.4;    /// A float setting the isolation cone size higher limit
Float_t  kConeSizeMin =  -1;     /// A float setting the isolation cone size lower limit
Float_t  kIsoCut      =  1.5;    /// A float setting the isolation pT threshold (sum of particles in cone or leading particle)
Int_t    kPartInCone  = AliIsolationCut::kNeutralAndCharged; ///  Type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, kOnlyNeutral, kOnlyCharged
Int_t    kIsoMethod   = AliIsolationCut::kPtThresIC;         ///  An int setting the isolation method: AliIsolationCut::kPtThresIC, ...

///
/// Configure the class handling the events and cluster/tracks filtering.
/// Set it off for most of its utilities, no data filtering just access kine MC.
///
AliCaloTrackReader * ConfigureReader()
{
  AliCaloTrackReader * reader = new AliCaloTrackMCReader();
  
  reader->SetDebug(kDebug);//10 for lots of messages
    
//  reader->SetPtHardAndJetPtComparison(kTRUE);
//  reader->SetPtHardAndJetPtFactor(4);
//  
//  reader->SetPtHardAndClusterPtComparison(kTRUE);
//  reader->SetPtHardAndClusterPtFactor(1.);

  reader->SwitchOffWriteDeltaAOD()  ;
  
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

  reader->SwitchOnFiducialCut();
  //reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;

  reader->SwitchOffCTS();
  reader->SwitchOffEMCALCells();  
  reader->SwitchOffEMCAL();
  reader->SwitchOffPHOSCells();  
  reader->SwitchOffPHOS();
  
  reader->SetZvertexCut(10.);               // Open cut
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex

  reader->SwitchOffPileUpEventRejection();  // remove pileup by default
  reader->SwitchOffV0ANDSelection() ;       // and besides v0 AND
        
  if ( kPrintConf ) reader->Print("");
  
  return reader;

}

///
/// Configure the class handling the calorimeter clusters specific methods.
/// Set it off for most of its utilities, no data no calorimeter clusters used, access kine MC.
///
AliCalorimeterUtils* ConfigureCaloUtils()
{
  
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(kDebug);
  
  cu->SwitchOffClusterPlot();
  
  cu->SwitchOffRecalculateClusterTrackMatching(); // Done in clusterization
  
  cu->SwitchOffBadChannelsRemoval() ;
  
  cu->SwitchOffLoadOwnEMCALGeometryMatrices();
  
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  
  cu->SwitchOffRecalibration(); 
  cu->SwitchOffTimeRecalibration();     
  cu->SwitchOffRunDepCorrection(); 
  cu->SwitchOffCorrectClusterLinearity();
  cu->SwitchOffEMCALOADB() ; 
  
  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
  
  if ( kPrintConf ) cu->Print("");
  
  return cu;
  
}

///
/// Set common histograms binning and ranges
///
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges)
{
  histoRanges->SetHistoPtRangeAndNBins(0, 200, 400) ; // Energy and pt histograms
  
  if ( kCalorimeter == "EMCAL" )
  {    
    if ( kYears == 2010 )
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
    }
    else if ( kYears < 2014 )
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 104) ;
    }
    else // Run2
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 189*TMath::DegToRad(), 111) ;
    }
    
    histoRanges->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  }
  else if ( kCalorimeter=="DCAL" )
  {
    histoRanges->SetHistoPhiRangeAndNBins(260*TMath::DegToRad(), 327*TMath::DegToRad(), 67) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  } 
  else if ( kCalorimeter == "PHOS" )
  {
    histoRanges->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 320*TMath::DegToRad(), 70) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
  }
  else if ( kCalorimeter == "CTS" )
  {
    histoRanges->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    histoRanges->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  if ( kCalorimeter == "FullCalo" ) 
  {
    histoRanges->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 327*TMath::DegToRad(), 247) ;
  }
}


///
/// Configure the isolation cuts 
/// used in ConfigureIsolationAnalysis() and ConfigureHadronCorrelationAnalysis()
///
/// \param ic : Pointer to task doing the isolation
///
void ConfigureIsolationCut(AliIsolationCut * ic)
{
  ic->SetDebug(kDebug);
  ic->SetParticleTypeInCone(kPartInCone);
  ic->SetICMethod(kIsoMethod);
  ic->SetPtFraction(0.1);
  ic->SetPtThreshold(0.5); // default, change in next lines
  ic->SetSumPtThreshold(1.0); // default, change in next lines
  
  if ( kConeSize > 0 && kIsoCut > 0 )
  {
    ic->SetConeSize(kConeSize);
    ic->SetMinDistToTrigger(kConeSizeMin);    
    ic->SetPtThresholdMax(10000);
    
    if ( kIsoMethod == AliIsolationCut::kPtThresIC )
    {
      printf("ConfigureIsolationCuts() *** PtThresMin = %1.1f GeV/c *** R = %1.2f *** R min %1.2f\n",kIsoCut,kConeSize,kConeSizeMin);
      ic->SetPtThreshold(kIsoCut);
    }
    
    if ( kIsoMethod == AliIsolationCut::kSumPtIC || 
         kIsoMethod >= AliIsolationCut::kSumBkgSubIC )
    {
      printf("ConfigureIsolationCuts() *** SumPtMin = %1.1f GeV/c *** R = %1.1f *** R min %1.2f\n",kIsoCut,kConeSize,kConeSizeMin);
      ic->SetSumPtThreshold(kIsoCut);
    }
  }
  else
  {
    printf("ConfigureIsolationCuts() *** Careful, use old hardcoded values\n");
    if ( kCollisName == "pp" )
    {
      ic->SetPtThreshold(0.5);
      ic->SetSumPtThreshold(1.0) ;
      ic->SetConeSize(0.4);
      ic->SetMinDistToTrigger(-1);    
    }
    if ( kCollisName == "PbPb" )
    {
      ic->SetPtThreshold(3.);
      ic->SetSumPtThreshold(3.0) ;
      ic->SetConeSize(0.3);
      ic->SetMinDistToTrigger(-1);    
    }
  }
  
  //ic->SwitchOnFillEtaPhiHistograms();
  
//  if ( kAnaCutsString.Contains("FixIsoConeExcess") ) 
//    ic->SwitchOnConeExcessCorrectionHistograms();
//  else                                         
//    ic->SwitchOffConeExcessCorrectionHistograms();
}

///
/// Configure the task filling generated particle kinematics histograms
/// \param makePartonAna : Bool to de/activate parton/jet related analysis
///
AliAnaGeneratorKine* ConfigureGenKineAnalysis ( Bool_t makePartonAna = kTRUE )
{
  AliAnaGeneratorKine *ana = new AliAnaGeneratorKine();
  
  // Trigger detector, acceptance and pT cut
  ana->SetTriggerDetector(kCalorimeter);
  if ( kCalorimeter == "DCAL" ) 
  {
    TString calo = "EMCAL";
    ana->SetTriggerDetector(calo);
  }
  
  ana->SetMinPt(2); // Trigger photon, pi0 minimum pT
  if      ( kCalorimeter == "EMCAL" )
  {
    if      ( kYears > 2014 ) ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.67,  81.2, 185.8) ; //12 SM
    else if ( kYears > 2010 ) ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.67,  81.2, 178.8) ; //10 SM
    else                    ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.67,  81.2, 118.8) ; // 4 SM
  }
  else if ( kCalorimeter == "DCAL"  ) ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.67, 261.2, 325.8) ; 
  else if ( kCalorimeter == "PHOS"  )
  {
    if      ( kYears > 2014 ) ana->GetFiducialCutForTrigger()->SetSimplePHOSFiducialCut (0.125, 250.5, 319.5) ; 
    else                    ana->GetFiducialCutForTrigger()->SetSimplePHOSFiducialCut (0.125, 260.5, 319.5) ;
  }
  
  if ( kCalorimeter == "FullCalo" ) 
  {
    ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.67,  81.2, 325.8) ; 
  }
  
  // Particles associated to trigger or isolation cone acceptance and pT cut
  ana->SetCalorimeter(kCalorimeter);
  if ( kCalorimeter == "DCAL" ||  kCalorimeter == "FullCalo" ) 
  {
    TString calo = "EMCAL";
    ana->SetCalorimeter(calo);
  }
  
  ana->SetMinChargedPt(0.2);
  ana->SetMinNeutralPt(0.5);
  
  if      ( kCalorimeter == "EMCAL" )
  {
    if      ( kYears > 2014 ) ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.67,  81.2, 185.8) ; //12 SM
    else if ( kYears > 2010 ) ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.67,  81.2, 178.8) ; //10 SM
    else                      ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.67,  81.2, 118.8) ; // 4 SM
  }
  else if ( kCalorimeter == "DCAL"  ) ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.67, 261.2, 325.8) ; 
  else if ( kCalorimeter == "PHOS"  )
  {
    if      ( kYears > 2014 ) ana->GetFiducialCut()->SetSimplePHOSFiducialCut (0.125, 250.5, 319.5) ; 
    else                      ana->GetFiducialCut()->SetSimplePHOSFiducialCut (0.125, 260.5, 319.5) ;
  }
  
  if ( kCalorimeter == "FullCalo" ) 
  {
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.67,  81.2, 325.8) ; 
  }
  
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360);
  
  // Isolation paramters
  AliIsolationCut * ic = ana->GetIsolationCut();
  ConfigureIsolationCut(ic);

  ana->SwitchOnPartonAnalysis();
   if ( !makePartonAna ) 
     ana->SwitchOffPartonAnalysis();
  
  ana->AddToHistogramsName("AnaGenKine_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
  return ana;
}


///
/// Main method calling all the configuration.
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// \param makePartonAna : Bool to de/activate parton/jet related analysis
/// \param partInCone : An int setting the type of particles inside the isolation cone: AliIsolationCut::kNeutralAndCharged, AliIsolationCut::kOnlyNeutral, AliIsolationCut::kOnlyCharged
/// \param isoMethod : An int setting the isolation method: AliIsolationCut::kPtThresIC, ...
/// \param cone : A float setting the isolation cone size higher limit
/// \param coneMin : A float setting the isolation cone size lower limit
/// \param isoCut: A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param calorimeter : A string with the detector: EMCAL, DCAL, PHOS, CTS, FullCalo
/// \param year : An int with the data year
/// \param col : A string with the colliding system
/// \param scaloeFactor : double with pT hard bin scale factor
/// \param debug : An int to define the debug level of all the tasks
/// \param printConf : A bool to print configuration settings
/// \param outputfile: string with output file name
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskGenKine
 ( Bool_t makePartonAna, Int_t partInCone, Int_t isoMethod, Float_t cone, Float_t coneMin, Float_t isoCut, 
   TString calorimeter, Int_t year, TString col,  Double_t scaleFactor   = -1,
   Int_t debug = -1, Bool_t printConf = kFALSE, TString outputfile = "")
{
  // Change global variables
  kCalorimeter = calorimeter; 
  kYears       = year;
  kCollisName  = col;
  kDebug       = debug;
  kPrintConf   = printConf;
  kConeSize    = cone;
  kConeSizeMin = coneMin;
  kIsoCut      = isoCut;
  kPartInCone  = partInCone;
  kIsoMethod   = isoMethod;
  
  // Get the pointer to the existing analysis manager via the static access method.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }  
    
  // #### Configure analysis ####
    
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  printf("SCALE FACTOR %e\n",scaleFactor);
  maker->SetScaleFactor(scaleFactor); // for MC, negative (not scaled) by default
  
  // General frame setting and configuration
  maker->SetReader   (ConfigureReader()   ); 
  maker->SetCaloUtils(ConfigureCaloUtils()); 

  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important  
  maker->AddAnalysis(ConfigureGenKineAnalysis(makePartonAna), n++); // Photon cluster selection
   
  maker->SetAnaDebug(kDebug)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOffAODsMaker() ;
 
  // Calculate the cross section weights, apply them to all histograms 
  // and fill xsec and trial histo. Sumw2 must be activated.
  //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionCalculation(); 
  //maker->SwitchOnSumw2Histograms();
  
  // For recent productions where the cross sections and trials are not stored in separate file
  //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionFromEventHeader() ;
  
  // Just fill cross section and trials histograms.
  maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionHistoFill(); 
  
  if ( kPrintConf ) maker->Print("");
  
  // Create task
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation ("GeneratorKineAnalysis");
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(kDebug);
  task->SetBranches("ESD:AliESDRun.,AliESDHeader"); 
  task->SetAnalysisMaker(maker);
  mgr->AddTask(task);
  
  //Create containers

  if(outputfile.Length()==0) outputfile = AliAnalysisManager::GetCommonFileName(); 
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer("GenKine", TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s",outputfile.Data()));
  
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer("Param_GemKine", TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             "AnalysisParameters.root");
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
      
  return task;
}
