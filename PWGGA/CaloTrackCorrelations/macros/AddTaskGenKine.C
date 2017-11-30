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

/// Global variables to be accessed by the different methods.
TString  kCalorimeter = "EMCAL"; /// detector acceptance of trigger particle.
Int_t    kYears       = 2011;    /// year configuration for geometry setting.
Int_t    kDebug       = -1;      /// debug level.
Bool_t   kPrint       = kFALSE;  /// print configuration settings.

///
/// Main method calling all the configuration.
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskGenKine(TString outputfile,
                                                    const Double_t scaleFactor   = -1)
{
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
  maker->AddAnalysis(ConfigureGenKine(), n++); // Photon cluster selection
   
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
  

  
  if(kPrint) maker->Print("");
  
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
        
  if(kPrint)reader->Print("");
  
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
  
  if(kPrint) cu->Print("");
  
  return cu;
  
}


///
/// Configure the task doing analysis at the generator level
///
AliAnaGeneratorKine* ConfigureGenKine()
{
  AliAnaGeneratorKine *ana = new AliAnaGeneratorKine();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // cluster selection cuts
  
  ana->SwitchOffFiducialCut();

    
  ana->AddToHistogramsName("AnaGenKine_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
  if(kPrint)ana->Print("");
  
  return ana;
}

///
/// Set common histograms binning and ranges
///
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges)
{
  histoRanges->SetHistoPtRangeAndNBins(0, 200, 400) ; // Energy and pt histograms
  
  if(kCalorimeter=="EMCAL")
  {
    if(kYears==2010)
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
  else
  {
    histoRanges->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 320*TMath::DegToRad(), 70) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
  }
  
  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 4.9, 500);
  
  // Invariant mass histoRangeslysis
  histoRanges->SetHistoMassRangeAndNBins(0., 1., 200) ;
  histoRanges->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  
  // check if time calibration is on
  histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistodRRangeAndNBins(0.,0.15,150);//QA
  
  // QA, electron, charged
  histoRanges->SetHistoPOverERangeAndNBins(0,10.,100);
  histoRanges->SetHistodEdxRangeAndNBins(0.,200.,200);
  
  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,100);
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  histoRanges->SetHistoNClusterCellRangeAndNBins(0,500,500);
  histoRanges->SetHistoZRangeAndNBins(-350,350,175);
  histoRanges->SetHistoRRangeAndNBins(430,460,30);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
}

