/// \file AddTaskCalorimeterQA.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Configuration of EMCal QA
///
/// Configuration macro of EMCal detector QA analysis, although it can be
/// also used for PHOS.
///
/// Wagon contacts: EMCAL Gustavo.Conesa.Balbastre@cern.ch, Marie.Germain@cern.ch
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///


///
/// Main method calling all the configuration
/// It creates a CaloTrackCorrelations task for calorimeters performance studies,
/// configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param suffix : A string with the type of trigger (default: MB, EMC)
/// \param simulation : A bool identifying the data as simulation
/// \param outputFile : A string to change the name of the histograms output file, default is AnalysisResults.root
/// \param year: The year the data was taken, used to configure some histograms
/// \param printSettings : A bool to enable the print of the settings per task
/// \param calibrate: if OADB was updated with calibration parameters not used in reconstruction, apply them here.
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskCalorimeterQA(const char *suffix="default",
                                                          Bool_t  simulation = kFALSE,
                                                          TString outputFile = "",
                                                          Int_t   year = 2015, 
                                                          Bool_t  printSettings = kFALSE,
                                                          Bool_t  calibrate = kTRUE)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPartCorr", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPartCorr", "This task requires an input event handler");
    return NULL;
  }
  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  Bool_t kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
  
  TString ssuffix = suffix;
  if(kUseKinematics || simulation)
  {
    simulation = kTRUE;
    printf("AddTaskCalorimeterQA - CAREFUL : Triggered events not checked in simulation!! \n");
    if(!ssuffix.Contains("default")) return;
  }
  
  // Configure analysis
  //===========================================================================
  
  //Reader
  //For this particular analysis few things done by the reader.
  //Nothing else needs to be set.
  
  AliCaloTrackReader * reader = 0x0;
  if     (inputDataType.Contains("AOD")) reader = new AliCaloTrackAODReader();
  else if(inputDataType.Contains("ESD")) reader = new AliCaloTrackESDReader();
  //reader->SetDebug(10);//10 for lots of messages
  
  reader->SwitchOnEMCALCells(); 
  reader->SwitchOnEMCAL();
  reader->SwitchOnPHOSCells(); // For correlation plots
  reader->SwitchOnPHOS();      // For correlation plots
  reader->SetEMCALPtMin(0.); 
  reader->SwitchOnCTS();
  reader->SetCTSPtMin  (0.);
  reader->SetZvertexCut(10.);
  
  if(simulation)
  {
    if(inputDataType == "ESD")
    {
      reader->SwitchOnStack();          
      reader->SwitchOffAODMCParticles(); 
    }
    else if(inputDataType == "AOD")
    {
      reader->SwitchOffStack();          
      reader->SwitchOnAODMCParticles(); 
    }
  }
  
  reader->SetDeltaAODFileName(""); //Do not create deltaAOD file, this analysis do not create branches.
  reader->SwitchOffWriteDeltaAOD()  ;
  
  if(!ssuffix.Contains("default"))
  {
    reader->SwitchOnTriggerPatchMatching();
    reader->SwitchOffBadTriggerEventsRemoval();
    reader->SetTriggerPatchTimeWindow(8,9);
    //reader->SetEventTriggerL0Threshold(2.);
  }
  
  if(!simulation) reader->AnalyzeOnlyPhysicsEvents(); // in case physics selection was not on
  
  if(calibrate && !simulation) reader->SwitchOnClusterRecalculation();
  else                         reader->SwitchOffClusterRecalculation();
  
  if(printSettings) reader->Print("");
  
  // *** Calorimeters Utils	***
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);

  if      (year == 2010) cu->SetNumberOfSuperModulesUsed(4); //EMCAL first year
  else if (year <  2014) cu->SetNumberOfSuperModulesUsed(10);
  else                   cu->SetNumberOfSuperModulesUsed(20);
  
  
  cu->SwitchOffRecalibration(); // Check the reader if it is taken into account during filtering
  cu->SwitchOffRunDepCorrection();
  
  cu->SwitchOnCorrectClusterLinearity();
  
  Bool_t bExotic  = kTRUE;
  Bool_t bNonLin  = kTRUE;
  Bool_t bBadMap  = kTRUE;
  
  Bool_t bEnCalib = kFALSE;
  Bool_t bTiCalib = kFALSE;
  
  if(calibrate && !simulation)
  {
    cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
    cu->SwitchOffRunDepCorrection();    
    cu->SwitchOnRecalculateClusterPosition() ;
    
    bEnCalib = kTRUE;
    bTiCalib = kTRUE;
  }
  
  AliEMCALRecoUtils* recou = cu->GetEMCALRecoUtils();

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          simulation,
                          bExotic,
                          bNonLin,
                          bEnCalib,
                          bBadMap,
                          bTiCalib);  
  
  if(bBadMap) 
    cu->SwitchOnBadChannelsRemoval();
  
  cu->SetDebug(-1);
  if(printSettings) cu->Print("");	
  
  // ##### Analysis algorithm settings ####
  
  AliAnaCalorimeterQA *emcalQA = new AliAnaCalorimeterQA();
  //emcalQA->SetDebug(10); //10 for lots of messages
  emcalQA->SetCalorimeter("EMCAL");
  if(simulation)
  {
    // Access MC stack and fill more histograms
    emcalQA->SwitchOnDataMC() ;
    
    emcalQA->SwitchOffStudyBadClusters();
    emcalQA->SwitchOffFillAllCellTimeHisto();
  }
  else
  {
    emcalQA->SwitchOffDataMC() ;
    emcalQA->SwitchOffStudyBadClusters();
    emcalQA->SwitchOnFillAllCellTimeHisto();
  }
  
  emcalQA->AddToHistogramsName("EMCAL_"); //Begining of histograms name
  emcalQA->SwitchOffFiducialCut();
  emcalQA->SwitchOnCorrelation();
  emcalQA->SwitchOffFillAllTH3Histogram();
  emcalQA->SwitchOffFillAllPositionHistogram();
  emcalQA->SwitchOffFillAllPositionHistogram2();
  
  //Set Histrograms bins and ranges
  emcalQA->GetHistogramRanges()->SetHistoPtRangeAndNBins(0, 50, 100) ;
  emcalQA->GetHistogramRanges()->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  emcalQA->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-0.70, 0.70, 140) ;
  
  if     ( year==2010 )
  {  
    emcalQA->GetHistogramRanges()->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 120*TMath::DegToRad(), 48) ;
    emcalQA->GetHistogramRanges()->SetHistoXRangeAndNBins(-230,90,120);
    emcalQA->GetHistogramRanges()->SetHistoYRangeAndNBins(370,450,40);
  }
  else if ( year < 2014 )
  {            
    emcalQA->GetHistogramRanges()->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 180*TMath::DegToRad(), 120) ;
    emcalQA->GetHistogramRanges()->SetHistoXRangeAndNBins(-460,90,200);
    emcalQA->GetHistogramRanges()->SetHistoYRangeAndNBins(100,450,100);
  }
  else // Run2
  {
    emcalQA->GetHistogramRanges()->SetHistoPhiRangeAndNBins(80*TMath::DegToRad(), 327*TMath::DegToRad(), 250) ; 
    emcalQA->GetHistogramRanges()->SetHistoXRangeAndNBins(-460,460,230);
    emcalQA->GetHistogramRanges()->SetHistoYRangeAndNBins(-450,450,225);          
  }
  
  emcalQA->GetHistogramRanges()->SetHistoMassRangeAndNBins(0., 0.65, 325) ;
  emcalQA->GetHistogramRanges()->SetHistoAsymmetryRangeAndNBins(0., 1. , 10 );
  emcalQA->GetHistogramRanges()->SetHistoPOverERangeAndNBins(0,2.,50);
  emcalQA->GetHistogramRanges()->SetHistodEdxRangeAndNBins(0.,200.,100);
  emcalQA->GetHistogramRanges()->SetHistodRRangeAndNBins(0.,0.10,50);
  //emcalQA->GetHistogramRanges()->SetHistoTimeRangeAndNBins( 400,900,250);
  //emcalQA->GetHistogramRanges()->SetHistoTimeRangeAndNBins(-275,275,250);
  emcalQA->GetHistogramRanges()->SetHistoTimeRangeAndNBins(-275,975,250);  
  emcalQA->GetHistogramRanges()->SetHistoRatioRangeAndNBins(0.,2.,100);
  emcalQA->GetHistogramRanges()->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  emcalQA->GetHistogramRanges()->SetHistoNClusterCellRangeAndNBins(0,50,50);
  emcalQA->GetHistogramRanges()->SetHistoZRangeAndNBins(-350,350,175);
  emcalQA->GetHistogramRanges()->SetHistoRRangeAndNBins(430,460,30);
  emcalQA->GetHistogramRanges()->SetHistoV0SignalRangeAndNBins(0,5000,100);
  emcalQA->GetHistogramRanges()->SetHistoV0MultiplicityRangeAndNBins(0,5000,100);
  emcalQA->GetHistogramRanges()->SetHistoTrackMultiplicityRangeAndNBins(0,2500,100);
  emcalQA->GetHistogramRanges()->SetHistoShowerShapeRangeAndNBins(0, 3, 120);
  emcalQA->GetHistogramRanges()->SetHistoDiffTimeRangeAndNBins(-300, 300, 120);
  emcalQA->GetHistogramRanges()->SetHistoTrackResidualEtaRangeAndNBins(-0.075,0.075,50);
  emcalQA->GetHistogramRanges()->SetHistoTrackResidualPhiRangeAndNBins(-0.075,0.075,50);

  if(printSettings) emcalQA->Print("");
  
  // #### Configure Maker ####
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  if(ssuffix.Contains("default")) maker->SwitchOffDataControlHistograms();
  else                            maker->SwitchOnDataControlHistograms();
  maker->SetReader(reader);//pointer to reader
  maker->SetCaloUtils(cu); //pointer to calorimeter utils
  maker->AddAnalysis(emcalQA,0);
  maker->SetAnaDebug(-1)  ; // 0 to at least print the event number
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOffAODsMaker()  ;
  
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
  
  printf("======================== \n");
  printf(" End Configuration of Calorimeter QA \n");
  printf("======================== \n");
  
  // Create task
  //===========================================================================
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("CalorimeterPerformance_%s",suffix));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  task->SetAnalysisMaker(maker);	
  task->SetBranches("ESD:AliESDRun.,AliESDHeader"); //just a trick to get Constantin's analysis to work
  mgr->AddTask(task);
  
  //Create containers
  //  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer("Calo.Performance",TList::Class(),
  //							   AliAnalysisManager::kOutputContainer, "Calo.Performance.root");
	
  TString cname;
  if(outputFile.Length()==0)outputFile = AliAnalysisManager::GetCommonFileName(); 
  
  
  cname = Form("CaloQA_%s", suffix);
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(cname, TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s:%s",outputFile.Data(),cname.Data()));
  
//  cname = Form("CaloQACuts_%s", suffix);
//  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(cname, TList::Class(), 
//                                                             AliAnalysisManager::kParamContainer, 
//                                                             Form("%s:%s",outputFile.Data(),cname.Data()));
  
	//Form("%s:PartCorrCuts",outputfile.Data()));	
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
//  mgr->ConnectOutput (task, 2, cout_cuts);
  
  return task;
}


