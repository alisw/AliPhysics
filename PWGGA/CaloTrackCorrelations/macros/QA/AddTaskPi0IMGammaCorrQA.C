
AliAnalysisTaskCaloTrackCorrelation *AddTaskPi0IMGammaCorrQA(const TString  calorimeter   = "EMCAL",
                                                             const Bool_t   simulation    = kFALSE,
                                                             const TString  collision     = "pp",
                                                             const TString  suffix        = "default",
                                                             const Bool_t   qaan          = kFALSE,
                                                             const Bool_t   hadronan      = kFALSE,
                                                             const Int_t    minCen        = -1,
                                                             const Int_t    maxCen        = -1,
                                                             const Int_t    debugLevel    = -1
                                                          )
{
  // Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
  
  
  if(simulation)
  {
    printf("AddTaskPi0IMGammaCorrQA - CAREFUL : Triggered events not checked in simulation!! \n");
    TString ssuffix = suffix;
    if(!ssuffix.Contains("default")) return;
  }

  // Get the pointer to the existing analysis manager via the static access method.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTaskPi0IMGammaCorrQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  
  if (!mgr->GetInputEventHandler()) 
  {
    ::Error("AddTaskPi0IMGammaCorrQA", "This task requires an input event handler");
    return NULL;
  }
  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Name for containers
  
  TString containerName = Form("%s_%s",calorimeter.Data(), suffix.Data());
  
  if(collision!="pp" && maxCen>=0) containerName+=Form("Cen%d_%d",minCen,maxCen);
    
  printf("AddTaskPi0IMGammaCorrQA - Container NAME: %s \n",containerName.Data());
  
  // #### Configure analysis ####
    
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();

  // General frame setting and configuration
  maker->SetReader   ( ConfigureReader   (inputDataType,minCen,maxCen,simulation,debugLevel) );
  maker->SetCaloUtils( ConfigureCaloUtils(simulation,debugLevel) );
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important
  
  // Photon analysis
  maker->AddAnalysis(ConfigurePhotonAnalysis(calorimeter,          containerName,simulation     ,debugLevel), n++); // Photon cluster selection
  maker->AddAnalysis(ConfigurePi0Analysis   (calorimeter,collision,containerName,simulation,qaan,debugLevel) ,n++); // Previous photon invariant mass

  if(hadronan)
  {
    maker->GetReader()->SwitchOnCTS();
    maker->AddAnalysis(ConfigureChargedAnalysis(containerName,simulation,debugLevel), n++); // charged tracks plots
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",calorimeter,collision,containerName,simulation,debugLevel), n++); // Gamma hadron correlation
  }
  
  if(qaan) maker->AddAnalysis(ConfigureQAAnalysis(calorimeter,simulation,debugLevel),n++);
  
  maker->SetAnaDebug(debugLevel)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker() ;
  if(simulation || !suffix.Contains("EMC"))
    maker->SwitchOffDataControlHistograms();
  else
    maker->SwitchOnDataControlHistograms();

  if(debugLevel > 0) maker->Print("");
  
  // Create task
  
  TString taskName =Form("Pi0IM_GammaTrackCorr_%s",containerName.Data());
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (taskName);
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(debugLevel);
  //task->SetBranches("ESD:AliESDRun.,AliESDHeader");
  //task->SetBranches("AOD:header,tracks,vertices,emcalCells,caloClusters");
  task->SetAnalysisMaker(maker);
  mgr->AddTask(task);
  
  //Create containers
  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(taskName, TList::Class(),
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s:%s",outputfile.Data(),taskName.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",taskName.Data()), TList::Class(),
                                                             AliAnalysisManager::kParamContainer, 
                                                             Form("%s_Parameters.root",taskName.Data()));
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  
  return task;
}

//___________________________________________________________________________
AliCaloTrackReader * ConfigureReader(TString inputDataType,
                                     Int_t minCen,          Int_t maxCen,
                                     Bool_t  simulation,    Int_t debugLevel)
{
  
  AliCaloTrackReader * reader = 0;
  if     (inputDataType=="AOD")
    reader = new AliCaloTrackAODReader();
  else if(inputDataType=="ESD")
    reader = new AliCaloTrackESDReader();
  else 
    printf("AliCaloTrackReader::ConfigureReader() - Data combination not known input Data=%s\n",
           inputDataType.Data());
  
  reader->SetDebug(debugLevel);//10 for lots of messages

  // MC settings
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

  // Time cut off
  reader->SwitchOffUseTrackTimeCut();
  reader->SwitchOffUseParametrizedTimeCut();
  reader->SwitchOffUseEMCALTimeCut() ;
  reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
  reader->SetTrackTimeCut(-1e10,1e10);

  reader->SwitchOnFiducialCut();
  
  // Tracks
  reader->SwitchOffCTS();
  reader->SwitchOffRejectNoTrackEvents();
  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();
  reader->SwitchOffUseTrackDCACut();
  //reader->SetTrackDCACut(0,0.0105);
  //reader->SetTrackDCACut(1,0.035);
  //reader->SetTrackDCACut(2,1.1);
  
  if(inputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");
//    if(year > 2010)
//    {
      //Hybrids 2011
      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
      reader->SetTrackCuts(esdTrackCuts);
      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
      reader->SetTrackComplementaryCuts(esdTrackCuts2);
//    }
//    else
//    {
//      //Hybrids 2010
//      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001006);
//      reader->SetTrackCuts(esdTrackCuts);
//      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10041006);
//      reader->SetTrackComplementaryCuts(esdTrackCuts2);
//    }
  }
  else if(inputDataType=="AOD")
  {
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SetTrackStatus(AliVTrack::kITSrefit);    
  }
  
  // Calorimeter
  
  reader->SetEMCALClusterListName("");
  reader->SwitchOffClusterRecalculation();
  
  reader->SwitchOnEMCALCells();
  reader->SwitchOnEMCAL();

  reader->SwitchOnPHOSCells();
  reader->SwitchOnPHOS();
  
  //-----------------
  // Event selection
  //-----------------
  
  reader->SwitchOnEventTriggerAtSE();
  
  reader->SetZvertexCut(10.);
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
  reader->SwitchOffPileUpEventRejection();  // remove pileup
  reader->SwitchOffV0ANDSelection() ;       // and besides v0 AND
  
  reader->SetCentralityBin(minCen,maxCen); // Accept all events, if not select range
  reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)

  if(debugLevel > 0) reader->Print("");
  
  return reader;
  
}

//__________________________________________________________________________
AliCalorimeterUtils* ConfigureCaloUtils(Bool_t simulation, Int_t debugLevel)
{
  
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(debugLevel);
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  // Search of local maxima in cluster
  cu->SetLocalMaximaCutE(0.1);
  cu->SetLocalMaximaCutEDiff(0.03);

  cu->SwitchOffClusterPlot();
  cu->SwitchOffRecalculateClusterTrackMatching();
  
  cu->SwitchOnBadChannelsRemoval() ;
  
  //EMCAL settings

  if(!simulation)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  cu->SwitchOffRecalibration(); // Check the reader if it is taken into account during filtering
  cu->SwitchOffRunDepCorrection();

  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          simulation,
                          kTRUE,//kExotic,
                          kTRUE,//kNonLinearity,
                          kFALSE,//kCalibE,
                          kTRUE,//kBadMap,
                          kFALSE);//kCalibT
  recou->SetExoticCellDiffTimeCut(50.);

  cu->SwitchOnCorrectClusterLinearity();

  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  if(debugLevel > 0) cu->Print("");
  
  return cu;
  
}

//_______________________________________________________________________________
AliAnaPhoton* ConfigurePhotonAnalysis(TString calorimeter, TString containerName,
                                      Bool_t simulation, Int_t debugLevel)
{
  
  AliAnaPhoton *ana = new AliAnaPhoton();
  ana->SetDebug(debugLevel); //10 for lots of messages
  
  // cluster selection cuts
  
  ana->SwitchOffFiducialCut();

  ana->SetCalorimeter(calorimeter);
  
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
    ana->SetMaxEnergy(1000); 
    ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off 
    // restrict to less than 100 ns when time calibration is on 
    ana->SetMinDistanceToBadChannel(2, 4, 6); 
    // Not useful if M02 cut is already strong
    ana->SetNLMCut(1, 2) ;
  }
  
  ana->SwitchOnTrackMatchRejection() ;
  ana->SwitchOffTMHistoFill() ;

  
  //PID cuts (shower shape)
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
  //Not used in bayesian
  
  //EMCAL
  caloPID->SetEMCALLambda0CutMax(0.4); // Rather open
  caloPID->SetEMCALLambda0CutMin(0.10);
  
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
    
  //PHOS
  caloPID->SetPHOSDispersionCut(2.5);
  caloPID->SetPHOSRCut(2.);
  
  ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection
  //if(!simulation)ana->SwitchOnFillPileUpHistograms();

  // Input / output delta AOD settings
  ana->SetOutputAODName(Form("Photon%s",containerName.Data()));
  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  ana->SetInputAODName (Form("Photon%s",containerName.Data()));

  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaPhoton_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms(7);
  ana->FillNPrimaryHistograms(4);
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0 ) ana->Print("");
  
  return ana;
  
}

//_________________________________________________________________________
AliAnaPi0* ConfigurePi0Analysis(TString calorimeter, TString collision,
                                TString containerName, Bool_t simulation,
                                Bool_t qaan, Int_t debugLevel)
{
  
  AliAnaPi0 *ana = new AliAnaPi0();
  
  ana->SetDebug(debugLevel);//10 for lots of messages
  
  // Input delta AOD settings
  ana->SetInputAODName(Form("Photon%s",containerName.Data()));
  
  // Calorimeter settings
  ana->SetCalorimeter(calorimeter);
  if(calorimeter=="PHOS") ana->SetNumberOfModules(3); //PHOS first year
  else 
  {                   
//    if     (year == 2010) ana->SetNumberOfModules( 4); // EMCAL first year
//    else if(year == 2011) ana->SetNumberOfModules(10); // Second year
//    else                    ana->SetNumberOfModules(12);
    ana->SetNumberOfModules(12);
  }
  
  //settings for pp collision mixing
  ana->SwitchOnOwnMix(); //Off when mixing done with general mixing frame
  
  // Cuts 
  if(calorimeter=="EMCAL") ana->SetPairTimeCut(70);
  
  if     (collision == "pp"  )
  {
    ana->SetNCentrBin(1);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(1);
    ana->SetNMaxEvMix(100);    
  }
  else if(collision =="PbPb")
  {
    ana->SetNCentrBin(10);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(4);
    ana->SetNMaxEvMix(10);
  }
  else if(collision =="pPb")
  {
    ana->SetNCentrBin(1);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(4);
    ana->SetNMaxEvMix(100);
  }

  ana->SwitchOffMultipleCutAnalysis();
  ana->SwitchOnSMCombinations();
  ana->SwitchOffFillAngleHisto();
  ana->SwitchOffFillOriginHisto();

  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaPi0_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below

  if(simulation) ana->SwitchOnDataMC();

  if(debugLevel > 0) ana->Print("");
  
  return ana;
  
}

//___________________________________________________________________________________
AliAnaChargedParticles* ConfigureChargedAnalysis(TString containerName,
                                                 Bool_t simulation, Int_t debugLevel)
{
  
  AliAnaChargedParticles *ana = new AliAnaChargedParticles();
  ana->SetDebug(debugLevel); //10 for lots of messages
  
  // selection cuts
  
  ana->SetMinPt(0.5);
  ana->SwitchOnFiducialCut();
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ; //more restrictive cut in reader and after in isolation
  
  ana->SwitchOffFillVertexBC0Histograms() ;
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  // Input / output delta AOD settings
  
  ana->SetOutputAODName(Form("Hadron%s",containerName.Data()));
  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  ana->SetInputAODName(Form("Hadron%s",containerName.Data()));

  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaHadrons_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),""); // see method below
  
  ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
  ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  
  if(debugLevel > 0) ana->Print("");
  
  return ana;
  
}

//__________________________________________________________________________________________________________
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle,  TString calorimeter,
                                                                    TString collision, TString containerName,
                                                                    Bool_t simulation, Int_t debugLevel)
{
  
  AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
  ana->SetDebug(debugLevel);
  
  ana->SetMinimumTriggerPt(5);
  ana->SetAssociatedPtRange(0.2,100); 
  ana->SetDeltaPhiCutRange( TMath::Pi()/2,3*TMath::Pi()/2 ); //[90 deg, 270 deg]
  
  ana->SetNAssocPtBins(4);
  
  ana->SetAssocPtBinLimit(0, 0.5) ;
  ana->SetAssocPtBinLimit(1, 2) ;
  ana->SetAssocPtBinLimit(2, 5) ;
  ana->SetAssocPtBinLimit(3, 10) ;
  ana->SetAssocPtBinLimit(4, 20) ;

  ana->SelectIsolated(kFALSE); // do correlation with isolated photons

  //if(!simulation) ana->SwitchOnFillPileUpHistograms();
  
  ana->SwitchOnAbsoluteLeading();  // Select trigger leading particle of all the selected tracks
  ana->SwitchOffNearSideLeading(); // Select trigger leading particle of all the particles at +-90 degrees, default
  
  // Mixing with own pool
  ana->SwitchOffOwnMix();
  
  ana->SetNZvertBin(20);
  
  if(collision=="pp")
  {
    ana->SetNMaxEvMix(100);    
    ana->SwitchOnTrackMultBins();
    ana->SetNCentrBin(9); // Fixed track mult values
    ana->SetNRPBin(1);
  }
  else 
  {
    ana->SetNMaxEvMix(10);    
    ana->SwitchOffTrackMultBins(); // centrality bins
    ana->SetNCentrBin(10); 
    ana->SetNRPBin(3);
  }
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s",particle.Data(),containerName.Data()));
  ana->SetAODObjArrayName(Form("%sHadronCorr_%s",particle.Data(),containerName.Data()));
  
  ana->SwitchOffPi0TriggerDecayCorr();
  ana->SwitchOffDecayTriggerDecayCorr();
  ana->SwitchOffNeutralCorr(); // Do only correlation with TPC
  ana->SwitchOffHMPIDCorrelation();
  ana->SwitchOffFillBradHistograms();
  
  // Underlying event
  ana->SwitchOffEventSelection();
  ana->SwitchOffSeveralUECalculation();
  ana->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  ana->SetMultiBin(1);
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("Ana%sHadronCorr_",particle.Data()));
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0) ana->Print("");
 
  return ana;
  
}

//________________________________________________________________________________
AliAnaCalorimeterQA* ConfigureQAAnalysis(TString calorimeter,  Bool_t simulation,
                                         Int_t debugLevel)
{
  
  AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
  ana->SetDebug(debugLevel); //10 for lots of messages
  ana->SetCalorimeter(calorimeter);
  
  ana->SetTimeCut(-1e10,1e10); // Open time cut
  ana->SwitchOnCorrelation();
  ana->SwitchOffStudyBadClusters() ;
  ana->SwitchOffFiducialCut();
  ana->SwitchOffFillAllTH3Histogram();
  ana->SwitchOffFillAllPositionHistogram();
  ana->SwitchOffFillAllPositionHistogram2();
  ana->SwitchOffStudyBadClusters();
  ana->SwitchOffStudyClustersAsymmetry();
  ana->SwitchOffStudyWeight();
  ana->SwitchOnFillAllTrackMatchingHistogram();
  ana->SwitchOnFillAllCellTimeHisto() ;
  
  if(calorimeter=="EMCAL")
  {
//    if     (year==2010)  ana->SetNumberOfModules(4); 
//    else if(year==2011)  ana->SetNumberOfModules(10);
//    else
    ana->SetNumberOfModules(12);
  }
  else 
  {//PHOS
    ana->SetNumberOfModules(3); 
  }
  
  ana->AddToHistogramsName("QA_"); //Begining of histograms name
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0) ana->Print("");
  
  return ana;
  
}


//________________________________________________________
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges,
                            TString calorimeter)
{
  // Set common bins for all analysis and MC histograms filling
    
  histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
  
  if(calorimeter=="EMCAL")
  {
//    if(year==2010)
//    {
//      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
//      histoRanges->SetHistoXRangeAndNBins(-230,90,120); // QA
//      histoRanges->SetHistoYRangeAndNBins(370,450,40);  // QA
//    }
//    else if(year==2011)
//    {           
//      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
//      histoRanges->SetHistoXRangeAndNBins(-600,90,200); // QA
//      histoRanges->SetHistoYRangeAndNBins(100,450,100); // QA
//    }
//    else
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
  histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,2000);
  //histoRanges->SetHistoTimeRangeAndNBins(-400.,400,400);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistodRRangeAndNBins(0.,0.15,150);//QA

  // QA, electron, charged
  histoRanges->SetHistoPOverERangeAndNBins(0,  2.5 ,500);
  histoRanges->SetHistodEdxRangeAndNBins  (0.,250.0,500);
  
  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,100);
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  histoRanges->SetHistoNClusterCellRangeAndNBins(0,50,50);
  histoRanges->SetHistoNClustersRangeAndNBins(0,100,100);
  histoRanges->SetHistoZRangeAndNBins(-400,400,200);
  histoRanges->SetHistoRRangeAndNBins(400,450,25);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
  
}

