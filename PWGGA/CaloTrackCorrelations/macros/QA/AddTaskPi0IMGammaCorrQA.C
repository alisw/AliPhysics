/// \file AddTaskPi0IMGammaCorrQA.C
/// \brief Configuration of the analysis QA wagon for PWG-GA EMCal analysis
///
/// Configuration macro of EMCal related PWG-GA  analysis, although it can be
/// also used for PHOS. It does:
/// * a simple photon cluster selection
/// * invariant mass analysis
/// * cluster-charged track correlation analysis (optional)
/// * detector general QA analysis (optional)
/// * track general QA analysis (optional)
///
/// Wagon responsible: Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)


///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle.
/// \param simulation : A bool identifying the data as simulation.
/// \param collision: A string with the colliding system.
/// \param suffix : A string with the type of trigger (default: MB, EMC).
/// \param qaan: execute the detector QA analysis.
/// \param hadronan: execute the track QA and cluster-track correlation analysis.
/// \param calibrate: if OADB was updated with calibration parameters not used in reconstruction, apply them here.
/// \param minTime: minimum time cut, leave it open by default even if calibration available, ns
/// \param maxTime: maximum time cut, leave it open by default even if calibration available, ns
/// \param minCen : An int to select the minimum centrality, -1 means no selection.
/// \param maxCen : An int to select the maximum centrality, -1 means no selection.
/// \param debugLevel : An int to define the debug level of all the tasks.
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskPi0IMGammaCorrQA(const TString  calorimeter   = "EMCAL",
                                                             const Bool_t   simulation    = kFALSE,
                                                             const TString  collision     = "pp",
                                                             const TString  suffix        = "default",
                                                             const Bool_t   qaan          = kFALSE,
                                                             const Bool_t   hadronan      = kFALSE,
                                                             const Bool_t   calibrate     = kTRUE,
                                                             const Int_t    minTime       = -1000,
                                                             const Int_t    maxTime       =  1000,
                                                             const Int_t    minCen        = -1,
                                                             const Int_t    maxCen        = -1,
                                                             const Int_t    debugLevel    = -1,
                                                             const Int_t    year          = 2015
                                                          )
{
  if(simulation && !suffix.Contains("default"))
  {
    printf("AddTaskPi0IMGammaCorrQA - CAREFUL : Triggered events not checked in simulation!! \n");
    return 0x0;
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
  
  // Make sure the B field is enabled for track selection, some cuts need it
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);
  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Name for containers
  
  TString containerName = Form("%s_%s",calorimeter.Data(), suffix.Data());
  
  if(collision!="pp" && maxCen>=0) containerName+=Form("Cen%d_%d",minCen,maxCen);
    
  printf("AddTaskPi0IMGammaCorrQA - Container NAME: %s \n",containerName.Data());
  
  // #### Configure analysis ####
    
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();

  // General frame setting and configuration
  maker->SetReader   ( ConfigureReader   (inputDataType,collision,calibrate,minTime,maxTime,minCen,maxCen,simulation,year,debugLevel) );
  maker->SetCaloUtils( ConfigureCaloUtils(calorimeter,simulation,calibrate,year,debugLevel) );
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important
  
  // Photon analysis
  maker->AddAnalysis(ConfigurePhotonAnalysis(calorimeter,collision,containerName,simulation,year     ,debugLevel), n++); // Photon cluster selection
  maker->AddAnalysis(ConfigurePi0Analysis   (calorimeter,collision,containerName,simulation,year,qaan,debugLevel) ,n++); // Previous photon invariant mass

  if(hadronan)
  {
    maker->GetReader()->SwitchOnCTS();
    maker->AddAnalysis(ConfigureChargedAnalysis(collision,containerName,simulation,year,debugLevel), n++); // charged tracks plots
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",calorimeter,collision,containerName,simulation,year,debugLevel), n++); // Gamma hadron correlation
  }
  
  if(qaan) maker->AddAnalysis(ConfigureQAAnalysis(calorimeter,collision,simulation,year,debugLevel),n++);
  
  maker->SetAnaDebug(debugLevel)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker() ;
  maker->SwitchOffDataControlHistograms();
  if(suffix.Contains("EMC"))
    maker->SwitchOnDataControlHistograms();

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

///
/// Configure the class handling the events and cluster/tracks filtering.
///
AliCaloTrackReader * ConfigureReader(TString inputDataType, TString collision, Bool_t calibrate,
                                     Int_t   minTime,       Int_t maxTime,
                                     Int_t   minCen,        Int_t maxCen,
                                     Bool_t  simulation,    Int_t year,        Int_t debugLevel)
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

  // Time cut 
  reader->SwitchOffUseParametrizedTimeCut();
  reader->SwitchOnUseEMCALTimeCut() ;
  reader->SetEMCALTimeCut(minTime,maxTime);
  
  reader->SwitchOffUseTrackTimeCut();
  reader->SetTrackTimeCut(-1e10,1e10);

  reader->SwitchOffFiducialCut();
  
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
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
    if(year > 2010)
    {
      //Hybrids 2011
      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
      reader->SetTrackCuts(esdTrackCuts);
      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
      reader->SetTrackComplementaryCuts(esdTrackCuts2);
    }
    else
    {
      //Hybrids 2010
      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001006);
      reader->SetTrackCuts(esdTrackCuts);
      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10041006);
      reader->SetTrackComplementaryCuts(esdTrackCuts2);
    }
  }
  else if(inputDataType=="AOD")
  {
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SetTrackStatus(AliVTrack::kITSrefit);    
  }
  
  // Calorimeter
  
  reader->SetEMCALClusterListName("");
  
  if(calibrate && !simulation) reader->SwitchOnClusterRecalculation();
  else                         reader->SwitchOffClusterRecalculation();
  
  reader->SwitchOnEMCALCells();
  reader->SwitchOnEMCAL();

  reader->SwitchOnPHOSCells();
  reader->SwitchOnPHOS();
  
  //-----------------
  // Event selection
  //-----------------
  
  reader->SwitchOnEventTriggerAtSE();
  
  reader->SetZvertexCut(10.);
  reader->SwitchOnPrimaryVertexSelection();  // and besides primary vertex
  reader->SwitchOffPileUpEventRejection();   // remove pileup
  reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
  
  if(collision=="PbPb")
  {
    reader->SetCentralityBin(minCen,maxCen); // Accept all events, if not select range
    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
  }
  
  if(debugLevel > 0) reader->Print("");
  
  return reader;
}

///
/// Configure the class handling the calorimeter clusters specific methods
///
AliCalorimeterUtils* ConfigureCaloUtils(TString calorimeter, Bool_t simulation, Bool_t calibrate,
                                        Int_t year, Int_t debugLevel)
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
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          simulation,
                          bExotic,
                          bNonLin,
                          bEnCalib,
                          bBadMap,
                          bTiCalib);
  recou->SetExoticCellDiffTimeCut(50.);

  if(calorimeter=="PHOS")
  {
    cu->SetNumberOfSuperModulesUsed(3);
  }
  else
  {
    if      (year == 2010) cu->SetNumberOfSuperModulesUsed(4); // EMCAL first year
    else if (year <  2014) cu->SetNumberOfSuperModulesUsed(10);
    else                   cu->SetNumberOfSuperModulesUsed(20);
  }

  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  if(debugLevel > 0) cu->Print("");
  
  return cu;
}

///
/// Configure the task doing the first photon cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection ...
///
AliAnaPhoton* ConfigurePhotonAnalysis(TString calorimeter,   TString collision,
                                      TString containerName, Bool_t simulation, 
                                      Int_t year,            Int_t debugLevel)
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
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter,collision,year); // see method below
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms(7);
  ana->FillNPrimaryHistograms(4);
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0 ) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the 2 cluster invariant mass analysis
///
AliAnaPi0* ConfigurePi0Analysis(TString calorimeter, TString collision,
                                TString containerName, Bool_t simulation, Int_t year,
                                Bool_t qaan, Int_t debugLevel)
{
  AliAnaPi0 *ana = new AliAnaPi0();
  
  ana->SetDebug(debugLevel);//10 for lots of messages
  
  // Input delta AOD settings
  ana->SetInputAODName(Form("Photon%s",containerName.Data()));
  
  // Calorimeter settings
  ana->SetCalorimeter(calorimeter);
  
  // Acceptance plots
  //  ana->SwitchOnFiducialCut(); // Needed to fill acceptance plots with predefined calorimeter acceptances
  //  ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.7, 100, 180) ; 
  //  ana->GetFiducialCut()->DoEMCALFiducialCut(kTRUE);
  
  ana->SwitchOffFiducialCut();
  ana->SwitchOnRealCaloAcceptance();
  
  // Settings for pp collision mixing
  ana->SwitchOnOwnMix(); //Off when mixing done with general mixing frame
  
  // Cuts 
  if(calorimeter=="EMCAL") ana->SetPairTimeCut(70);
    
  ana->SetNPIDBits(1);
  ana->SetNAsymCuts(1); // no asymmetry cut, previous studies showed small effect.
  // In EMCAL assymetry cut prevents combination of assymetric decays which is the main source of pi0 at high E.
    
  if     (collision == "pp"  )
  {
    ana->SetNCentrBin(1);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(1);
    ana->SetNMaxEvMix(100);
    ana->SetMinPt(0.5);
  }
  else if(collision =="PbPb")
  {
    ana->SetNCentrBin(10);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(4);
    ana->SetNMaxEvMix(10);
    ana->SetMinPt(1.5);
  }
  else if(collision =="pPb")
  {
    ana->SetNCentrBin(1);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(4);
    ana->SetNMaxEvMix(100);
    ana->SetMinPt(0.5);
  }

  ana->SwitchOffMultipleCutAnalysis();
  ana->SwitchOnSMCombinations();
  ana->SwitchOffFillAngleHisto();
  ana->SwitchOffFillOriginHisto();
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaPi0_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter,collision,year); // see method below

  if(simulation) ana->SwitchOnDataMC();

  if(debugLevel > 0) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing charged track selection
///
AliAnaChargedParticles* ConfigureChargedAnalysis(TString collision,TString containerName,
                                                 Bool_t simulation, Int_t year, Int_t debugLevel)
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
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),"",collision,year); // see method below
  
  ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
  ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  
  if(debugLevel > 0) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the trigger particle hadron correlation
///
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle,  TString calorimeter,
                                                                    TString collision, TString containerName,
                                                                    Bool_t simulation, Int_t year, Int_t debugLevel)
{
  AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
  ana->SetDebug(debugLevel);
  
  ana->SetTriggerPtRange(5,100);
  ana->SetAssociatedPtRange(0.2,100);
  //ana->SetDeltaPhiCutRange( TMath::Pi()/2,3*TMath::Pi()/2 ); //[90 deg, 270 deg]
  ana->SetDeltaPhiCutRange  (TMath::DegToRad()*120.,TMath::DegToRad()*240.);
  ana->SetUeDeltaPhiCutRange(TMath::DegToRad()*60. ,TMath::DegToRad()*120.);
  ana->SwitchOffFillEtaGapHistograms();

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
  
  //ana->SwitchOnLeadHadronSelection();
  //ana->SetLeadHadronPhiCut(TMath::DegToRad()*100., TMath::DegToRad()*260.);
  //ana->SetLeadHadronPtCut(0.5, 100);

  
  // Mixing with own pool
  ana->SwitchOffOwnMix();
  
  ana->SetNZvertBin(20);
  ana->SwitchOffCorrelationVzBin() ;
  ana->SwitchOffFillHighMultiplicityHistograms();
  
  if(collision=="pp")
  {
    ana->SetNMaxEvMix(100);    
    ana->SwitchOnTrackMultBins();
    ana->SetNTrackMultBin(10); // same as SetNCentrBin(10);
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
  ana->SwitchOffSeveralUECalculation();
  ana->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("Ana%sHadronCorr_",particle.Data()));
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter, collision,year); // see method below
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0) ana->Print("");
 
  return ana;
}

///
/// Configure the task doing standard calorimeter QA
///
AliAnaCalorimeterQA* ConfigureQAAnalysis(TString calorimeter,  TString collision,
                                         Bool_t simulation,    Int_t year,    Int_t debugLevel)
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
  
  ana->AddToHistogramsName("QA_"); //Begining of histograms name
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter, collision,year); // see method below
  
  if(simulation) ana->SwitchOnDataMC();
  
  if(debugLevel > 0) ana->Print("");
  
  return ana;
}


///
/// Configure the selection of MC events
///
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges, TString calorimeter, 
                            TString collision,               Int_t year)
{
  histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
  
  if(calorimeter=="EMCAL")
  {
    if ( year == 2010 )
    {
      histoRanges->SetHistoPhiRangeAndNBins(79*TMath::DegToRad(), 121*TMath::DegToRad(), 42) ;
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
  histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,500);
  //histoRanges->SetHistoTimeRangeAndNBins(-400.,400,400);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistodRRangeAndNBins(0.,0.15,150);//QA

  // QA, electron, charged
  histoRanges->SetHistoPOverERangeAndNBins(0,  2. ,200);
  histoRanges->SetHistodEdxRangeAndNBins  (0.,200.0,200);
  
  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  histoRanges->SetHistoZRangeAndNBins(-350,350,175);
  histoRanges->SetHistoRRangeAndNBins(430,460,30);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  
  // QA, correlation
  if(collision=="PbPb")
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

