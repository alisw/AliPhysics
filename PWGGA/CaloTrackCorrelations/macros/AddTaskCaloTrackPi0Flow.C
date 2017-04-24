
/// \file AddTaskPi0Flow.C
/// \ingroup CaloTrackCorrMacros
/// \brief Configuration of pi0 flow analysis with calorimeters
///
/// Configuration macro for analysis pi0 flow.
///
/// \author Qyie Shou, <qiye.shou@cern.ch> Wuhan 

// Global variables to be accessed by the different methods
TString kData          = "";                                    ///< Declare data MC or deltaAOD
TString kInputDataType = "AOD";                                 ///< Declare data ESD/AOD
TString kCalorimeter   = "EMCAL";                               ///< Use main analysis detector EMCal or PHOS or CTS
Bool_t  kSimulation    = kFALSE;                                ///< Declare the analysis simulation
Bool_t  kEventSelection= kFALSE;                                ///< Remove bad events
Bool_t  kExotic        = kTRUE;                                 ///< Remove exotic clusters
Bool_t  kNonLinearity  = kFALSE;                                ///< Correct cluster non linearity
Int_t   kYears         = 2011;                                  ///< Declare the year of the data
TString kCollisions    = "PbPb";                                ///< Declare the collision type of the data
TString kClusterArray  = "V1_Ecell150_Eseed300_DT0_WT0";        ///< Name of branch with clusters, from AliAnalysisTaskEMCALClusterize
Bool_t  kRecalTM       = kFALSE;                                ///< Recalculate track-cluster matching
Bool_t  kTM            = kTRUE;                                 ///< Remove matched clusters to tracks
Int_t   kMinCen        = -1;                                    ///< Set the minimum centrality to be analyzed
Int_t   kMaxCen        = -1;                                    ///< Set the maximum centrality to be analyzed
Bool_t  kCalibE        = kTRUE;                                 ///< Calibrate energy of clusters
Bool_t  kBadMap        = kTRUE;                                 ///< Reject bad cells/clusters
Bool_t  kCalibT        = kTRUE;                                 ///< Calibrate time of clusters
Bool_t  kTender        = kFALSE;                                ///< Declare that tender was executed
Bool_t  kOutputAOD     = kFALSE;                                ///< Create output AOD with generated particle AOD objects
Bool_t  kPrint         = kFALSE;                                ///< Print setted parameters when configuring
Int_t   kRunNumber     = -1;                                    //< Declare the run number
Bool_t  kPhosCali      = kTRUE;                                 ////< Switch on EP flattening by phos 
Bool_t  kCentFlat      = kFALSE;                                ///< Switch on Centrality flattening
Int_t   kDebug         = 0;                                     ///< Do the analysis with this debug level

Bool_t  kUseKinematics = kFALSE;                                ///< Use the MC information
TString kName          = "";                                    ///< Name of the analysis, used in created AOD branches and histo container

///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskCaloTrackPi0Flow(const TString  data          = "",
                                                             const TString  dataType      = "AOD",
                                                             const TString  calorimeter   = "EMCAL", 
                                                             const Bool_t   simulation    = kFALSE,
                                                             const Bool_t   eventsel      = kFALSE,
                                                             const Bool_t   exotic        = kTRUE,
                                                             const Bool_t   nonlin        = kFALSE,
                                                             const Int_t    year          = 2011,
                                                             const TString  col           = "PbPb", 
                                                             AliVEvent::EOfflineTriggerTypes trig = AliVEvent::kCentral + AliVEvent::kSemiCentral + AliVEvent::kMB + AliVEvent::kEMCEGA,
                                                             const TString  clustersArray = "V1_Ecell150_Eseed300_DT0_WT0",
                                                             const Int_t    nlmMin        = 1,
                                                             const Int_t    nlmMax        = 2,
                                                             const Bool_t   simpleM02Cut  = kFALSE,
                                                             const Bool_t   simpleMassCut = kFALSE,
                                                             const Double_t massPi0Min    = 0.11,
                                                             const Double_t massPi0Max    = 0.18,
                                                             const Bool_t   recaltm       = kFALSE,
                                                             const Bool_t   tm            = kTRUE,
                                                             const Int_t    minCen        = -1,
                                                             const Int_t    maxCen        = -1,
                                                             const Bool_t   calibE        = kTRUE,
                                                             const Bool_t   badmap        = kTRUE,
                                                             const Bool_t   calibT        = kTRUE,
                                                             const Bool_t   tender        = kFALSE,
                                                             const Bool_t   outputAOD     = kFALSE, 
                                                             const Bool_t   printSettings = kFALSE,
                                                             const Int_t    runNumber     = -1,
                                                             const Bool_t   isPhosCali    = kTRUE, 
                                                             const Bool_t   isCentFlat    = kFALSE,
                                                             const Int_t    debug         = 0
                                                             )
{
  kData           = data;
  kInputDataType  = dataType;
  kCalorimeter    = calorimeter;
  kSimulation     = simulation;
  kEventSelection = eventsel;
  kExotic         = exotic;
  kNonLinearity   = nonlin;
  kYears          = year;
  kCollisions     = col;
  kClusterArray   = clustersArray;
  kRecalTM        = recaltm;
  kTM             = tm;
  kMinCen         = minCen;
  kMaxCen         = maxCen;
  kCalibE         = calibE;
  kBadMap         = badmap;
  kCalibT         = calibT;
  kTender         = tender;
  kOutputAOD      = outputAOD;
  kPrint          = printSettings;
  kRunNumber      = runNumber;
  kPhosCali       = isPhosCali;
  kCentFlat       = isCentFlat;
  kDebug          = debug;

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
  kName = Form("%s_Cl%s_TM%d",kCalorimeter.Data(), kClusterArray.Data(), kTM);
  if (kCollisions=="PbPb" && kMaxCen>=0) kName+=Form("Cen%d_%d",kMinCen,kMaxCen);
  printf("<<<< NAME: %s >>>>>\n",kName.Data());
  
  //
  // Create maker
  //    
  AliAnaCaloTrackCorrMaker* maker = new AliAnaCaloTrackCorrMaker();
  // General frame setting and configuration
  maker->SetReader   (ConfigureReader()   ); 
  maker->SetCaloUtils(ConfigureCaloUtils()); 
  
  //
  // Add and configure analysis
  //
  Int_t n = 0; // Analysis number, order is important
  
  // Split cluster analysis
  if (kCalorimeter == "EMCAL") {
    // Pi0 event by event selection, cluster splitting
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kSSCalo, kTRUE, kTRUE, nlmMin, nlmMax, 
                                               simpleM02Cut, simpleMassCut, massPi0Min, massPi0Max), n++);
  }
  maker->AddAnalysis(ConfigurePi0Flow(), n++);
  
  maker->SetAnaDebug(kDebug)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(kData.Contains("delta")) maker->SwitchOffAODsMaker() ;
  else                        maker->SwitchOnAODsMaker()  ;
  if(kPrint) maker->Print("");
  if(kSimulation) maker->SwitchOffDataControlHistograms();
  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, kCalorimeter.Data());
 
  //
  // Create task
  //
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation(Form("CaloTrackCorr%s",kName.Data()));
  task->SetConfigFileName(""); // Don't configure the analysis via configuration file.
  task->SetDebugLevel(kDebug);
  // task->SetBranches("ESD:AliESDRun.,AliESDHeader");
  task->SelectCollisionCandidates(trig);
  task->SetAnalysisMaker(maker);
  mgr->AddTask(task);

  //
  // Create and set containers
  //
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(kName, TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s",outputfile.Data()));
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",kName.Data()), TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             "AnalysisParameters.root");
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);

  return task;
}

///
/// Configure the class handling the events and cluster/tracks filtering.
///
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
        printf("Set time parameters for LHC11c");
        reader->SetEMCALParametrizedMinTimeCut(0,-5  ); reader->SetEMCALParametrizedMinTimeCut(1,-1 ); reader->SetEMCALParametrizedMinTimeCut(2, 1.87); reader->SetEMCALParametrizedMinTimeCut(3, 0.4);   
        reader->SetEMCALParametrizedMaxTimeCut(0, 3.5); reader->SetEMCALParametrizedMaxTimeCut(1, 50); reader->SetEMCALParametrizedMaxTimeCut(2, 0.15); reader->SetEMCALParametrizedMaxTimeCut(3, 1.6);   
      }
      else if(kRunNumber >= 156447 && kRunNumber <= 159635 )
      {
        printf("Set time parameters for LHC11d");
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
  
  reader->SwitchOnFiducialCut();
  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;

  // Tracks
  reader->SwitchOnCTS();
  reader->SwitchOffRejectNoTrackEvents();

  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();
  
  reader->SwitchOffUseTrackTimeCut();
  reader->SetTrackTimeCut(0,50);
  
  reader->SwitchOffUseTrackDCACut();
  //reader->SetTrackDCACut(0,0.0105);
  //reader->SetTrackDCACut(1,0.035);
  //reader->SetTrackDCACut(2,1.1);
  
  if(kInputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
    //AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
    //reader->SetTrackCuts(esdTrackCuts);
    //reader->SwitchOnConstrainTrackToVertex();
    
    if(kYears>2010)
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
  else if(kInputDataType=="AOD")
  {
    //reader->SetTrackFilterMask(128);           // Filter bit, not mask, use if off hybrid
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SetTrackStatus(AliVTrack::kITSrefit);
    //reader->SwitchOnTrackHitSPDSelection();    // Check that the track has at least a hit on the SPD, not much sense to use for hybrid or TPC only tracks
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
  
  if(!kNonLinearity) reader->SwitchOffClusterELinearityCorrection();
  
  if(kCalorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();  
    reader->SwitchOnEMCAL();
  }
  if(kCalorimeter == "PHOS") { 
    reader->SwitchOnPHOSCells();  
    reader->SwitchOnPHOS();
  }
  
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
  
  //reader->RejectFastClusterEvents()  ;

  // Event triggered selection settings
  reader->SwitchOffTriggerPatchMatching();
  //reader->SwitchOffTriggerClusterTimeRecal() ;

  //redefine for other periods, triggers

  //if(!kUseKinematics) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
  
  reader->SwitchOnEventTriggerAtSE();  
  
  reader->SetZvertexCut(10.);               // Open cut
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex

  if(kEventSelection)
  {
    reader->SwitchOnEventPileUpRejection(); // remove pileup by default
    reader->SwitchOnV0ANDSelection() ;      // and besides v0 AND
  }
  else 
  {
    reader->SwitchOffPileUpEventRejection();// remove pileup by default
    reader->SwitchOffV0ANDSelection() ;     // and besides v0 AND
  }
    
  if(kCollisions=="PbPb") 
  {
    // Centrality
    reader->SwitchOnAliCentrality();
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
    reader->SetCentralityBin(kMinCen,kMaxCen); // Accept all events, if not select range
    
    // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
    reader->SetEventPlaneMethod("V0");
  }
  
  if(kPrint) reader->Print("");
  
  return reader;
}

///
/// Configure the class handling the calorimeter clusters specific methods
///
AliCalorimeterUtils* ConfigureCaloUtils()
{
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(kDebug);
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
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

  if (kRecalTM) cu->SwitchOnRecalculateClusterTrackMatching(); // Done in clusterization
  else          cu->SwitchOffRecalculateClusterTrackMatching();
  
  cu->SwitchOnBadChannelsRemoval() ;
  
  //EMCAL settings

  if(!kSimulation)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  
  if(!kSimulation)
  {
    cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
    if(kClusterArray == "" && !kTender) cu->SwitchOnRunDepCorrection();
  }
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          kSimulation,                             
                          kExotic,
                          kTRUE,//kNonLinearity,
                          kCalibE, 
                          kBadMap,
                          kCalibT);
  //recou->SetExoticCellDiffTimeCut(50.);

  
  if( kNonLinearity ) 
  { 
    printf("ConfigureCaloUtils() - Apply non linearity to EMCAL\n");
    cu->SwitchOnCorrectClusterLinearity();
  }
    
  printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
  printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
  
  
  if(kCalorimeter=="PHOS")
  {
    if (kYears <  2014) cu->SetNumberOfSuperModulesUsed(3);
    else                cu->SetNumberOfSuperModulesUsed(4);
  }
  else
  {
    if      (kYears == 2010) cu->SetNumberOfSuperModulesUsed(4); //EMCAL first year
    else if (kYears <  2014) cu->SetNumberOfSuperModulesUsed(10);
    else                     cu->SetNumberOfSuperModulesUsed(20);
  }
  
  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  if(kPrint) cu->Print("");
  
  return cu;
}

///
/// Configure the task doing the pi0 even by event selection via the split method.
/// Here the pairs, clusters, are added to an AOD branch to be used by other analysis
/// unlike in ConfigurePi0Analysis.
///
AliAnaPi0EbE* ConfigurePi0EbEAnalysis(TString particle,
                                      Int_t analysis, Bool_t useSS = kTRUE, Bool_t useAsy = kTRUE,
                                      Int_t nlmMin = 1, Int_t nlmMax = 2, 
                                      Bool_t simpleSplitM02Cut = kFALSE, Bool_t simpleSplitMassCut = kFALSE,
                                      Double_t massPi0Min, Double_t massPi0Max)
{
  AliAnaPi0EbE *ana = new AliAnaPi0EbE();
  ana->SetDebug(kDebug);//10 for lots of messages
  
  ana->SetAnalysisType(analysis);
  TString opt = "";
  if (analysis == AliAnaPi0EbE::kIMCaloTracks) opt = "Conv";
  if (analysis == AliAnaPi0EbE::kSSCalo)       opt = "SS";
  if (analysis == AliAnaPi0EbE::kIMCalo && kCalorimeter=="EMCAL" && !kSimulation) ana->SetPairTimeCut(100);
  if (analysis == AliAnaPi0EbE::kIMCaloTracks) ana->SetInputAODGammaConvName("PhotonsCTS");

  // Common settings for all 3 type of analysis
  ana->SwitchOnSelectedClusterHistoFill();
  ana->SetCalorimeter(kCalorimeter);
  
  //Set Histograms name tag, bins and ranges
  ana->AddToHistogramsName(Form("Ana%s%sEbE_TM%d_",particle.Data(),opt.Data(),kTM));
  
  // Specific settings for different type of analysis
  
  ana->SwitchOffFillWeightHistograms();
  if (!kSimulation) ana->SwitchOnFillPileUpHistograms();
  
  if (kTM) {
    //printf("--->>>REMOVE MATCHED Pi0\n");
    ana->SwitchOnTrackMatchRejection() ;
    ana->SwitchOffTMHistoFill() ;
  } else {
    //printf("---->>>ACCEPT MATCHED Pi0\n");
    ana->SwitchOffTrackMatchRejection() ;
    ana->SwitchOnTMHistoFill() ;
  }
  
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  //ana->SwitchOnFillEMCALBCHistograms();

  if (!kInputDataType.Contains("delta")) {
    ana->SetOutputAODName(Form("%s%s%s",particle.Data(), opt.Data(), kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4Particle");
  } else
    ana->SetInputAODName(Form("%s%s%s",particle.Data(),opt.Data(),kName.Data()));
  
  // cluster splitting settings
  ana->SetMinEnergy(6);
  ana->SetMaxEnergy(200.);
  ana->SetNLMMinEnergy(0, 10);
  ana->SetNLMMinEnergy(1, 6);
  ana->SetNLMMinEnergy(2, 6);
  ana->SetMinDistanceToBadChannel(2, 4, 6); // only use the first one
  ana->SwitchOnSplitClusterDistToBad();
  ana->SetTimeCut(-1e10,1e10); // Open time cut

  // NLM cut, used in all, exclude clusters with more than 2 maxima
  ana->SetNLMCut(nlmMin, nlmMax) ;
    
  AliCaloPID* caloPID = ana->GetCaloPID();
  caloPID->SetSplitWidthSigma(3.); // cut at 3 sigma of the mean pi0 peak.

  if (!useSS) {
    printf("Do not apply SS cut on merged pi0 analysis \n");
    caloPID->SwitchOffSplitShowerShapeCut() ;
    ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_TM%d_",particle.Data(),opt.Data(),kTM));
    ana->SetOutputAODName(Form("%s%s%s_OpenSS",particle.Data(), opt.Data(), kName.Data()));
    caloPID->SetClusterSplittingM02Cut(0.1,10000); 
  } else {
    caloPID->SetClusterSplittingM02Cut(0.3,5); // Do the selection in the analysis class and not in the PID method to fill SS histograms
    caloPID->SwitchOnSplitShowerShapeCut() ;
    if (simpleSplitM02Cut) caloPID->SwitchOnSimpleSplitM02Cut();
  }

  if (simpleSplitMassCut) caloPID->SwitchOnSimpleSplitMassCut();

  if (useAsy) caloPID->SwitchOnSplitAsymmetryCut() ;
  else {
    caloPID->SwitchOffSplitAsymmetryCut() ;
    if (!useSS) {
      ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_OpenAsy_TM%d_",particle.Data(),opt.Data(),kTM));
      ana->SetOutputAODName(Form("%s%s%s_OpenSS_OpenAsy",particle.Data(), opt.Data(), kName.Data()));
    } else {
      ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenAsy_TM%d_",particle.Data(),opt.Data(),kTM));
      ana->SetOutputAODName(Form("%s%s%s_OpenAsy",particle.Data(), opt.Data(), kName.Data()));
    }
  }

  //For Pi0 only if  SwitchOnSimpleSplitMassCut()
  caloPID->SetPi0MassRange(massPi0Min, massPi0Max);
  caloPID->SetEtaMassRange(0.40, 0.60);
  caloPID->SetPhotonMassRange(0.00, 0.08);
  caloPID->SetClusterSplittingMinNCells(6);
  //caloPID->SetSplitEnergyFractionMinimum(0, 0.95);
  //caloPID->SetSplitEnergyFractionMinimum(1, 0.95);
  //caloPID->SetSplitEnergyFractionMinimum(2, 0.8);
  if (kCollisions=="PbPb" || kName.Contains("150")) {
    caloPID->SetClusterSplittingMinNCells(4);
    caloPID->SetPi0MassShiftHighECell(0.005);
  }
  
  if (kPrint) ana->Print("");
  return  ana;
}

///
/// Configure the task doing the pi0 flow.
/// Input AOD branch is from task AliAnaPi0EbE.
///
AliAnaPi0Flow* ConfigurePi0Flow()
{
  AliAnaPi0Flow *ana = new AliAnaPi0Flow();
  ana->SetDebug(kDebug);
      
  TString particle = "Pi0";
  TString opt = "SS";
  ana->SetInputAODName(Form("%s%s%s",particle.Data(),opt.Data(),kName.Data()));
  ana->IsPHOSCali(kPhosCali);
  ana->IsCentFlat(kCentFlat);

  if (kPrint) ana->Print("");
  return ana;
}

///
/// Set common histograms binning and ranges
///
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges)
{
  // Set common bins for all analysis and MC histograms filling
    
  histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
  
  if(kCalorimeter=="EMCAL")
  {
    if ( kYears == 2010 )
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
      histoRanges->SetHistoXRangeAndNBins(-230,90,120); // QA
      histoRanges->SetHistoYRangeAndNBins(370,450,40);  // QA
    }
    else if ( kYears < 2014 )
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
