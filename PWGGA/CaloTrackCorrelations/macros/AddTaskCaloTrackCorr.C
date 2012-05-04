
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
TString kData          = "ESD";
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

AliAnalysisTaskCaloTrackCorrelation *AddTaskCaloTrackCorr(const TString  data          = "AOD",
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
                                                          const Bool_t   recaltm       = kTRUE,
                                                          const Bool_t   tm            = kTRUE,
                                                          const Int_t    minCen        = -1,
                                                          const Int_t    maxCen        = -1,
                                                          const Bool_t   qaan          = kFALSE,
                                                          const Bool_t   hadronan      = kFALSE,
                                                          const Bool_t   calibE        = kTRUE,
                                                          const Bool_t   badmap        = kTRUE,
                                                          const Bool_t   calibT        = kTRUE,
                                                          const Bool_t   outputAOD     = kFALSE, 
                                                          const Bool_t   printSettings = kFALSE,
                                                          const Double_t scaleFactor   = -1
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
  
  kName = Form("%s_Trig%s_Cl%s_TM%d",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data(),kTM);

  if(kCollisions=="PbPb" && kMaxCen>=0) kName+=Form("Cen%d_%d",kMinCen,kMaxCen);
    
  printf("<<<< NAME: %s >>>>>\n",kName.Data());
  
  // #### Configure analysis ####
    
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  
  maker->SetScaleFactor(scaleFactor); // for MC, negative (not scaled) by default
  
  // General frame setting and configuration
  maker->SetReader   (ConfigureReader()   ); 
  maker->SetCaloUtils(ConfigureCaloUtils()); 
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important

  // Isolation settings
  Int_t partInCone = AliIsolationCut::kNeutralAndCharged; // kOnlyCharged;
  Int_t thresType  = AliIsolationCut::kPtThresIC;//  AliIsolationCut::kSumPtFracIC ; 
  
  // Photon analysis
  
  maker->AddAnalysis(ConfigurePhotonAnalysis(), n++); // Photon cluster selection
  maker->AddAnalysis(ConfigureIsolationAnalysis("Photon", partInCone,thresType), n++); // Photon isolation   
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kFALSE), n++); // Gamma hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kTRUE) , n++); // Isolated gamma hadron correlation
  //maker->AddAnalysis(ConfigureIsolationAnalysis("Photon", partInCone,thresType,kTRUE), n++); // Photon multi isolation, leave it the last   

  // Split cluster analysis
  if(kCalorimeter == "EMCAL")
  {
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kSSCalo), n++); // Pi0 event by event selection, cluster splitting   
    maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SS", partInCone,thresType), n++);       // Pi0 isolation, cluster splits   
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS" ,kFALSE), n++); // Pi0 hadron correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS" ,kTRUE) , n++); // Isolated pi0 hadron correlation
    //maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SS",  partInCone,thresType,kTRUE), n++); // Pi0 multi isolation, split cluster  
    //maker->AddAnalysis(ConfigureInClusterIMAnalysis(0.5,3), n++); 
  }
  
  // Invariant mass analysis
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kIMCalo), n++); // Pi0 event by event selection, invariant mass and photon tagging from decay    
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0SideBand", AliAnaPi0EbE::kIMCalo), n++); // Pi0 event by event selection, and photon tagging from decay    
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Eta", AliAnaPi0EbE::kIMCalo), n++); // Eta event by event selection, invariant mass and photon tagging from decay
  maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0", partInCone,thresType), n++);         // Pi0 isolation, invariant mass   
  maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SideBand", partInCone,thresType), n++); // Pi0 isolation, side band   
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0"   ,kFALSE), n++); // Pi0 hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0"   ,kTRUE) , n++); // Isolated pi0 hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SideBand" ,kFALSE), n++); // Pi0 hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SideBand" ,kTRUE) , n++); // Isolated pi0 hadron correlation
  //maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0",    partInCone,thresType,kTRUE), n++); // Pi0 multi isolation, invariant mass, leave it the last   

  if(kHadronAN)
  {
    maker->AddAnalysis(ConfigureChargedAnalysis(), n++);                                // track selection
    maker->AddAnalysis(ConfigureIsolationAnalysis("Hadron",AliIsolationCut::kOnlyCharged,thresType), n++); // track isolation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Hadron",kFALSE), n++);       // track-track correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Hadron",kTRUE) , n++);       // Isolated track-track correlation
    //maker->AddAnalysis(ConfigureIsolationAnalysis("Hadron",partInCone,thresType,kTRUE), n++);// Hadron multi isolation  
  }
  
  // Analysis with ghost triggers, only for Min Bias like events
  if( kTrig.Contains("INT") || kTrig.Contains("Central") || kTrig.Contains("MB")  )
  {
    maker->AddAnalysis(ConfigureRandomTriggerAnalysis(), n++); 
    maker->AddAnalysis(ConfigureIsolationAnalysis(Form("RandomTrigger%s",kCalorimeter.Data()), partInCone,thresType), n++); // Ghost trigger isolation   
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis(Form("RandomTrigger%s",kCalorimeter.Data()),kFALSE), n++); // Ghost trigger hadron correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis(Form("RandomTrigger%s",kCalorimeter.Data()),kTRUE) , n++); // Isolated ghost hadron correlation
    //maker->AddAnalysis(ConfigureIsolationAnalysis(Form("RandomTrigger%s",kCalorimeter.Data()), partInCone,thresType,kTRUE), n++); // Ghost multi isolation   
    
    if(kHadronAN)
    {
      maker->AddAnalysis(ConfigureRandomTriggerAnalysis("CTS"), n++);                                // track selection
      maker->AddAnalysis(ConfigureIsolationAnalysis("RandomTriggerCTS",AliIsolationCut::kOnlyCharged,thresType), n++); // track isolation
      maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("RandomTriggerCTS",kFALSE), n++);       // track-track correlation
      maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("RandomTriggerCTS",kTRUE) , n++);       // Isolated track-track correlation
      //maker->AddAnalysis(ConfigureIsolationAnalysis("RandomTriggerCTS",AliIsolationCut::kOnlyCharged,thresType,kTRUE), n++); // Ghost multi isolation   

    }
  }
  
  if(kQA)  maker->AddAnalysis(ConfigureQAAnalysis(),n++);
  
  maker->SetAnaDebug(kDebug)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(kData.Contains("delta")) maker->SwitchOffAODsMaker() ;
  else                        maker->SwitchOnAODsMaker()  ;
  
  if(kPrint) maker->Print("");
  
  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, kCalorimeter.Data());
  // CAREFUL
  //kName = Form("%s_Trig%s_Cl%s_TM%d",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data(),kFALSE);
  kName = Form("%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data());
  if(kCollisions=="PbPb" && kMaxCen>=0) kName+=Form("Cen%d_%d",kMinCen,kMaxCen);

  // Create task
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("CaloTrackCorr%s",kName.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(kDebug);
  task->SetBranches("ESD:AliESDRun.,AliESDHeader"); 
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
  
  
  if(kTrig=="EMC7")
  {
    printf("CaloTrackCorr trigger EMC7\n");
    task->SelectCollisionCandidates(AliVEvent::kEMC7);
  }
  else if (kTrig=="INT7")
  {
    printf("CaloTrackCorr trigger INT7\n");
    task->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  else if(kTrig=="EMC1")
  {
    printf("CaloTrackCorr trigger EMC1\n");
    task->SelectCollisionCandidates(AliVEvent::kEMC1);
  }
  else if(kTrig=="MB")
  {
    printf("CaloTrackCorr trigger MB\n");
    task->SelectCollisionCandidates(AliVEvent::kMB);
  }  
  else if(kTrig=="PHOS")
  {
    printf("CaloTrackCorr trigger PHOS\n");
    task->SelectCollisionCandidates(AliVEvent::kPHI7);
  }  
  else if(kTrig=="PHOSPb")
  {
    printf("CaloTrackCorr trigger PHOSPb\n");
    task->SelectCollisionCandidates(AliVEvent::kPHOSPb);
  }
  else if(kTrig=="AnyINT")
  {
    printf("CaloTrackCorr trigger AnyINT\n");
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  }  
  else if(kTrig=="INT")
  {
    printf("CaloTrackCorr trigger AnyINT\n");
    task->SelectCollisionCandidates(AliVEvent::kAny);
  }
  else if(kTrig=="EMCEGA")
  {
    printf("CaloTrackCorr trigger EMC Gamma\n");
    task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  } 
  else if(kTrig=="EMCEJE")
  {
    printf("CaloTrackCorr trigger EMC Jet\n");
    task->SelectCollisionCandidates(AliVEvent::kEMCEJE);
  }
  else if(kTrig=="Central")
  {
    printf("CaloTrackCorr trigger Central\n");
    task->SelectCollisionCandidates(AliVEvent::kCentral);
  } 
  else if(kTrig=="SemiCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral\n");
    task->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  }
  
  return task;
}

//____________________________________
AliCaloTrackReader * ConfigureReader()
{
  
  AliCaloTrackReader * reader = 0;
  if     (kData.Contains("AOD"))   reader = new AliCaloTrackAODReader();
  else if(kData=="ESD")            reader = new AliCaloTrackESDReader();
  else if(kData=="MC" && 
          kInputDataType == "ESD") reader = new AliCaloTrackMCReader();
  
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

  if(!kSimulation && kCalibT) reader->SetEMCALTimeCut(-30,30); 
  else                        reader->SetEMCALTimeCut(-1000,1000); // Open time cut

  reader->SwitchOnFiducialCut();
  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;

  // Tracks
  reader->SwitchOnCTS();
  if(kInputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");
    AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
    reader->SetTrackCuts(esdTrackCuts);
    reader->SwitchOnConstrainTrackToVertex();
  }
  else if(kInputDataType=="AOD")
  {
    reader->SwitchOnAODHybridTrackSelection();
    reader->SetTrackFilterMask(128); // Filter bit, not mask, use if off hybrid
  }
  
  // Calorimeter
  
  reader->SetEMCALClusterListName(kClusterArray);
  if(kClusterArray == "") 
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
  
  //if(kCalorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();  
    reader->SwitchOnEMCAL();
  //}
  //if(kCalorimeter == "PHOS") { 
    reader->SwitchOnPHOSCells();  
    reader->SwitchOnPHOS();
  //}
  
  // for case data="deltaAOD", no need to fill the EMCAL/PHOS cluster lists
  if(kData.Contains("delta")){
    reader->SwitchOffEMCAL();
    reader->SwitchOffPHOS();
    reader->SwitchOffEMCALCells(); 
    reader->SwitchOffPHOSCells(); 
  }
  
  //-----------------
  // Event selection
  //-----------------
  
  //if(!kUseKinematics) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
  
  // For mixing with AliAnaParticleHadronCorrelation switch it off
  reader->SwitchOnEventTriggerAtSE();
  reader->SetEventTriggerMaks(AliVEvent::kEMC7); // Only for mixing and SwitchOffEventTriggerAtSE();
  
  reader->SetZvertexCut(10.);                // Open cut
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex

  if(kEventSelection)
  {
    reader->SwitchOnEventSelection();         // remove pileup by default
    reader->SwitchOnV0ANDSelection() ;        // and besides v0 AND
  }
  else 
  {
    reader->SwitchOffEventSelection();         // remove pileup by default
    reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
  }
    
  if(kCollisions=="PbPb") 
  {
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(10);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
    reader->SetCentralityBin(kMinCen,kMaxCen); // Accept all events, if not select range
    
    // Event plane (only used in AliAnaPi0 for the moment)
    reader->SetEventPlaneMethod("Q");
  }
  
  reader->SetImportGeometryFromFile(kTRUE);
  
  if(kPrint) reader->Print("");
  
  return reader;
  
}

//_______________________________________
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
  

  cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
  
  
  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          kSimulation,                             
                          kExotic,
                          kNonLinearity,
                          kCalibE, 
                          kBadMap,
                          kCalibT);   
  
  if( kNonLinearity ) 
  { 
    printf("ConfigureCaloUtils() - Apply non linearity to EMCAL\n");
    cu->SwitchOnCorrectClusterLinearity();
  }
    
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
  
  ana->SwitchOffFiducialCut();

  ana->SetCalorimeter(kCalorimeter);
  
  if(kCalorimeter == "PHOS")
  {
    ana->SetNCellCut(2);// At least 3 cells
    ana->SetMinPt(0.3);
    ana->SetMinDistanceToBadChannel(2, 4, 5);
    ana->SetTimeCut(-2000,2000); // open cut
  }
  else 
  {//EMCAL
    ana->SetNCellCut(1);// At least 2 cells
    ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
    ana->SetMaxEnergy(1000); 
    ana->SetTimeCut(-1000,1000); // open cut, usual time window of [425-825] ns if time recalibration is off 
    // restrict to less than 100 ns when time calibration is on 
    ana->SetMinDistanceToBadChannel(2, 4, 6); 
  }
  
  if(kTM)
  {
    ana->SwitchOnTrackMatchRejection() ;
    ana->SwitchOffTMHistoFill() ;
  }
  else
  {
    ana->SwitchOffTrackMatchRejection() ;
    ana->SwitchOnTMHistoFill() ;
  }
  
  //PID cuts (shower shape)
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
  //Not used in bayesian
  
  //EMCAL
  caloPID->SetEMCALLambda0CutMax(0.27);
  caloPID->SetEMCALLambda0CutMin(0.10);
  
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
    
  //PHOS
  caloPID->SetPHOSDispersionCut(2.5);
  caloPID->SetPHOSRCut(2.);
  if(kData=="AOD") caloPID->SetPHOSRCut(2000.); // Open cut since dX, dZ not stored
      
  ana->SwitchOffFillShowerShapeHistograms();  // Filled before photon shower shape selection
  
  // Input / output delta AOD settings
  
  if(!kData.Contains("delta")) 
  {
    ana->SetOutputAODName(Form("Photon%s",kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    //ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
  }
  else ana->SetInputAODName(Form("Photon%s",kName.Data()));
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaPhoton_TM%d_",kTM));
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms(20);
  ana->FillNPrimaryHistograms(20);
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana;
  
}

//________________________________________________________________________________
AliAnaElectron* ConfigureElectronAnalysis()
{
  
  AliAnaElectron *ana = new AliAnaElectron();
  ana->SetDebug(kDebug); //10 for lots of messages
    
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
  if(kSimulation)
  { // LHC11a
    ana->SetdEdxCut(76, 85);
  }
  else // LHC11c
  {
    ana->SetdEdxCut(56, 64);
  }
  
  ana->SetCalorimeter(kCalorimeter);
  
  ana->SwitchOffCaloPID();
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

//__________________________________________________________________________________________
AliAnaRandomTrigger* ConfigureRandomTriggerAnalysis(TString detector = "")
{
  
  AliAnaRandomTrigger *ana = new AliAnaRandomTrigger();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  if(detector=="") detector = kCalorimeter;
  ana->SetDetector(detector);

  // selection cuts
  ana->SetMinPt(4.); 
  ana->SetMaxPt(51.);   
  
  if     (detector=="EMCAL")
  {
    ana->SetEtaCut(-0.71,0.71);
    ana->SetPhiCut(100*TMath::DegToRad(), 160*TMath::DegToRad());
  }
  else if(detector=="PHOS")
  {
    ana->SetEtaCut(-0.13,0.13);
    ana->SetPhiCut(260*TMath::DegToRad(), 320*TMath::DegToRad());
  }
  else if(detector=="CTS")
  {
    ana->SetEtaCut(-0.9,0.9);
    ana->SetPhiCut(0, TMath::TwoPi());
  }
  
  // AOD branch
  if(!kData.Contains("delta")) 
  {
    ana->SetOutputAODName(Form("RandomTrigger%s%s",detector.Data(),kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else 
    ana->SetInputAODName(Form("RandomTrigger%s%s",detector.Data(),kName.Data()));
  
  printf("Set RandomTrigger%s%s\n",detector.Data(),kName.Data());
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaRandomTrigger%s_",detector.Data()));
  
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  if(detector=="CTS")
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  if(kPrint) ana->Print("");
  
  return ana;
  
}

//__________________________________________________________________________________________
AliAnaInsideClusterInvariantMass* ConfigureInClusterIMAnalysis(Float_t l0min, Float_t l0max)
{
  AliAnaInsideClusterInvariantMass *ana = new AliAnaInsideClusterInvariantMass();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // selection cuts
  
  ana->SetMinEnergy(5); 
  ana->SetMaxEnergy(200.);   
  ana->SetMinNCells(3);
  ana->SetM02Cut(l0min,l0max);
  ana->SetCalorimeter(kCalorimeter);
  
  //ana->AddToHistogramsName(Form("AnaInClusterIM_%1.2f_%1.2f_",l0min,l0max));
  ana->AddToHistogramsName("AnaInClusterIM_");

  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
    
  AliCaloPID* caloPID = ana->GetCaloPID();
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
  caloPID->SetClusterSplittingM02Cut(0,100); // Do the selection in the analysis class and not in the PID method to fill SS histograms

  caloPID->SetPi0MassRange(0.10, 0.18);
  caloPID->SetEtaMassRange(0.40, 0.60);
  caloPID->SetPhotonMassRange(0.00, 0.08);
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana;
  
}

//_______________________________________________
AliAnaChargedParticles* ConfigureChargedAnalysis()
{
  
  AliAnaChargedParticles *ana = new AliAnaChargedParticles();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // selection cuts
  
  ana->SetMinPt(0.5);
  ana->SwitchOnFiducialCut();
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ; //more restrictive cut in reader and after in isolation

  // Input / output delta AOD settings
  
  if(!kData.Contains("delta")) 
  {
    ana->SetOutputAODName(Form("Hadron%s",kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    //ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
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


//_______________________________
AliAnaPi0* ConfigurePi0Analysis()
{
  
  AliAnaPi0 *ana = new AliAnaPi0();
  
  ana->SetDebug(kDebug);//10 for lots of messages
  
  // Input delta AOD settings
  ana->SetInputAODName(Form("Photon%s",kName.Data()));
  
  // Calorimeter settings
  ana->SetCalorimeter(kCalorimeter);
  if(kCalorimeter=="PHOS") ana->SetNumberOfModules(3); //PHOS first year
  else 
  {                   
    if     (kYears == 2010) ana->SetNumberOfModules( 4); // EMCAL first year
    else if(kYears == 2011) ana->SetNumberOfModules(10); // Second year
    else                    ana->SetNumberOfModules(12);
  }
  
  //settings for pp collision mixing
  ana->SwitchOnOwnMix(); //Off when mixing done with general mixing frame
  
  // Cuts 
  if(kCalorimeter=="EMCAL") ana->SetPairTimeCut(70);
  
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

  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana;
  
}

//_____________________________________________________
AliAnaPi0EbE* ConfigurePi0EbEAnalysis(TString particle, 
                                      Int_t analysis)
{
  
  AliAnaPi0EbE *ana = new AliAnaPi0EbE();
  ana->SetDebug(kDebug);//10 for lots of messages
  
  ana->SetAnalysisType(analysis);
  TString opt = "";
  if(analysis==AliAnaPi0EbE::kIMCaloTracks) opt = "Conv";
  if(analysis==AliAnaPi0EbE::kSSCalo)       opt = "SS";
  
  ana->SetMinPt(0.5);
  
  if(kCalorimeter=="EMCAL")ana->SetPairTimeCut(15); // More strict than in pi0 inv mass analysis
  
  ana->SetCalorimeter(kCalorimeter);
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("Photon%s",kName.Data()));
  if(!kInputDataType.Contains("delta")) 
  {
    ana->SetOutputAODName(Form("%s%s%s",particle.Data(), opt.Data(), kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else  
    ana->SetInputAODName(Form("%s%s%s",particle.Data(),opt.Data(),kName.Data()));
  
  if(analysis == AliAnaPi0EbE::kIMCaloTracks) ana->SetInputAODGammaConvName("PhotonsCTS");
  
  if(analysis!=AliAnaPi0EbE::kSSCalo)
  {
    AliNeutralMesonSelection *nms = ana->GetNeutralMesonSelection();
    nms->SetParticle(particle);
    
    // Tighten a bit mass cut with respect to default window
    if(particle=="Pi0") nms->SetInvMassCutRange(0.120,0.150);
    if(particle=="Eta") nms->SetInvMassCutRange(0.520,0.580);    
    
    nms->SwitchOnAngleSelection();
    nms->KeepNeutralMesonSelectionHistos(kTRUE);
    //nms->SetAngleMaxParam(2,0.2);
    nms->SetHistoERangeAndNBins(0, 20, 80) ;
    //nms->SetHistoIMRangeAndNBins(0, 1, 400);
  }
  else
  { // cluster splitting settings
    ana->SetTimeCut(-1000,1000); // Open time cut
    AliCaloPID* caloPID = ana->GetCaloPID();
    caloPID->SetPi0MassRange(0.10, 0.18);
    caloPID->SetEtaMassRange(0.40, 0.60);
    caloPID->SetPhotonMassRange(0.00, 0.08);
    caloPID->SetClusterSplittingM02Cut(0.5,100); // Do the selection in the analysis class and not in the PID method to fill SS histograms
  }
  
  ana->SwitchOffSelectedClusterHistoFill(); 
  ana->SwitchOffFillWeightHistograms();
  
  if(!kTM) ana->SwitchOnTMHistoFill();
  else     ana->SwitchOffTMHistoFill();
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("Ana%s%sEbE_TM%d_",particle.Data(),opt.Data(),kTM));
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return  ana;
  
}

//____________________________________________________________________________________________________
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle="Photon", 
                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                    Bool_t multi      = kFALSE)
{
  
  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
  //ana->SetDebug(kDebug);
  ana->SetDebug(kDebug);
  
  ana->SwitchOnFiducialCut();
  //Avoid borders of EMCal
  if(kCalorimeter=="EMCAL")
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;

  // Same Eta as EMCal, cut in phi if EMCAL was triggering
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    if(kTrig.Contains("EMC"))
      ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
    else
      ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;    
  }
  
  ana->SetMinPt(4);
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s",particle.Data(),kName.Data()));
  ana->SetAODObjArrayName(Form("IC%sTM%d",particle.Data(),kTM)); 
  
  ana->SetCalorimeter(kCalorimeter);
  
  if(!kTM)  ana->SwitchOnTMHistoFill();
  else      ana->SwitchOffTMHistoFill();
  
  ana->SwitchOffSSHistoFill();
  
  //Do settings for main isolation cut class
  AliIsolationCut * ic =  ana->GetIsolationCut();	
  ic->SetDebug(kDebug);
  
  
  if(kCollisions=="pp") 
  {
    ic->SetPtThreshold(0.5);
    ic->SetConeSize(0.5);
  }
  if(kCollisions=="PbPb")
  {
    ic->SetPtThreshold(2.);
  //ic->SetPtThreshold(1.);
    ic->SetConeSize(0.3);
  }
  
  ic->SetPtFraction(0.1);
  ic->SetSumPtThreshold(1.0) ;
  ic->SetParticleTypeInCone(partInCone);
  ic->SetICMethod(thresType);
  
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
     
    ana->SetNCones(4);
    ana->SetNPtThresFrac(4);
    ana->SetConeSizes(0,0.3);       ana->SetConeSizes(1,0.4);
    ana->SetConeSizes(2,0.5);       ana->SetConeSizes(3,0.6);
    ana->SetPtThresholds(0, 0.5);   ana->SetPtThresholds(1, 1);     ana->SetPtThresholds(2, 2);
    ana->SetPtFractions (0, 0.05) ; ana->SetPtFractions (1, 0.1);   ana->SetPtFractions (2, 0.2) ;  ana->SetPtFractions (3, 0.3) ;
    ana->SetSumPtThresholds(0, 1) ; ana->SetSumPtThresholds(1, 3) ; ana->SetSumPtThresholds(2, 5);  ana->SetSumPtThresholds(3, 7)  ;
    
    ana->SwitchOffTMHistoFill();
    ana->SwitchOffSSHistoFill();
  }
  else      
    ana->SwitchOffSeveralIsolation() ;
  
  AliCaloPID* caloPID = ana->GetCaloPID();
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
  
  //Set Histograms name tag, bins and ranges
  
  if(!multi)ana->AddToHistogramsName(Form("AnaIsol%s_TM%d_",particle.Data(),kTM));
  else      ana->AddToHistogramsName(Form("AnaMultiIsol%s_TM%d_",particle.Data(),kTM));

  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  ana->SetHistoPtInConeRangeAndNBins(0, 50 , 250);
  ana->SetHistoPtSumRangeAndNBins   (0, 100, 250);
  
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
                                                                    Int_t bIsolated)
{
  
  AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
  ana->SetDebug(kDebug);
  
  ana->SetMinimumTriggerPt(4);
  ana->SetAssociatedPtRange(0.2,200); 
  ana->SetDeltaPhiCutRange( TMath::Pi()/2,3*TMath::Pi()/2 ); //[90 deg, 270 deg]
  
  ana->SelectIsolated(bIsolated); // do correlation with isolated photons
  
  ana->SwitchOnAbsoluteLeading();  // Select trigger leading particle of all the selected tracks
  ana->SwitchOffNearSideLeading(); // Select trigger leading particle of all the particles at +-90 degrees, default
  
  // Mixing with own pool
  ana->SwitchOffOwnMix();
  ana->SetNCentrBin(9);
  ana->SetNZvertBin(10);
  if(kCollisions=="pp")
  {
    ana->SetNMaxEvMix(10);    
    ana->SwitchOnTrackMultBins();
  }
  else 
  {
    ana->SetNMaxEvMix(2);    
    ana->SwitchOffTrackMultBins();
  }  
  
  ana->SwitchOnFiducialCut();
  
  //Avoid borders of EMCal, same as for isolation
  if(kCalorimeter=="EMCAL")
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;
  
  // Same Eta as EMCal, cut in phi if EMCAL was triggering
  if(particle=="Hadron" || particle.Contains("CTS"))
  {
    if(kTrig.Contains("EMC"))
      ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
    else
      ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;    
  }
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s",particle.Data(),kName.Data()));
  ana->SetAODObjArrayName(Form("%sHadronCorrIso%dTM%d",particle.Data(),bIsolated,kTM)); 
  
  // Fill extra plots on tagged decay photons
  // If trigger is pi0/eta found with invariant mass, get the decays
  // If trigger is photon, check if it was tagged as decay previously
  if(particle!="Hadron" )
  {
    if(particle=="Pi0" || particle =="Eta") 
    {
      if(!particle.Contains("SS")) ana->SwitchOnPi0TriggerDecayCorr();
      else                         ana->SwitchOffPi0TriggerDecayCorr();
      ana->SwitchOffDecayTriggerDecayCorr();
    }
    else
    {
      ana->SwitchOffPi0TriggerDecayCorr();
      ana->SwitchOnDecayTriggerDecayCorr(); // Make sure pi0 decay tagging runs before this task
    }
  }
  else
  {
    ana->SwitchOffPi0TriggerDecayCorr();
    ana->SwitchOffDecayTriggerDecayCorr(); 
  }
  
  // if triggering on PHOS and EMCAL is on
  //if(kCalorimeter=="PHOS") ana->SwitchOnNeutralCorr();
  ana->SwitchOffNeutralCorr(); // Do only correlation with TPC
  
  ana->SwitchOffHMPIDCorrelation();
  
  ana->SwitchOffFillBradHistograms();
  
  // Underlying event
  ana->SwitchOffEventSelection();
  ana->SwitchOnSeveralUECalculation();
  ana->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  ana->SetMultiBin(1);
  
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
  
  return ana;
  
}

//________________________________________
AliAnaCalorimeterQA* ConfigureQAAnalysis()
{
  
  AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
  ana->SetDebug(kDebug); //10 for lots of messages
  ana->SetCalorimeter(kCalorimeter);
  
  ana->SetTimeCut(-1000,1000); // Open time cut
  
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
  
  
  ana->SwitchOffFiducialCut();
  ana->SwitchOffFillAllTH3Histogram();
  ana->SwitchOffFillAllPositionHistogram();
  ana->SwitchOffFillAllPositionHistogram2();
  if(!kExotic)ana->SwitchOnStudyBadClusters();
  else        ana->SwitchOffStudyBadClusters();
  ana->SwitchOffStudyClustersAsymmetry();
  ana->SwitchOffStudyWeight();
  ana->SwitchOnFillAllTrackMatchingHistogram();
  
  if(kCalorimeter=="EMCAL")
  {
    if     (kYears==2010)  ana->SetNumberOfModules(4); 
    else if(kYears==2011)  ana->SetNumberOfModules(10);
    else                   ana->SetNumberOfModules(12); 
  }
  else 
  {//PHOS
    ana->SetNumberOfModules(3); 
  }
  
  ana->AddToHistogramsName("QA_"); //Begining of histograms name
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  ConfigureMC(ana);
  
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
  histoRanges->SetHistoZRangeAndNBins(-400,400,200);
  histoRanges->SetHistoRRangeAndNBins(400,450,25);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
  
}

