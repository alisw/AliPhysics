
Bool_t  kPrint         = kFALSE;
Bool_t  kSimulation    = kFALSE;
Bool_t  kUseKinematics = kFALSE;
Bool_t  kOutputAOD     = kFALSE;
Int_t   kRunNumber     = 0;
Int_t   kYears         = 2011;
TString kCollisions    = "pp";
TString kTrig          = "EMC7" ;
TString kClusterArray  = "";
TString kData          = "ESD";
TString kInputDataType = "ESD";
TString kCalorimeter   = "EMCAL";

AliAnalysisTaskParticleCorrelation *AddTaskPartCorr(
                                                    const TString data          = "AOD",
                                                    const TString calorimeter   = "EMCAL", 
                                                    const Bool_t  printSettings = kFALSE,
                                                    const Bool_t  simulation    = kFALSE,
                                                    const Bool_t  outputAOD     = kFALSE, 
                                                    const TString outputfile    = "",
                                                    const Int_t   year          = 2010,
                                                    const Int_t   run           = 0,
                                                    const TString col           = "pp", 
                                                    const TString trigger       = "MB", 
                                                    const TString clustersArray = "V1" 
                                                    )
{
  // Creates a PartCorr task, configures it and adds it to the analysis manager.
  
  kPrint         = printSettings;
  kSimulation    = simulation;
  kRunNumber     = run;
  kYears         = year;
  kCollisions    = col;
  kTrig          = trigger;
  kClusterArray  = clustersArray;
  kData          = data;
  kCalorimeter   = calorimeter;
  kOutputAOD     = outputAOD;
  
  // Get the pointer to the existing analysis manager via the static access method.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }
  kInputDataType = "AOD";
  if(!kData.Contains("delta"))
    kInputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  Bool_t kUseKinematics = kFALSE; 
  if(kSimulation) { 
    kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE; 
    if (!kUseKinematics && data=="AOD" && kInputDataType != "ESD") kUseKinematics = kTRUE; //AOD primary should be available ... 
  } 
  
  cout<<"********* ACCESS KINE? "<<kUseKinematics<<endl;
  
  // #### Configure analysis ####
    
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  
  // General frame setting and configuration
  maker->SetReader   (ConfigureReader()   ); 
  maker->SetCaloUtils(ConfigureCaloUtils()); 
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important
  
  maker->AddAnalysis(ConfigureQAAnalysis()    , n++);  

  // Isolation settings
  Int_t partInCone = AliIsolationCut::kOnlyCharged;
  Int_t thresType  = AliIsolationCut::kPtThresIC;
  
  if(kClusterArray!=""){
    
    maker->AddAnalysis(ConfigurePhotonAnalysis(), n++); // Photon cluster selection
    maker->AddAnalysis(ConfigurePi0Analysis()   , n++); // Pi0 invariant mass analysis
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kIMCalo), n++); // Pi0 event by event selection, and photon tagging from decay
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Eta", AliAnaPi0EbE::kIMCalo), n++); // Eta event by event selection, and photon tagging from decay
    
    //maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kIMCaloTracks), n++); // Pi0 (calo+conversion) event by event selection, 
    // and photon tagging from decay, need to execute at the same time conversions analysis
    
    maker->AddAnalysis(ConfigureIsolationAnalysis("Photon",partInCone,thresType), n++); // Photon isolation
    maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0",   partInCone,thresType), n++); // Pi0 isolation
    
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kFALSE), n++); // Gamma hadron correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kTRUE) , n++); // Isolated gamma hadron correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0"   ,kFALSE), n++); // Pi0 hadron correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0"   ,kTRUE) , n++); // Isolated pi0 hadron correlation
  }
  else  
  {
    //Trigger on tracks, do only once, tracks do not depend on clusterizer
    maker->AddAnalysis(ConfigureChargedAnalysis(), n++);                                // track selection
    maker->AddAnalysis(ConfigureIsolationAnalysis("Hadron",partInCone,thresType), n++); // track isolation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Hadron",kFALSE), n++);       // track-track correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Hadron",kTRUE) , n++);       // Isolated track-track correlation
  }  
  
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(kData.Contains("delta")) maker->SwitchOffAODsMaker() ;
  else                        maker->SwitchOnAODsMaker()  ;
	
  if(kPrint) maker->Print("");
  
  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, kCalorimeter.Data());
  
  // Create task
  
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation (Form("PartCorr%s_Trig%s_Cl%s",kCalorimeter.Data(),kTrig.Data(),kClusterArray.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  task->SetBranches("ESD:AliESDRun.,AliESDHeader"); //just a trick to get Constantin's analysis to work
  task->SetAnalysisMaker(maker);
  mgr->AddTask(task);
  
  //Create containers
  
  if(outputfile.Length()==0)outputfile = AliAnalysisManager::GetCommonFileName(); 
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(Form("%s_Trig%s_Cl%s",kCalorimeter.Data(),kTrig.Data(),kClusterArray.Data()), TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s",outputfile.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("%sCuts_Trig%s_Cl%s",kCalorimeter.Data(),kTrig.Data(), kClusterArray.Data()), TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             Form("%s",outputfile.Data()));
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  if(!kData.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  
  
  if(kTrig=="EMC7"){
    printf("PartCorr trigger EMC7\n");
    task->SelectCollisionCandidates(AliVEvent::kEMC7);
  }
  else if (kTrig=="INT7"){
    printf("PartCorr trigger INT7\n");
    task->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  else if(kTrig=="EMC1"){
    printf("PartCorr trigger EMC1\n");
    task->SelectCollisionCandidates(AliVEvent::kEMC1);
  }
  else if(kTrig=="MB"){
    printf("PartCorr trigger MB\n");
    task->SelectCollisionCandidates(AliVEvent::kMB);
  }  
  else if(kTrig=="PHOS"){
    printf("PartCorr trigger PHOS\n");
    task->SelectCollisionCandidates(AliVEvent::kPHI7);
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
  
  reader->SetDebug(-1);//10 for lots of messages
  
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
  reader->SetEMCALEMin(0.5); 
  reader->SetEMCALEMax(1000); 
  reader->SetPHOSEMin(0.3);
  reader->SetPHOSEMax(1000);
  reader->SetCTSPtMin(0.1);
  reader->SetCTSPtMax(1000);
  
  reader->SwitchOffFiducialCut();
  
  // Tracks
  reader->SwitchOnCTS();
  gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/CreateTrackCutsPWG4.C"); 
  AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWG4(10011004); 
  reader->SetTrackCuts(esdTrackCuts);
  
  // Calorimeter
  
  reader->SetEMCALClusterListName(kClusterArray);
  if(kClusterArray == "") {
    printf("**************** Normal analysis **************** \n");
    reader->SwitchOnClusterRecalculation();
  }
  else {
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
  
  // Settings for LHC11a
  if(kRunNumber > 140000 && kRunNumber < = 146860){
    if(kClusterArray == "") reader->SwitchOnLEDEventsRemoval();
    reader->RejectFastClusterEvents();
    printf("Reader: Reject LED events and Fast cluster\n");
  }  
  
  //if(!kUseKinematics) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
  
  if     (kCollisions=="pp"  ) {
    if(kRunNumber < 140000) reader->SwitchOnEventSelection(); // remove pileup by default
    else                    reader->SwitchOffEventSelection(); 
    reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
    reader->SwitchOffPrimaryVertexSelection(); // and besides primary vertex
    reader->SetZvertexCut(50.);                // Open cut
  }
  else if(kCollisions=="PbPb") {
    reader->SwitchOffEventSelection();         // remove pileup by default
    reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
    reader->SwitchOffPrimaryVertexSelection(); // and besides primary vertex
    reader->SetZvertexCut(10.);                // Centrality defined in this range.
    
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt("10");  // 10 centrality bins
    reader->SetCentralityBin(-1,-1); // Accept all events, if not select range
    
    // Event plane (only used in AliAnaPi0 for the moment)
    reader->SetEventPlaneMethod("Q");
  }
  
  if(kPrint) reader->Print("");
  
  return reader;
  
}

//_______________________________________
AliCalorimeterUtils* ConfigureCaloUtils()
{
  
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(-1);
  
  if(kYears==2010) cu->SetEMCALGeometryName("EMCAL_FIRSTYEARV1");
  else             cu->SetEMCALGeometryName("EMCAL_COMPLETEV1");
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  if(kClusterArray == "") 
    cu->SwitchOffRecalculateClusterTrackMatching(); // Done in clusterization
  else            
    cu->SwitchOnRecalculateClusterTrackMatching();
  
  //EMCAL only settings
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  
  TGeoHMatrix* matrix[10];
  gROOT->LoadMacro("ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(
                          recou,
                          kSimulation, 
                          matrix,
                          "",//AODB path, default
                          kRunNumber, 
                          kPass
                          );   
  
  if(kCollisions=="pp"  ) { // Do only for pp for the moment
    cu->SwitchOnCorrectClusterLinearity();
    if(!kSimulation) recou->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrected);
    else             recou->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MC);
  }
  
  recou->SwitchOnRejectExoticCell();
  if(kClusterArray == "") recou->SwitchOnRejectExoticCluster();
  else                    recou->SwitchOffRejectExoticCluster();
  
  if(kPrint) cu->Print("");
  
  return cu;
  
}

//_____________________________________
AliAnaPhoton* ConfigurePhotonAnalysis()
{
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  
  // cluster selection cuts
  
  anaphoton->SwitchOffFiducialCut();
  
  anaphoton->SetCalorimeter(kCalorimeter);
  
  if(kCalorimeter == "PHOS"){
    anaphoton->SetNCellCut(2);// At least 2 cells
    anaphoton->SetMinPt(0.3);
    anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
  }
  else {//EMCAL
    anaphoton->SetNCellCut(1);// At least 2 cells
    anaphoton->SetMinPt(0.5); // avoid mip peak at E = 260 MeV
    anaphoton->SetMaxPt(1000); 
    if(!kUseKinematics) anaphoton->SetTimeCut(-1000,1000);// Time window of [400-900] ns if time recalibration is off, 
    // restrict to less than 100 ns when time calibration is on 
    anaphoton->SetMinDistanceToBadChannel(1, 2, 3); // For filtered AODs, new releases.
  }
  
  anaphoton->SwitchOnTrackMatchRejection() ;
  
  //PID cuts (shower shape)
  
  AliCaloPID* caloPID = anaphoton->GetCaloPID();
  anaphoton->SwitchOnCaloPID(); // if nothing else specified bayesian
  anaphoton->SwitchOnCaloPIDRecalculation(); // off, get bayesian weights, on use simple cut
  caloPID->SetLambda0CutMax(0.30);
  caloPID->SetLambda0CutMin(0.10);
  anaphoton->SwitchOffFillShowerShapeHistograms();  // Filled before photon shower shape selection
  
  // Input / output delta AOD settings
  
  if(!kData.Contains("delta")) {
    anaphoton->SetOutputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
    anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    //anaphoton->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
  }
  else anaphoton->SetInputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  
  //Set Histograms name tag, bins and ranges
  
  anaphoton->AddToHistogramsName("AnaPhoton_");
  SetHistoRangeAndNBins(anaphoton); // see method below
  
  // Number of particle type MC histograms
  anaphoton->FillNOriginHistograms(8);
  anaphoton->FillNPrimaryHistograms(4);
  
  if(kPrint) anaphoton->Print("");
  
  return anaphoton;
  
}

//_______________________________________________
AliAnaChargedParticles* ConfigureChargedAnalysis()
{
  
  AliAnaChargedParticles *anatrack = new AliAnaChargedParticles();
  anatrack->SetDebug(-1); //10 for lots of messages
  
  // selection cuts
  
  anatrack->SetDebug(-1);//10 for lots of messages
  anatrack->SetMinPt(5.);
  anatrack->SwitchOffFiducialCut();
  
  // Input / output delta AOD settings
  
  if(!kData.Contains("delta")) {
    anatrack->SetOutputAODName(Form("Hadron%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
    anatrack->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    //anaphoton->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
  }
  else anatrack->SetInputAODName(Form("Hadron%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  printf("Set Hadron%s_Trig%s_Cl%s\n",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data());
  //Set Histograms name tag, bins and ranges
  
  anatrack->AddToHistogramsName("AnaHadrons_");
  SetHistoRangeAndNBins(anatrack); // see method below

  if(kPrint) anatrack->Print("");
  
  return anatrack;
  
}


//_______________________________
AliAnaPi0* ConfigurePi0Analysis()
{
  
  AliAnaPi0 *anapi0 = new AliAnaPi0();
  
  anapi0->SetDebug(-1);//10 for lots of messages
  if(kPrint) anapi0->Print("");
  
  // Input delta AOD settings
  anapi0->SetInputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  
  // Calorimeter settings
  anapi0->SetCalorimeter(kCalorimeter);
  if(kCalorimeter=="PHOS") anapi0->SetNumberOfModules(3); //PHOS first year
  else {                   
    if(kYears==2010) anapi0->SetNumberOfModules(4); //EMCAL first year
    else             anapi0->SetNumberOfModules(10);
  }
  
  //settings for pp collision mixing
  anapi0->SwitchOnOwnMix(); //Off when mixing done with general mixing frame
  
  // Cuts 
  if(kCalorimeter=="EMCAL") anapi0->SetPairTimeCut(70);
  
  if     (kCollisions=="pp"  ) {
    anapi0->SetNCentrBin(1);
    anapi0->SetNZvertBin(10);
    anapi0->SetNRPBin(1);
    anapi0->SetNMaxEvMix(100);    
    anapi0->SwitchOffSMCombinations();
    anapi0->SwitchOffMultipleCutAnalysis();
  }
  else if(kCollisions=="PbPb") {
    anapi0->SetNCentrBin(10);
    anapi0->SetNZvertBin(10);
    anapi0->SetNRPBin(4);
    anapi0->SetNMaxEvMix(10);
    anapi0->SwitchOffSMCombinations();
  }
  
  //Set Histograms name tag, bins and ranges
  
  anapi0->AddToHistogramsName("AnaPi0_");
  SetHistoRangeAndNBins(anapi0); // see method below
  
  return anapi0;
  
}

//_____________________________________________________
AliAnaPi0EbE* ConfigurePi0EbEAnalysis(TString particle, 
                                      Int_t analysis)
{
  
  AliAnaPi0EbE *anapi0ebe = new AliAnaPi0EbE();
  anapi0ebe->SetDebug(-1);//10 for lots of messages
  
  anapi0ebe->SetAnalysisType(analysis);
  TString opt = "";
  if(analysis==AliAnaPi0EbE::kIMCaloTracks) opt = "Conv";
  if(analysis==AliAnaPi0EbE::kSSCalo)       opt = "SS";
  
  anapi0ebe->SwitchOffFillWeightHistograms();
  anapi0ebe->SetMinPt(0.5);
  if(kCalorimeter=="EMCAL")anapi0ebe->SetPairTimeCut(20); 
  anapi0ebe->SetCalorimeter(kCalorimeter);
  
  // Input / output delta AOD settings
  
  anapi0ebe->SetInputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  if(!kInputDataType.Contains("delta")) {
    anapi0ebe->SetOutputAODName(Form("%s%s%s_Trig%s_Cl%s",particle.Data(), opt.Data(), kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
    anapi0ebe->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else  anapi0ebe->SetInputAODName(Form("%s%s%s_Trig%s_Cl%s",particle.Data(),opt.Data(),kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  
  if(analysis == AliAnaPi0EbE::kIMCaloTracks) anapi0ebe->SetInputAODGammaConvName("PhotonsCTS");
  
  if(analysis!=AliAnaPi0EbE::kSSCalo){
    
    AliNeutralMesonSelection *nms = anapi0ebe->GetNeutralMesonSelection();
    nms->SetParticle(particle);
    nms->SwitchOnAngleSelection();
    nms->KeepNeutralMesonSelectionHistos(kTRUE);
    //nms->SetAngleMaxParam(2,0.2);
    nms->SetHistoERangeAndNBins(0, 20, 100) ;
    //nms->SetHistoIMRangeAndNBins(0, 1, 400);
  }
  
  //Set Histograms name tag, bins and ranges
  
  anapi0ebe->AddToHistogramsName(Form("Ana%s%sEbE_",particle.Data(),opt.Data()));
  SetHistoRangeAndNBins(anapi0ebe); // see method below
  
  if(kPrint) anapi0ebe->Print("");
  
  return  anapi0ebe;
  
}

//___________________________________________________________________
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle="Photon", 
                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                    Bool_t multi      = kFALSE)
{
  
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  
  anaisol->SetMinPt(5);
  
  // Input / output delta AOD settings
  
  anaisol->SetInputAODName(Form("%s%s_Trig%s_Cl%s",particle.Data(),kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  anaisol->SetAODObjArrayName(Form("IC%s",particle.Data())); 
  
  anaisol->SetCalorimeter(kCalorimeter);
  
  //Do settings for main isolation cut class
  AliIsolationCut * ic =  anaisol->GetIsolationCut();	
  ic->SetConeSize(0.4);
  ic->SetPtThreshold(0.5);
  ic->SetPtFraction(0.1);
  ic->SetSumPtThreshold(1.0) ;
  ic->SetParticleTypeInCone(partInCone);
  ic->SetICMethod(thresType);
  
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisol->SwitchOffReIsolation();
  //Multiple IC
  if(multi) anaisol->SwitchOnSeveralIsolation() ;
  else      anaisol->SwitchOffSeveralIsolation() ;
  
  //Set Histograms name tag, bins and ranges
  
  anaisol->AddToHistogramsName(Form("AnaIsol%s_",particle.Data()));
  SetHistoRangeAndNBins(anaisol); // see method below
  
  if(kPrint) ic     ->Print("");
  if(kPrint) anaisol->Print("");
  
  return anaisol;
  
}

//___________________________________________________________________________________
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle, 
                                                                    Int_t bIsolated)
{
  
  AliAnaParticleHadronCorrelation *anacorrhadron = new AliAnaParticleHadronCorrelation();
  anacorrhadron->SetDebug(-1);
  
  anacorrhadron->SetPtCutRange(5,200);
  
  // Input / output delta AOD settings
  
  anacorrhadron->SetInputAODName(Form("%s%s_Trig%s_Cl%s",particle.Data(),kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  anacorrhadron->AddToHistogramsName(Form("Ana%sHadronCorrIso%d_",particle.Data(),bIsolated));
  anacorrhadron->SetAODObjArrayName(Form("%sHadronCorrIso%d",particle.Data(),bIsolated)); 
  
  anacorrhadron->SelectIsolated(bIsolated); // do correlation with isolated photons
  
  anacorrhadron->SwitchOnDecayCorr();
  anacorrhadron->SetMultiBin(1);
  anacorrhadron->SwitchOffNeutralCorr();
  anacorrhadron->SwitchOffEventSelection();
  anacorrhadron->SetDeltaPhiCutRange(1.5,4.5);
  
  anacorrhadron->SwitchOnSeveralUECalculation();
  anacorrhadron->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  
  //if(kCalorimeter=="PHOS"){
  //Correlate with particles in EMCAL
  //anacorrhadron->SwitchOnCaloPID();
  //anacorrhadron->SwitchOnCaloPIDRecalculation(); 
  //}
  
  //Set Histograms name tag, bins and ranges
  
  anacorrhadron->AddToHistogramsName(Form("Ana%sHadronCorr_Iso%d_",particle.Data(),bIsolated));
  SetHistoRangeAndNBins(anacorrhadron); // see method below
  
  if(kPrint) anacorrhadron->Print("");  
  
  return anacorrhadron;
  
}

//________________________________________
AliAnaCalorimeterQA* ConfigureQAAnalysis()
{
  
  AliAnaCalorimeterQA *anaQA = new AliAnaCalorimeterQA();
  //anaQA->SetDebug(10); //10 for lots of messages
  anaQA->SetCalorimeter(kCalorimeter);
  
  anaQA->SwitchOffFiducialCut();
  anaQA->SwitchOnCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
  anaQA->SwitchOffFillAllTH3Histogram();
  anaQA->SwitchOffFillAllPositionHistogram();
  
  anaQA->SwitchOnStudyBadClusters() ;
  anaQA->SwitchOffStudyClustersAsymmetry();
  anaQA->SwitchOffStudyWeight();
  anaQA->SwitchOffFillAllTrackMatchingHistogram();
  
  if(kCalorimeter=="EMCAL"){
    if(kYears==2010)  
      anaQA->SetNumberOfModules(4); 
    else{           
      anaQA->SetNumberOfModules(10); 
    }
  } 
  anaQA->AddToHistogramsName("QA_"); //Begining of histograms name
  SetHistoRangeAndNBins(anaQA); // see method below
  
  if(kPrint) anaQA->Print("");	
  
  return anaQA;
  
}

//________________________________________________________
void SetHistoRangeAndNBins (AliAnaPartCorrBaseClass* ana)
{
  // Set common bins for all analysis and MC histograms filling
  
  if(kSimulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else            ana->SwitchOffDataMC() ;
  
  ana->SetHistoPtRangeAndNBins(0, 100, 250) ; // Energy and pt histograms
  
  if(kCalorimeter=="EMCAL"){
    if(kYears==2010){
      ana->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
      ana->SetHistoXRangeAndNBins(-230,90,120); // QA
      ana->SetHistoYRangeAndNBins(370,450,40);  // QA
    }
    else {           
      ana->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
      ana->SetHistoXRangeAndNBins(-600,90,200); // QA
      ana->SetHistoYRangeAndNBins(100,450,100); // QA
    }
    
    ana->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  }
  else{
    ana->SetHistoPhiRangeAndNBins(260*TMath::DegToRad(), 320*TMath::DegToRad(), 60) ;
    ana->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
    
  }
  
  ana->SetHistoShowerShapeRangeAndNBins(0, 3, 300);
  
  // Invariant mass analysis
  ana->SetHistoMassRangeAndNBins(0., 1., 200) ;
  ana->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  
  // check if time calibration is on
  ana->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  ana->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // QA, electron, charged
  ana->SetHistoPOverERangeAndNBins(0,10.,100);
  ana->SetHistodEdxRangeAndNBins(0.,200.,200);
  
  // QA
  ana->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  ana->SetHistodRRangeAndNBins(0.,TMath::Pi(),150);
  ana->SetHistoRatioRangeAndNBins(0.,2.,100);
  ana->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  ana->SetHistoNClusterCellRangeAndNBins(0,500,500);
  ana->SetHistoZRangeAndNBins(-400,400,200);
  ana->SetHistoRRangeAndNBins(400,450,25);
  ana->SetHistoV0SignalRangeAndNBins(0,5000,500);
  ana->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  ana->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
  
}

