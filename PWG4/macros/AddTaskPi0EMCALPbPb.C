
AliAnalysisTaskParticleCorrelation *AddTaskPi0EMCALPbPb(
                                                    const TString data          = "AOD",
                                                    const TString calorimeter   = "EMCAL", 
                                                    const Bool_t  printSettings = kFALSE,
                                                    const Bool_t  simulation    = kFALSE,
                                                    const Bool_t  outputAOD     = kFALSE, 
                                                    const TString outputfile    = "",
                                                    const Int_t   year          = 2011,
                                                    const TString col           = "PbPb",
                                                    const TString trig          = "",
                                                    const TString clustersArray = "" 
                                                    )
{
  // Creates a PartCorr task, configures it and adds it to the analysis manager.
    
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
  Bool_t kInputDataType = "AOD";
  if(!data.Contains("delta"))
    kInputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  Bool_t kUseKinematics = kFALSE; 
  if(simulation) { 
    kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE; 
    if (!kUseKinematics && data=="AOD" && kInputDataType != "ESD") kUseKinematics = kTRUE; //AOD primary should be available ... 
  } 
  
  //cout<<"********* ACCESS KINE? "<<kUseKinematics<<endl;
  
  // #### Configure analysis ####
    
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  
  // General frame setting and configuration
  maker->SetReader   (ConfigureReader(data,kInputDataType,calorimeter,clustersArray,col,outputAOD,kUseKinematics,printSettings)   ); 
  maker->SetCaloUtils(ConfigureCaloUtils(year,col,clustersArray, simulation, printSettings)); 
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important
  maker->AddAnalysis(ConfigurePhotonAnalysis(data,calorimeter,year,col,clustersArray,simulation,"",printSettings), n++); // Photon cluster selection
  maker->AddAnalysis(ConfigurePi0Analysis   (data,calorimeter,year,col,clustersArray,simulation,"",printSettings), n++); // Pi0 invariant mass analysis
  
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(data.Contains("delta")) maker->SwitchOffAODsMaker() ;
  else                       maker->SwitchOnAODsMaker()  ;
	
  if(printSettings) maker->Print("");
  
  //printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, calorimeter.Data());
  
  // Create task
  
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation (Form("PartCorr%s_Trig%s_Cl%s",calorimeter.Data(),trig.Data(),clustersArray.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  task->SetBranches("ESD:AliESDRun.,AliESDHeader"); //just a trick to get Constantin's analysis to work
  task->SetAnalysisMaker(maker);
  mgr->AddTask(task);
  
  //Create containers
  
  if(outputfile.Length()==0)outputfile = AliAnalysisManager::GetCommonFileName(); 
  TString name(Form("%s_Trig%s_Cl%s",calorimeter.Data(),trig.Data(),clustersArray.Data()));
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(name.Data(), TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s:%s",outputfile.Data(),name.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("%sCuts",name.Data()), TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             Form("%s:%sCuts",outputfile.Data(),name.Data()));

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  if(!data.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
    
  return task;
}

//____________________________________
AliCaloTrackReader * ConfigureReader(TString kData,TString kInputDataType, TString kCalorimeter, 
                                     TString kClusterArray, TString kCollisions,
                                     Bool_t kOutputAOD, Bool_t kUseKinematics, Bool_t kPrint)
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
  reader->SetEMCALPtMin(0.5); 
  reader->SetEMCALPtMax(1000); 
  reader->SetPHOSPtMin(0.3);
  reader->SetPHOSPtMax(1000);
  reader->SetCTSPtMin(0.1);
  reader->SetCTSPtMax(1000);
  
  reader->SwitchOffFiducialCut();
  
  // Tracks
  reader->SwitchOffCTS();
  
  //gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/CreateTrackCutsPWG4.C"); 
  //AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWG4(10011004); 
  //reader->SetTrackCuts(esdTrackCuts);
  
  // Calorimeter
  
  reader->SetEMCALClusterListName(kClusterArray);
//  if(kClusterArray == "") {
//    //printf("**************** Normal analysis **************** \n");
//    reader->SwitchOnClusterRecalculation(); // Bad map removal
//  }
//  else {
//    printf("**************** Input for analysis is Clusterizer %s **************** \n", kClusterArray.Data());
//    reader->SwitchOffClusterRecalculation();
//  }  
  
  if(kCalorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();  
    reader->SwitchOnEMCAL();
  }
  if(kCalorimeter == "PHOS") { 
    reader->SwitchOnPHOSCells();  
    reader->SwitchOnPHOS();
  }
  
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
  
  if     (kCollisions=="pp"  ) {
    reader->SwitchOffEventSelection(); 
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
    reader->SetCentralityOpt(10);  // 10 centrality bins
    reader->SetCentralityBin(-1,-1); // Accept all events, if not select range
    
    // Event plane (only used in AliAnaPi0 for the moment)
    reader->SetEventPlaneMethod("Q");
  }
  
  if(kPrint) reader->Print("");
  
  return reader;
  
}

//_______________________________________
AliCalorimeterUtils* ConfigureCaloUtils(Int_t kYears, TString kCollisions, TString kClusterArray, 
                                        Bool_t kSimulation, Bool_t kPrint)
{
  
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(-1);
  
  if(kYears==2010) cu->SetEMCALGeometryName("EMCAL_FIRSTYEARV1");
  else             cu->SetEMCALGeometryName("EMCAL_COMPLETEV1");
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  //if(kClusterArray == "") 
  //  cu->SwitchOffRecalculateClusterTrackMatching(); // Done in clusterization
  //else            
  //  cu->SwitchOnRecalculateClusterTrackMatching();
  
  //EMCAL only settings
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  
  
  if(kCollisions=="pp"  ) { // Do only for pp for the moment
    cu->SwitchOnCorrectClusterLinearity();
    if(!kSimulation) recou->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrected);
    else             recou->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MC);
  }
  
  //recou->SwitchOnRejectExoticCell();
  //if(kClusterArray == "") recou->SwitchOnRejectExoticCluster();
  //else                    recou->SwitchOffRejectExoticCluster();
  
  if(kPrint) cu->Print("");
  
  return cu;
  
}

//_____________________________________
AliAnaPhoton* ConfigurePhotonAnalysis(TString kData, TString kCalorimeter, Int_t kYears, TString kCollisions, TString kClusterArray, 
                                      Bool_t kSimulation, TString kTrig = "", Bool_t kPrint)
{
  
  printf("Photon kData %s, kCalorimeter %s, kYears %d, kCollisions %s, kClusterArray %s, kSimulation %d, kTrig %s, kPrint %d \n",
         kData.Data(), kCalorimeter.Data(), kYears, kCollisions.Data(), kClusterArray.Data(), kSimulation, kTrig.Data(), kPrint);
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
    anaphoton->SetTimeCut(-1000,1000);// Time window of [400-900] ns if time recalibration is off, 
    // restrict to less than 100 ns when time calibration is on 
    anaphoton->SetMinDistanceToBadChannel(1, 2, 3); // For filtered AODs, new releases.
  }
  
  anaphoton->SwitchOnTrackMatchRejection() ;
  
  //PID cuts (shower shape)
  
  AliCaloPID* caloPID = anaphoton->GetCaloPID();
  anaphoton->SwitchOnCaloPID(); // if nothing else specified bayesian
  anaphoton->SwitchOnCaloPIDRecalculation(); // off, get bayesian weights, on use simple cut
  //caloPID->SetLambda0CutMax(0.30);
  //caloPID->SetLambda0CutMin(0.10);
  //anaphoton->SwitchOffFillShowerShapeHistograms();  // Filled before photon shower shape selection
  
  // Input / output delta AOD settings
  
  if(!kData.Contains("delta")) {
    anaphoton->SetOutputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
    anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    //anaphoton->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
  }
  else anaphoton->SetInputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  
  //Set Histograms name tag, bins and ranges
  
  anaphoton->AddToHistogramsName("AnaPhoton_");
  SetHistoRangeAndNBins(anaphoton, kCalorimeter, kYears, kCollisions, kSimulation); // see method below
  
  // Number of particle type MC histograms
  //anaphoton->FillNOriginHistograms(8);
  //anaphoton->FillNPrimaryHistograms(4);
  
  if(kPrint) anaphoton->Print("");
  
  return anaphoton;
  
}



//_______________________________
AliAnaPi0* ConfigurePi0Analysis(TString kData, TString kCalorimeter, Int_t kYears, TString kCollisions, TString kClusterArray, 
                                Bool_t kSimulation, TString kTrig = "", Bool_t kPrint)
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
    printf("****** Configure Pi0 for pp analysis\n");
    anapi0->SetNCentrBin(1);
    anapi0->SetNZvertBin(10);
    anapi0->SetNRPBin(1);
    anapi0->SetNMaxEvMix(100);    
  }
  else if(kCollisions=="PbPb") {
    printf("****** Configure Pi0 for PbPb analysis\n");
    anapi0->SetNCentrBin(5);
    anapi0->SetNZvertBin(3);
    anapi0->SetNRPBin(1);
    anapi0->SetNMaxEvMix(5);
  }
  
  anapi0->SwitchOffMultipleCutAnalysis();
  anapi0->SwitchOffSMCombinations();

  //Set Histograms name tag, bins and ranges
  
  anapi0->AddToHistogramsName("AnaPi0_");
  SetHistoRangeAndNBins(anapi0, kCalorimeter, kYears, kCollisions, kSimulation); // see method below
  
  return anapi0;
  
}

//________________________________________________________
void SetHistoRangeAndNBins (AliAnaPartCorrBaseClass* ana, TString kCalorimeter, 
                            Int_t kYears, TString kCollisions, Bool_t kSimulation)
{
  // Set common bins for all analysis and MC histograms filling
  
  if(kSimulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else            ana->SwitchOffDataMC() ;
  
  ana->SetHistoPtRangeAndNBins(0, 100, 250) ; // Energy and pt histograms
  
  if(kCalorimeter=="EMCAL"){
    if(kYears==2010){
      ana->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
      //ana->SetHistoXRangeAndNBins(-230,90,120); // QA
      //ana->SetHistoYRangeAndNBins(370,450,40);  // QA
    }
    else {           
      ana->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
      //ana->SetHistoXRangeAndNBins(-600,90,200); // QA
      //ana->SetHistoYRangeAndNBins(100,450,100); // QA
    }
    
    ana->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  }
  else{
    ana->SetHistoPhiRangeAndNBins(260*TMath::DegToRad(), 320*TMath::DegToRad(), 60) ;
    ana->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
    
  }
  
  //ana->SetHistoShowerShapeRangeAndNBins(0, 3, 300);
  
  // Invariant mass analysis
  ana->SetHistoMassRangeAndNBins(0., 1., 200) ;
  ana->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  
  // check if time calibration is on
  //ana->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  //ana->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // QA, electron, charged
  //ana->SetHistoPOverERangeAndNBins(0,10.,100);
  //ana->SetHistodEdxRangeAndNBins(0.,200.,200);
  
  // QA
  //ana->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  //ana->SetHistodRRangeAndNBins(0.,TMath::Pi(),150);
  //ana->SetHistoRatioRangeAndNBins(0.,2.,100);
  //ana->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  //ana->SetHistoNClusterCellRangeAndNBins(0,500,500);
  //ana->SetHistoZRangeAndNBins(-400,400,200);
  //ana->SetHistoRRangeAndNBins(400,450,25);
  //ana->SetHistoV0SignalRangeAndNBins(0,5000,500);
  //ana->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  //ana->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
  
}

