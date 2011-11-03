Bool_t  kPrint = kFALSE;
Bool_t  kSimulation    = kFALSE;
Bool_t  kUseKinematics = kFALSE;
Bool_t  kOutputAOD     = kFALSE;
Int_t   kRunNumber     = 0;
Int_t   kYears         = 2011;
TString kCollisions     = "pp";
TString kTrig          ="EMC7" ;
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
  
  kPrint = printSettings;
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
  
  maker->AddAnalysis(ConfigureQAAnalysis()         , n++); // QA analysis (careful, if no tender used, cells are not updated or corrected, clusters OK)
  maker->AddAnalysis(ConfigurePhotonAnalysis()     , n++); // Photon cluster selection
  maker->AddAnalysis(ConfigurePi0Analysis()        , n++); // Pi0 invariant mass analysis
  
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kIMCalo)      , n++); // Pi0 event by event selection, and photon tagging from decay
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Eta", AliAnaPi0EbE::kIMCalo)      , n++); // Eta event by event selection, and photon tagging from decay

  //maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kIMCaloTracks), n++); // Pi0 (calo+conversion) event by event selection, 
  // and photon tagging from decay, need to execute at the same time conversions analysis

  maker->AddAnalysis(ConfigureIsolationAnalysis("Photon"), n++);   // Photon isolation
  maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0")   , n++); // Pi0 isolation
  
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kFALSE), n++);   // Gamma hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kTRUE) , n++);   // Isolated gamma hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0",kFALSE)   , n++); // Pi0 hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0",kTRUE)    , n++); // Isolated pi0 hadron correlation
  
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(kData.Contains("delta")) maker->SwitchOffAODsMaker() ;
  else                        maker->SwitchOnAODsMaker()  ;
	
  if(kPrint) maker->Print("");
  
  printf("<< End Configuration of analysis %d for calorimeter %s >>\n",n, kCalorimeter.Data());
  
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
  if(kTrig=="EMC1"){
    printf("PartCorr trigger EMC1\n");
    task->SelectCollisionCandidates(AliVEvent::kEMC1);
  }
  else if(kTrig=="MB"){
    printf("PartCorr trigger MB\n");
    task->SelectCollisionCandidates(AliVEvent::kMB);
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
  
  reader->SetEMCALClusterListName(kClusterArray);
  if(kClusterArray == "") {
    printf("**************** Normal analysis **************** \n");
    //reader->SwitchOnClusterRecalculation();
  }
  else {
    printf("**************** Input for analysis is Clusterizer %s **************** \n", kClusterArray.Data());
    //reader->SwitchOffClusterRecalculation();
  }
  
  // Settings for LHC11a
  if(kRunNumber > 140000 && kRunNumber < = 146860){
    if(kClusterArray == "") reader->SwitchOnLEDEventsRemoval();
    reader->RejectFastClusterEvents();
    printf("Reader: Reject LED events and Fast cluster\n");
  }
  
  //reader->SetDeltaAODFileName("");
  //if(!kSimulation) reader->SetFiredTriggerClassName("CINT1B-ABCE-NOPF-ALL");
  
  // Detector input filling
  
  reader->SwitchOnCTS();
  
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
  
  //Event selection
  if     (kCollisions=="pp"  ) {
    if(kRunNumber < 140000) reader->SwitchOnEventSelection(); // remove pileup by default
    reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
    reader->SwitchOffPrimaryVertexSelection(); // and besides primary vertex
  }
  else if(kCollisions=="PbPb") {
    reader->SwitchOffEventSelection();         // remove pileup by default
    reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
    reader->SwitchOffPrimaryVertexSelection(); // and besides primary vertex
  }
  
  if     (kCollisions=="pp"  )   reader->SetZvertexCut(50.);  //Open cut
  else if(kCollisions=="PbPb")   reader->SetZvertexCut(10.);  //Centrality defined in this range.
  
  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  reader->SetEMCALPtMax(1000); 
  reader->SetPHOSPtMin(0.3);
  reader->SetCTSPtMin(0.);
  
  if(kOutputAOD)     reader->SwitchOnWriteDeltaAOD()  ;
  if(kPrint) reader->Print("");
  
  return reader;
  
}

//_______________________________________
AliCalorimeterUtils* ConfigureCaloUtils()
{
  
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  
  if(kYears==2010) cu->SetEMCALGeometryName("EMCAL_FIRSTYEARV1");
  else            cu->SetEMCALGeometryName("EMCAL_COMPLETEV1");
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  if(kClusterArray == "") 
    cu->SwitchOffRecalculateClusterTrackMatching(); // Done in clusterization
  else            
    cu->SwitchOnRecalculateClusterTrackMatching();
  
  //EMCAL only settings
  AliEMCALRecoUtils * reco = cu->GetEMCALRecoUtils();
  
  if(kCollisions=="pp"  ) {
    cu->SwitchOnCorrectClusterLinearity();
    if(!kSimulation) reco->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrected);
    else             reco->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MC);
  }
  
  if(kClusterArray == "") reco->SwitchOnRejectExoticCluster();
  else                    reco->SwitchOffRejectExoticCluster();
  
  cu->SetDebug(-1);
  if(kPrint) cu->Print("");
  
  return cu;
  
}

//_____________________________________
AliAnaPhoton* ConfigurePhotonAnalysis()
{
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  anaphoton->SwitchOnFillShowerShapeHistograms();  
  anaphoton->SwitchOnTrackMatchRejection() ;
  
  if(kCalorimeter == "PHOS"){
    anaphoton->SetNCellCut(2);// At least 2 cells
    anaphoton->SetMinPt(0.3);
    anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
  }
  else {//EMCAL
    anaphoton->SetNCellCut(1);// At least 2 cells
    anaphoton->SetMinPt(0.5); // no effect minium EMCAL cut.
    anaphoton->SetMaxPt(1000); 
    //if(!kUseKinematics) anaphoton->SetTimeCut(400,900);// Time window of [400-900] ns
    anaphoton->SetMinDistanceToBadChannel(1, 2, 3);//For filtered AODs, new releases.
  }
  
  anaphoton->SetCalorimeter(kCalorimeter);
  
  if(kUseKinematics) anaphoton->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else               anaphoton->SwitchOffDataMC() ;
  
  //PID cut
  AliCaloPID* caloPID = anaphoton->GetCaloPID();
  caloPID->SetLambda0CutMax(0.4);
  caloPID->SetLambda0CutMin(0.1);
  anaphoton->SwitchOnCaloPID();
  anaphoton->SwitchOnCaloPIDRecalculation(); // off, get bayesian weights, on use simple cut
  
  anaphoton->SwitchOffFiducialCut();
  
  if(!kData.Contains("delta")) {
    anaphoton->SetOutputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
    anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else anaphoton->SetInputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  
  anaphoton->AddToHistogramsName("AnaPhoton_");
  
  //Set Histograms bins and ranges
  anaphoton->SetHistoPtRangeAndNBins(0, 100, 250) ;
  if(kYears==2010)anaphoton->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
  else           anaphoton->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
  anaphoton->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  anaphoton->SetHistoShowerShapeRangeAndNBins(0, 1.5, 150);
  
  if(kPrint) anaphoton->Print("");
  
  return anaphoton;
  
}

//_______________________________
AliAnaPi0* ConfigurePi0Analysis()
{
  
  AliAnaPi0 *anapi0 = new AliAnaPi0();
  anapi0->SetDebug(-1);//10 for lots of messages
  anapi0->SetInputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  anapi0->SetCalorimeter(kCalorimeter);
  anapi0->SwitchOffMultipleCutAnalysis(); 
  anapi0->SetPairTimeCut(70);
  
  //settings for pp collision mixing
  anapi0->SwitchOnOwnMix(); //Off when mixing done with general mixing frame
  
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
  
  if(kUseKinematics)anapi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else              anapi0->SwitchOffDataMC() ;
  if(kCalorimeter=="PHOS") anapi0->SetNumberOfModules(3); //PHOS first year
  else {                   
    if(kYears==2010) anapi0->SetNumberOfModules(4); //EMCAL first year
    else            anapi0->SetNumberOfModules(10);
  }
  
  anapi0->SetHistoPtRangeAndNBins(0, 100, 250) ;    
  anapi0->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  if(kYears==2010)anapi0->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
  else           anapi0->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
  
  anapi0->SetHistoMassRangeAndNBins(0., 1., 200) ;
  anapi0->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  anapi0->SetHistoTrackMultiplicityRangeAndNBins(0, 200, 20); 
  anapi0->SetHistoShowerShapeRangeAndNBins(0, 3, 300);
  anapi0->SetHistoTimeRangeAndNBins(-200.,200,800);
  anapi0->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  if(kPrint) anapi0->Print("");
  
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
  anapi0ebe->SetPairTimeCut(20); 
  anapi0ebe->SetCalorimeter(kCalorimeter);
  
  anapi0ebe->AddToHistogramsName(Form("Ana%s%sEbE_",particle.Data(),opt.Data()));
  anapi0ebe->SetInputAODName(Form("Photon%s_Trig%s_Cl%s",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  if(!kInputDataType.Contains("delta")) {
    anapi0ebe->SetOutputAODName(Form("%s%s%s_Trig%s_Cl%s",particle.Data(), opt.Data(), kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
    anapi0ebe->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else  anapi0ebe->SetInputAODName(Form("%s%s%s_Trig%s_Cl%s",particle.Data(),opt.Data(),kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  
  if(analysis == AliAnaPi0EbE::kIMCaloTracks) anapi0ebe->SetInputAODGammaConvName("PhotonsCTS");
  
  if(kUseKinematics) anapi0ebe->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else               anapi0ebe->SwitchOffDataMC() ;  
  
  if(analysis!=AliAnaPi0EbE::kSSCalo){
    
    AliNeutralMesonSelection *nms = anapi0ebe->GetNeutralMesonSelection();
    nms->SetParticle(particle);
    nms->SwitchOnAngleSelection();
    nms->KeepNeutralMesonSelectionHistos(kTRUE);
    //nms->SetAngleMaxParam(2,0.2);
    nms->SetHistoERangeAndNBins(0, 20, 100) ;
    //nms->SetHistoIMRangeAndNBins(0, 1, 400);
  }
  
  //Set Histrograms bins and ranges
  if(kCalorimeter=="EMCAL" ){
    anapi0ebe->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
    if(kYears==2010)anapi0ebe->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
    else            anapi0ebe->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
    anapi0ebe->SetHistoShowerShapeRangeAndNBins(0, 3, 300);
  } 
  
  anapi0ebe->SetHistoTimeRangeAndNBins(-200.,200,800);
  anapi0ebe->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  if(kPrint) anapi0ebe->Print("");
  
  return  anapi0ebe;
  
}

//___________________________________________________________________
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle)
{
  
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  anaisol->SetMinPt(5);
  anaisol->SetInputAODName(Form("%s%s_Trig%s_Cl%s",particle.Data(),kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  anaisol->SetAODObjArrayName(Form("IC%s",particle.Data())); 
  
  anaisol->SetCalorimeter(kCalorimeter);
  
  if(kUseKinematics) anaisol->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else               anaisol->SwitchOffDataMC() ;
  
  //Do isolation cut
  AliIsolationCut * ic =  anaisol->GetIsolationCut();	
  ic->SetConeSize(0.4);
  ic->SetPtThreshold(0.7);
  ic->SetPtFraction(0.1);
  ic->SetSumPtThreshold(1.0) ;
  ic->SetParticleTypeInCone(AliIsolationCut::kOnlyCharged);
  ic->SetICMethod(AliIsolationCut::kSumPtFracIC);
  if(kPrint) ic->Print("");
  
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisol->SwitchOffReIsolation();
  //Multiple IC
  anaisol->SwitchOffSeveralIsolation() ;
  //Set Histograms bins and ranges
  anaisol->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  anaisol->AddToHistogramsName(Form("AnaIsol%s_",particle.Data()));
  if(kPrint) anaisol->Print("");
  
  return anaisol;
  
}

//___________________________________________________________________________________
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle, 
                                                                    Int_t bIsolated)
{
  
  AliAnaParticleHadronCorrelation *anacorrhadron = new AliAnaParticleHadronCorrelation();
  anacorrhadron->SetInputAODName(Form("%s%s_Trig%s_Cl%s",particle.Data(),kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  anacorrhadron->AddToHistogramsName(Form("Ana%sHadronCorrIso%d_",particle.Data(),bIsolated));
  anacorrhadron->SetAODObjArrayName(Form("%sHadronCorrIso%d",particle.Data(),bIsolated)); 
  anacorrhadron->SetDebug(-1);
  //  if(kSimulation){
  //    anacorrhadron->SwitchOnFiducialCut();
  //    AliFiducialCut * fidCut1stYear = anacorrhadron->GetFiducialCut();
  //    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
  //    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
  //    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
  //    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
  //    fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
  //    fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
  //  }
  
  anacorrhadron->SelectIsolated(bIsolated); // do correlation with isolated photons
  
  anacorrhadron->SwitchOnDecayCorr();
  anacorrhadron->SetMultiBin(1);
  anacorrhadron->SwitchOffNeutralCorr();
  anacorrhadron->SwitchOffEventSelection();
  anacorrhadron->SetPtCutRange(0.1,100);
  anacorrhadron->SetDeltaPhiCutRange(1.5,4.5);
  anacorrhadron->SwitchOnSeveralUECalculation();
  anacorrhadron->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  anacorrhadron->SelectIsolated(kFALSE); // do correlation with isolated photons
  if(kUseKinematics) anacorrhadron->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else               anacorrhadron->SwitchOffDataMC() ;
  //if(kCalorimeter=="PHOS"){
  //Correlate with particles in EMCAL
  //anacorrhadron->SwitchOnCaloPID();
  //anacorrhadron->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  //}
  //Set Histograms bins and ranges
  anacorrhadron->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrint) anacorrhadron->Print("");  
  
  return anacorrhadron;
  
}

//________________________________________
AliAnaCalorimeterQA* ConfigureQAAnalysis()
{
  
  AliAnaCalorimeterQA *anaQA = new AliAnaCalorimeterQA();
  //anaQA->SetDebug(10); //10 for lots of messages
  anaQA->SetCalorimeter(kCalorimeter);
  
  if(kSimulation) anaQA->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else            anaQA->SwitchOffDataMC() ;
  
  anaQA->AddToHistogramsName("QA_"); //Begining of histograms name
  
  anaQA->SwitchOffFiducialCut();
  anaQA->SwitchOnCorrelation();
  anaQA->SwitchOffFillAllTH3Histogram();
  anaQA->SwitchOffFillAllPositionHistogram();
  
  //anaQA->SwitchOnStudyBadClusters() ;
  //anaQA->SwitchOnStudyClustersAsymmetry();
  //anaQA->SwitchOnStudyWeight();
  anaQA->SwitchOffFillAllTrackMatchingHistogram();
  
  //Set Histrograms bins and ranges
  anaQA->SetHistoPtRangeAndNBins(0, 100, 200) ;
  anaQA->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  
  if(kCalorimeter=="EMCAL"){
    anaQA->SetHistoEtaRangeAndNBins(-0.71, 0.71, 200) ;
    if(kYears==2010){  
      anaQA->SetNumberOfModules(4); 
      anaQA->SetHistoPhiRangeAndNBins(79*TMath::DegToRad(), 121*TMath::DegToRad(), 200) ;
      anaQA->SetHistoXRangeAndNBins(-230,90,120);
      anaQA->SetHistoYRangeAndNBins(370,450,40);
    }
    else{            
      anaQA->SetNumberOfModules(10); 
      anaQA->SetHistoPhiRangeAndNBins(79*TMath::DegToRad(), 181*TMath::DegToRad(), 200) ;
      anaQA->SetHistoXRangeAndNBins(-600,90,200);
      anaQA->SetHistoYRangeAndNBins(100,450,100);
    }
  }
  anaQA->SetHistoMassRangeAndNBins(0., 1, 400) ;
  anaQA->SetHistoAsymmetryRangeAndNBins(0., 1. , 10 );
  anaQA->SetHistoPOverERangeAndNBins(0,10.,100);
  anaQA->SetHistodEdxRangeAndNBins(0.,200.,200);
  anaQA->SetHistodRRangeAndNBins(0.,TMath::Pi(),150);
  anaQA->SetHistoTimeRangeAndNBins(-500.,500,1000);
  anaQA->SetHistoRatioRangeAndNBins(0.,2.,100);
  anaQA->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  anaQA->SetHistoNClusterCellRangeAndNBins(0,500,500);
  anaQA->SetHistoZRangeAndNBins(-400,400,200);
  anaQA->SetHistoRRangeAndNBins(400,450,25);
  anaQA->SetHistoV0SignalRangeAndNBins(0,5000,500);
  anaQA->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  anaQA->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
	
  if(kPrint) anaQA->Print("");	
  
  return anaQA;
  
}


