AliAnalysisTaskParticleCorrelation *AddTaskPartCorr(Char_t * analysis, TString data, TString calorimeter)
{
  // Creates a PartCorr task, configures it and adds it to the analysis manager.
  
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
   //TString dataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
   // Configure analysis
   //===========================================================================
  
  //Reader
  AliCaloTrackReader * reader = 0x0;
  if(data=="AOD") reader = new AliCaloTrackAODReader();
  else if(data=="ESD") reader = new AliCaloTrackESDReader();
  else if(data=="MC") reader = new AliCaloTrackMCReader();
  reader->SetDebug(-1);//10 for lots of messages
  
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  switch (analysis) {
  case "Pi0":
    AnalysisPi0(maker,reader,calorimeter);
    break;
  case "GammaJetFinder":
    //JETAN must run before
    AnalysisGammaJetFinderCorrelation(maker,reader,calorimeter);		
    break;	
  case "GammaHadron":
    AnalysisGammaHadronCorrelation(maker,reader,calorimeter);		
    break;		
  }
  
  
   // Create task
   //===========================================================================
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation (analysis);
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetAnalysisMaker(maker);				
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form("%s", analysis),TList::Class(),
							   AliAnalysisManager::kOutputContainer, Form("%s_%s.root",analysis,calorimeter.Data()));
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  
  return task;
}

//_________________________________________________________________________________________________
void AnalysisPi0(AliAnaPartCorrMaker * maker, AliCaloTrackReader *  reader, TString calorimeter){
  //Configuration for pi0 invariant mass analysis
  
  // #### Reader ####
  //Switch on or off the detectors 
  if(calorimeter=="EMCAL"){
    reader->SwitchOnEMCAL();
    reader->SwitchOffPHOS();
  }
  else if(calorimeter=="PHOS"){
    reader->SwitchOffEMCAL();
    reader->SwitchOnPHOS();
  }
  else if {
    printf("ABORT analysis: Wrong calorimeter in configuration: %s\n",calorimeter.Data());
    abort();
  }	
  reader->SwitchOffCTS();
  
  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  reader->Print("");
  
  // ##### Analysis algorithm ####
  
  AliCaloPID * pid = new AliCaloPID();
  pid->SetDispersionCut(1.5);
  pid->SetTOFCut(5.e-9);
  pid->SetDebug(-1);
  pid->Print("");
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  //anaphoton->SetMinPt(0.5);
  anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
  anaphoton->SetCaloPID(pid);
  anaphoton->SetCalorimeter(calorimeter);
  anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anaphoton->SwitchOffCaloPID();
  anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
  anaphoton->SwitchOffFidutialCut();
  anaphoton->SetOutputAODName("Photons"+calorimeter);
  anaphoton->SetOutputAODClassName("AliAODPWG4Particle");
  
  AliAnaPi0 *anapi0 = new AliAnaPi0();
  anapi0->SetDebug(-1);//10 for lots of messages
  anapi0->SetInputAODName("Photons"+calorimeter);
  anapi0->SetCaloPID(pid);
  anapi0->SetCalorimeter(calorimeter);
  anapi0->SwitchOnFidutialCut();
  anapi0->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anapi0->Print("");
  
  // #### Configure Maker ####
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(anaphoton,0);
  maker->AddAnalysis(anapi0,1);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  maker->Print("");
  //
  printf("======================== \n");
  printf("End Configuration of AnalysisPi0() \n");
  printf("======================== \n");
  
  
}
//_________________________________________________________________________________________________
void AnalysisGammaJetFinderCorrelation(AliAnaPartCorrMaker * maker, AliCaloTrackReader *  reader, TString calorimeter){
  //Configuration for pi0 invariant mass analysis
  
  // #### Reader ####
  //Switch on or off the detectors 
  if(calorimeter=="PHOS")
    reader->SwitchOnPHOS();
  else if(calorimeter=="EMCAL")	
    reader->SwitchOffPHOS();
  else if {
    printf("ABORT analysis: Wrong calorimeter in configuration: %s\n",calorimeter.Data());
    abort();
  }	
  reader->SwitchOnEMCAL();
  reader->SwitchOnCTS();
  
  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  reader->SetCTSPtMin(0.1);
  reader->Print("");
  
  // ##### Analysis algorithm ####
  // ### Photon analysis ###
  //AliCaloPID * pid = new AliCaloPID();
  //pid->Print("");
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  anaphoton->SetMinPt(5);
  //anaphoton->SetCaloPID(pid);
  anaphoton->SetCalorimeter(calorimeter);
  anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anaphoton->SwitchOnCaloPID();
  if(calorimeter == "EMCAL") anaphoton->SwitchOnCaloPIDRecalculation();
  anaphoton->SwitchOffFidutialCut();
  anaphoton->SetOutputAODName("DirectPhotonsJet"+calorimeter);
  anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  
  // ### Isolation analysis ###	
  
  AliIsolationCut * ic = new AliIsolationCut();
  ic->SetConeSize(0.5);
  ic->SetPtThreshold(1.);
  ic->SetICMethod(AliIsolationCut::kPtThresIC);
  ic->Print("");
  
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  //anaisol->SetMinPt(5);
  anaisol->SetInputAODName("DirectPhotonsJet"+calorimeter);
  anaisol->SetCalorimeter(calorimeter);
  anaisol->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  //Select clusters with no pair, if both clusters with pi0 mass
  anaisol->SwitchOffInvariantMass();
  //anaisol->SetNeutralMesonSelection(nms);
  //Do isolation cut
  anaisol->SetIsolationCut(ic);	
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisol->SwitchOffReIsolation();
  //Multiple IC
  anaisol->SwitchOffSeveralIsolation() ;
  anaisol->Print("");
  
  // ### Correlatio with Jet Finder AOD output
  AliAnaParticleJetFinderCorrelation *anacorr = new AliAnaParticleJetFinderCorrelation();
  anacorr->SetInputAODName("DirectPhotonsJet"+calorimeter);
  anacorr->SwitchOffFidutialCut();
  anacorr->SetDebug(-1);
  anacorr->SetConeSize(1);  
  anacorr->SelectIsolated(kTRUE); // do correlation with isolated photons
  anacorr->SetPtThresholdInCone(0.2);
  anacorr->SetDeltaPhiCutRange(0.5,5.5);//Mostly Open Cuts 
  anacorr->SetRatioCutRange(0.01,3); //Mostly Open Cuts
  anacorr->UseJetRefTracks(kFALSE); //Not working now
  anacorr->Print("");
  
  // #### Configure Maker ####
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(anaphoton,0);
  maker->AddAnalysis(anaisol,1);
  maker->AddAnalysis(anacorr,2);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  maker->Print("");
  //
  printf("======================== \n");
  printf("End Configuration of AnalysisGammaJetFinderCorrelation() \n");
  printf("======================== \n");
  
  
}


//_________________________________________________________________________________________________
void AnalysisGammaHadronCorrelation(AliAnaPartCorrMaker * maker, AliCaloTrackReader *  reader, TString calorimeter){
  //Configuration for pi0 invariant mass analysis
  
  // #### Reader ####
  //Switch on or off the detectors 
  if(calorimeter=="PHOS")
    reader->SwitchOnPHOS();
  else if(calorimeter=="EMCAL")	
    reader->SwitchOffPHOS();
  else if {
    printf("ABORT analysis: Wrong calorimeter in configuration: %s\n",calorimeter.Data());
    abort();
  }	
  reader->SwitchOnEMCAL();
  reader->SwitchOnCTS();
  
  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  reader->SetCTSPtMin(0.1);
  reader->Print("");
  
  // ##### Analysis algorithm ####
  // ### Photon analysis ###
  //AliCaloPID * pid = new AliCaloPID();
  //pid->Print("");
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  anaphoton->SetMinPt(5);
  //anaphoton->SetCaloPID(pid);
  anaphoton->SetCalorimeter(calorimeter);
  anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anaphoton->SwitchOnCaloPID();
  if(calorimeter == "EMCAL") anaphoton->SwitchOnCaloPIDRecalculation();
  anaphoton->SwitchOffFidutialCut();
  anaphoton->SetOutputAODName("DirectPhotonsHadron"+calorimeter);
  anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  
  // ### Isolation analysis ###	
  
  AliIsolationCut * ic = new AliIsolationCut();
  ic->SetConeSize(0.5);
  ic->SetPtThreshold(1.);
  ic->SetICMethod(AliIsolationCut::kPtThresIC);
  ic->Print("");
  
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  //anaisol->SetMinPt(5);
  anaisol->SetInputAODName("DirectPhotonsHadron"+calorimeter);
  anaisol->SetCalorimeter(calorimeter);
  anaisol->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  //Select clusters with no pair, if both clusters with pi0 mass
  anaisol->SwitchOffInvariantMass();
  //anaisol->SetNeutralMesonSelection(nms);
  //Do isolation cut
  anaisol->SetIsolationCut(ic);	
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisol->SwitchOffReIsolation();
  //Multiple IC
  anaisol->SwitchOffSeveralIsolation() ;
  anaisol->Print("");
  
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorr = new AliAnaParticleHadronCorrelation();
  anacorr->SetInputAODName("DirectPhotonsHadron"+calorimeter);
  anacorr->SetDebug(-1);
  anacorr->SwitchOffFidutialCut();
  anacorr->SetPtCutRange(1,100);
  anacorr->SetDeltaPhiCutRange(1.5,4.5);
  if(calorimeter=="PHOS"){
    //Correlate with particles in EMCAL
    anacorr->SwitchOnCaloPID();
    anacorr->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  }
  anacorr->Print("");
  
  // #### Configure Maker ####
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(anaphoton,0);
  maker->AddAnalysis(anaisol,1);
  maker->AddAnalysis(anacorr,2);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  maker->Print("");
  //
  printf("======================== \n");
  printf("End Configuration of AnalysisGammaJetFinderCorrelation() \n");
  printf("======================== \n");
  
}
