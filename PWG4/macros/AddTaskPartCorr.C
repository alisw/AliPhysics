AliAnalysisTaskParticleCorrelation *AddTaskPartCorr(TString data, TString calorimeter)
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
  
  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  reader->SetCTSPtMin(0.2);
  reader->Print("");
  
  // ##### Analysis algorithm settings ####

  // --------------------
  // --- Pi0 Analysis ---
  // --------------------
  
  AliCaloPID * pid = new AliCaloPID();
  pid->SetDispersionCut(1.5);
  pid->SetTOFCut(5.e-9);
  pid->SetDebug(-1);
  pid->Print("");
  
  AliAnaPhoton *anaphoton1 = new AliAnaPhoton();
  anaphoton1->SetDebug(-1); //10 for lots of messages
  //anaphoton->SetMinPt(0.5);
  anaphoton1->SetMinDistanceToBadChannel(2, 4, 5);
  anaphoton1->SetCaloPID(pid);
  anaphoton1->SetCalorimeter(calorimeter);
  anaphoton1->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anaphoton1->SwitchOffCaloPID();
  anaphoton1->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
  anaphoton1->SwitchOffFidutialCut();
  anaphoton1->SetOutputAODName("PhotonsForPi0IM"+calorimeter);
  anaphoton1->Print("");

  AliAnaPi0 *anapi0 = new AliAnaPi0();
  anapi0->SetDebug(-1);//10 for lots of messages
  anapi0->SetInputAODName("PhotonsForPi0IM"+calorimeter);
  anapi0->SetCaloPID(pid);
  anapi0->SetCalorimeter(calorimeter);
  anapi0->SwitchOnFidutialCut();
  anapi0->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anapi0->Print("");
  
  
  // -------------------------------------------------
  // --- Photon Isolation and Correlation Analysis ---
  // -------------------------------------------------
  
  AliAnaPhoton *anaphoton2 = new AliAnaPhoton();
  anaphoton2->SetDebug(-1); //10 for lots of messages
  anaphoton2->SetMinPt(5);
  anaphoton2->SetCaloPID(pid);
  anaphoton2->SetCalorimeter(calorimeter);
  anaphoton2->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anaphoton2->SwitchOnCaloPID();
  if(calorimeter == "EMCAL") anaphoton2->SwitchOnCaloPIDRecalculation();
  anaphoton2->SwitchOffFidutialCut();
  anaphoton2->SetOutputAODName("DirectPhotons"+calorimeter);
  anaphoton2->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  anaphoton2->AddToHistogramsName("AnaPhotonCorr_");
  anaphoton2->Print("");
  // ### Isolation analysis ###	
  
  AliIsolationCut * ic = new AliIsolationCut();
  ic->SetConeSize(0.5);
  ic->SetPtThreshold(1.);
  ic->SetICMethod(AliIsolationCut::kPtThresIC);
  ic->Print("");
  
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  //anaisol->SetMinPt(5);
  anaisol->SetInputAODName("DirectPhotons"+calorimeter);
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
  
  // ### Correlation with Jet Finder AOD output
  AliAnaParticleJetFinderCorrelation *anacorrjet = new AliAnaParticleJetFinderCorrelation();
  anacorrjet->SetInputAODName("DirectPhotons"+calorimeter);
  anacorrjet->SwitchOffFidutialCut();
  anacorrjet->SetDebug(-1);
  anacorrjet->SetConeSize(1);  
  anacorrjet->SelectIsolated(kTRUE); // do correlation with isolated photons
  anacorrjet->SetPtThresholdInCone(0.2);
  anacorrjet->SetDeltaPhiCutRange(0.5,5.5);//Mostly Open Cuts 
  anacorrjet->SetRatioCutRange(0.01,3); //Mostly Open Cuts
  anacorrjet->UseJetRefTracks(kFALSE); //Not working now
  anacorrjet->Print("");
  
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorrhadron = new AliAnaParticleHadronCorrelation();
  anacorrhadron->SetInputAODName("DirectPhotons"+calorimeter);
  anacorrhadron->SetDebug(-1);
  anacorrhadron->SwitchOffFidutialCut();
  anacorrhadron->SetPtCutRange(1,100);
  anacorrhadron->SetDeltaPhiCutRange(1.5,4.5);
  anacorrhadron->SelectIsolated(kTRUE); // do correlation with isolated photons
	if(calorimeter=="PHOS"){
    //Correlate with particles in EMCAL
    anacorrhadron->SwitchOnCaloPID();
    anacorrhadron->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  }
  anacorrhadron->Print("");
  
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(anaphoton1,0);
  maker->AddAnalysis(anapi0,1);
  maker->AddAnalysis(anaphoton2,2);
  maker->AddAnalysis(anaisol,3);
  maker->AddAnalysis(anacorrjet,4);
  maker->AddAnalysis(anacorrhadron,5);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  maker->Print("");
  
  printf("======================== \n");
  printf(" End Configuration of PartCorr analysis with detector %s \n",calorimeter.Data());
  printf("======================== \n");
  
   // Create task
   //===========================================================================
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation ("PartCorr");
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  task->SetAnalysisMaker(maker);				
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form("PartCorr"),TList::Class(),
							   AliAnalysisManager::kOutputContainer, Form("PartCorr_%s.root",calorimeter.Data()));
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  
  return task;
}


