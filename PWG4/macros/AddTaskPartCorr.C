AliAnalysisTaskParticleCorrelation *AddTaskPartCorr(TString data, TString calorimeter, Bool_t kUseKinematics = kFALSE, Bool_t kPrintSettings = kFALSE)
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
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
   //cout<<"DATA TYPE :: "<<inputDataType<<endl;
   // inputDataType: data managed by the input handler
   // data: can be same as one managed by input handler, or the output AOD created by the filter. By default use AOD
	
   // Configure analysis
   //===========================================================================
  
  //Reader
  AliCaloTrackReader * reader = 0x0;
  if(data=="AOD") reader = new AliCaloTrackAODReader();
  else if(data=="ESD") reader = new AliCaloTrackESDReader();
  else if(data=="MC" && dataType == "ESD") reader = new AliCaloTrackMCReader();
  reader->SetDebug(-1);//10 for lots of messages
  if(calorimeter == "EMCAL") reader->SwitchOnEMCALCells();
  if(calorimeter == "PHOS")  reader->SwitchOnPHOSCells();
  if(kUseKinematics){
	if(inputDataType == "ESD"){
		reader->SwitchOnStack();          
		reader->SwitchOffAODMCParticles(); 
	}
	else if(inputDataType == "AOD"){
		reader->SwitchOffStack();          
		reader->SwitchOnAODMCParticles(); 
	}
  }
	
  //Min particle pT
  reader->SetEMCALPtMin(0.2); 
  reader->SetPHOSPtMin(0.2);
  reader->SetCTSPtMin(0.2);
  if(kPrintSettings) reader->Print("");
	
  //Needed line, do not clear standard output AODs, we are not writing there
  reader->SwitchOnWriteStdAOD();
	
  // ##### Analysis algorithm settings ####

  // --------------------
  // --- Pi0 Analysis ---
  // --------------------
  
  AliCaloPID * pid = new AliCaloPID();
  pid->SetDispersionCut(1.5);
  pid->SetTOFCut(5.e-9);
  pid->SetDebug(-1);
  if(kPrintSettings) pid->Print("");
	
  AliFidutialCut * fidCut = new AliFidutialCut();
  fidCut->DoCTSFidutialCut(kFALSE) ;
  fidCut->DoEMCALFidutialCut(kTRUE) ;
  fidCut->DoPHOSFidutialCut(kTRUE) ;
	
  AliAnaCalorimeterQA *qa = new AliAnaCalorimeterQA();
  qa->SetDebug(-1); //10 for lots of messages
  qa->SetCalorimeter(calorimeter);
  if(kUseKinematics && inputDataType!="AOD") qa->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else  qa->SwitchOffDataMC() ;
  //qa->AddToHistogramsName(Form("AnaCaloQA_%s",calorimeter.Data()));
  qa->AddToHistogramsName("AnaCaloQA_");
  qa->SetFidutialCut(fidCut);
  qa->SwitchOnFidutialCut();
  if(kPrintSettings) qa->Print("");	
	
  AliFidutialCut * fidCut1stYear = new AliFidutialCut();
  fidCut1stYear->DoCTSFidutialCut(kFALSE) ;
  fidCut1stYear->DoEMCALFidutialCut(kTRUE) ;
  fidCut1stYear->DoPHOSFidutialCut(kTRUE) ;
  fidCut1stYear->SetSimpleEMCALFidutialCut(0.7,80.,120.);
  fidCut1stYear->SetSimplePHOSFidutialCut(0.12,260.,320.);
	
  AliAnaPhoton *anaphoton1 = new AliAnaPhoton();
  anaphoton1->SetDebug(-1); //10 for lots of messages
  //anaphoton->SetMinPt(0.5);
  anaphoton1->SetMinDistanceToBadChannel(2, 4, 5);
  anaphoton1->SetCaloPID(pid);
  anaphoton1->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaphoton1->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaphoton1->SwitchOffDataMC() ;
  anaphoton1->SwitchOffCaloPID();
  anaphoton1->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
  anaphoton1->SwitchOnFidutialCut();
  anaphoton1->SetFidutialCut(fidCut1stYear);
  anaphoton1->SetOutputAODName(Form("PhotonsForIM%s",calorimeter.Data()));
  if(kPrintSettings) anaphoton1->Print("");

  AliAnaPi0 *anapi0 = new AliAnaPi0();
  anapi0->SetDebug(-1);//10 for lots of messages
  anapi0->SetInputAODName(Form("PhotonsForIM%s",calorimeter.Data()));
  anapi0->SetCaloPID(pid);
  anapi0->SetCalorimeter(calorimeter);
  anapi0->SwitchOnFidutialCut();
  anapi0->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  if(kPrintSettings) anapi0->Print("");
  
	
  // -------------------------------------------------
  // --- Photon Isolation and Correlation Analysis ---
  // -------------------------------------------------
  
  AliAnaPhoton *anaphoton2 = new AliAnaPhoton();
  anaphoton2->SetDebug(-1); //10 for lots of messages
  anaphoton2->SetMinPt(5);
  anaphoton2->SetCaloPID(pid);
  anaphoton2->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaphoton2->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaphoton2->SwitchOffDataMC() ;
  anaphoton2->SwitchOnCaloPID();
  if(calorimeter == "EMCAL") anaphoton2->SwitchOnCaloPIDRecalculation();
  anaphoton2->SwitchOffFidutialCut();
  anaphoton2->SetOutputAODName(Form("DirectPhotons%s",calorimeter.Data()));
  anaphoton2->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  anaphoton2->AddToHistogramsName("AnaPhotonCorr_");
  if(kPrintSettings) anaphoton2->Print("");
  // ### Isolation analysis ###	
  
  AliIsolationCut * ic = new AliIsolationCut();
  ic->SetConeSize(0.5);
  ic->SetPtThreshold(1.);
  ic->SetICMethod(AliIsolationCut::kPtThresIC);
  if(kPrintSettings) ic->Print("");
  
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  //anaisol->SetMinPt(5);
  anaisol->SetInputAODName(Form("DirectPhotons%s",calorimeter.Data()));
  anaisol->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaisol->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaisol->SwitchOffDataMC() ;
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
  if(kPrintSettings) anaisol->Print("");
  
  // ### Correlation with Jet Finder AOD output
  AliAnaParticleJetFinderCorrelation *anacorrjet = new AliAnaParticleJetFinderCorrelation();
  anacorrjet->SetInputAODName(Form("DirectPhotons%s",calorimeter.Data()));
  anacorrjet->SwitchOffFidutialCut();
  anacorrjet->SetDebug(-1);
  anacorrjet->SetConeSize(1);  
  anacorrjet->SelectIsolated(kTRUE); // do correlation with isolated photons
  anacorrjet->SetPtThresholdInCone(0.2);
  anacorrjet->SetDeltaPhiCutRange(0.5,5.5);//Mostly Open Cuts 
  anacorrjet->SetRatioCutRange(0.01,3); //Mostly Open Cuts
  anacorrjet->UseJetRefTracks(kFALSE); //Not working now
  if(kPrintSettings) anacorrjet->Print("");
  
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorrhadron = new AliAnaParticleHadronCorrelation();
  anacorrhadron->SetInputAODName(Form("DirectPhotons%s",calorimeter.Data()));
  anacorrhadron->SetOutputAODName(Form("CorrelatedPi0s%s",calorimeter.Data()));
  anacorrhadron->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  anacorrhadron->SetDebug(-1);
  anacorrhadron->SwitchOffCaloPID();
  anacorrhadron->SwitchOffFidutialCut();
  anacorrhadron->SetPtCutRange(1,100);
  anacorrhadron->SetDeltaPhiCutRange(1.5,4.5);
  anacorrhadron->SelectIsolated(kTRUE); // do correlation with isolated photons
  if(calorimeter=="PHOS"){
    //Correlate with particles in EMCAL
    anacorrhadron->SwitchOnCaloPID();
    anacorrhadron->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  }
  if(kPrintSettings) anacorrhadron->Print("");
  
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(qa,0);
  maker->AddAnalysis(anaphoton1,1);
  maker->AddAnalysis(anapi0,2);
  maker->AddAnalysis(anaphoton2,3);
  maker->AddAnalysis(anaisol,4);
  maker->AddAnalysis(anacorrjet,5);
  maker->AddAnalysis(anacorrhadron,6);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  if(kPrintSettings) maker->Print("");
  
  printf("======================== \n");
  printf(" End Configuration of PartCorr analysis with detector %s \n",calorimeter.Data());
  printf("======================== \n");
  
   // Create task
   //===========================================================================
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation (Form("PartCorr%s",calorimeter.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  task->SetAnalysisMaker(maker);				
  mgr->AddTask(task);
  
  char name[128];
  sprintf(name,"PartCorr_%s",calorimeter.Data());
  cout<<"Name of task "<<name<<endl;
  //AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form(name),TList::Class(),
  	//					   AliAnalysisManager::kOutputContainer, Form("PartCorr_%s.root",calorimeter.Data()));
  
  TString outputfile = AliAnalysisManager::GetCommonFileName(); 
  //  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form("PartCorr_%s",calorimeter.Data()),  TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PartCorr_%s",outputfile.Data(),calorimeter.Data()));
  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(calorimeter.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PartCorr",outputfile.Data()));

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  
  return task;
}


