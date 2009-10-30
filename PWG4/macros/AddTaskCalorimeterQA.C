AliAnalysisTaskParticleCorrelation *AddTaskCalorimeterQA(TString data)
{
  // Creates a PartCorr task for calorimeters performance studies, configures it and adds it to the analysis manager.
  
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
  reader->SetDebug(-1);//10 for lots of messages
  reader->SwitchOnEMCALCells();
  reader->SwitchOnPHOSCells();
  //Min particle pT
  reader->SetEMCALPtMin(0.2); 
  reader->SetPHOSPtMin(0.2);
  reader->SetCTSPtMin(0.2);
  //reader->Print("");
  
  // ##### Analysis algorithm settings ####

  AliFidutialCut * fidCut = new AliFidutialCut();
  fidCut->DoCTSFidutialCut(kFALSE) ;
  fidCut->DoEMCALFidutialCut(kTRUE) ;
  fidCut->DoPHOSFidutialCut(kTRUE) ;
		
  AliAnaCalorimeterQA *emcalQA = new AliAnaCalorimeterQA();
  emcalQA->SetDebug(-1); //10 for lots of messages
  emcalQA->SetCalorimeter("EMCAL");
  emcalQA->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  emcalQA->AddToHistogramsName("EMCAL"); //Begining of histograms name
  emcalQA->SetFidutialCut(fidCut);
  emcalQA->SwitchOnFidutialCut();
  //emcalQA->Print("");	
	
  AliAnaCalorimeterQA *phosQA = new AliAnaCalorimeterQA();
  phosQA->SetDebug(-1); //10 for lots of messages
  phosQA->SetCalorimeter("PHOS");
  phosQA->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  phosQA->AddToHistogramsName("PHOS");//Begining of histograms name
  phosQA->SetFidutialCut(fidCut);
  phosQA->SwitchOnFidutialCut();
  //phosQA->Print("");	
  
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(emcalQA,0);
  maker->AddAnalysis(phosQA,1);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  maker->Print("");
  
  printf("======================== \n");
  printf(" End Configuration of Calorimeter QA \n");
  printf("======================== \n");
  
   // Create task
   //===========================================================================
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation ("CalorimeterPerformance");
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  task->SetAnalysisMaker(maker);				
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer("Calo.Performance",TList::Class(),
							   AliAnalysisManager::kOutputContainer, "Calo.Performance.root");
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  
  return task;
}


