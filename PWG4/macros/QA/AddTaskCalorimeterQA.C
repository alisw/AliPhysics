AliAnalysisTaskParticleCorrelation *AddTaskCalorimeterQA(TString data, Bool_t kUseKinematics = kFALSE, Bool_t kPrintSettings = kFALSE)
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
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
   // Configure analysis
   //===========================================================================
  
  //Reader
  //For this particular analysis few things done by the reader.
  //Nothing else needs to be set.
  AliCaloTrackReader * reader = 0x0;
  if(data=="AOD")      reader = new AliCaloTrackAODReader();
  else if(data=="ESD") reader = new AliCaloTrackESDReader();
  //reader->SetDebug(10);//10 for lots of messages
   if(kPrintSettings) reader->Print("");
  
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
	
   reader->SetDeltaAODFileName(""); //Do not create deltaAOD file, this analysis do not create branches.

	
  // ##### Analysis algorithm settings ####
   //Only needed now for MC data
   //AliFiducialCut * fidCut = new AliFiducialCut();
   //fidCut->DoCTSFiducialCut(kFALSE) ;
   //fidCut->DoEMCALFiducialCut(kTRUE) ;
   //fidCut->DoPHOSFiducialCut(kTRUE) ;
		
  AliAnaCalorimeterQA *emcalQA = new AliAnaCalorimeterQA();
  //emcalQA->SetDebug(2); //10 for lots of messages
  emcalQA->SetCalorimeter("EMCAL");
  if(kUseKinematics) emcalQA->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else  emcalQA->SwitchOffDataMC() ;
  emcalQA->AddToHistogramsName("EMCAL_"); //Begining of histograms name
  //emcalQA->SetFiducialCut(fidCut);
  emcalQA->SwitchOffFiducialCut();
  if(kPrintSettings) emcalQA->Print("");	
  emcalQA->SwitchOnPlotsMaking();
  emcalQA->SwitchOnCalorimetersCorrelation();
  //Set Histrograms bins and ranges
  emcalQA->SetHistoPtRangeAndNBins(0, 10, 100) ;
  emcalQA->SetHistoPhiRangeAndNBins(75*TMath::DegToRad(), 125*TMath::DegToRad(), 100) ;
  emcalQA->SetHistoEtaRangeAndNBins(-0.8, 0.8, 160) ;
  emcalQA->SetNumberOfModules(4); //EMCAL first year
  //emcalQA->GetMCAnalysisUtils()->SetDebug(10);
	
   AliAnaCalorimeterQA *phosQA = new AliAnaCalorimeterQA();
   //phosQA->SetDebug(2); //10 for lots of messages
   phosQA->SetCalorimeter("PHOS");
   if(kUseKinematics) phosQA->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
   else  phosQA->SwitchOffDataMC() ;  
   phosQA->AddToHistogramsName("PHOS_");//Begining of histograms name
   //phosQA->SetFiducialCut(fidCut);
   phosQA->SwitchOffFiducialCut();
   if(kPrintSettings)phosQA->Print("");	
   //phosQA->GetMCAnalysisUtils()->SetDebug(10);
   phosQA->SwitchOnPlotsMaking();
    //Set Histrograms bins and ranges
  phosQA->SetHistoPtRangeAndNBins(0, 10, 100) ;
  phosQA->SetHistoPhiRangeAndNBins(215*TMath::DegToRad(), 325*TMath::DegToRad(), 100) ;
  phosQA->SetHistoEtaRangeAndNBins(-0.13, 0.13, 100) ;
  phosQA->SetNumberOfModules(3); //PHOS first year


  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(emcalQA,0);
  maker->AddAnalysis(phosQA,1);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOffAODsMaker()  ;
  if(kPrintSettings) maker->Print("");
 

  printf("======================== \n");
  printf(" End Configuration of Calorimeter QA \n");
  printf("======================== \n");
  
   // Create task
   //===========================================================================
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation ("CalorimeterPerformance");
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(0);
  task->SetAnalysisMaker(maker);				
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer("Calo.Performance",TList::Class(),
							   AliAnalysisManager::kOutputContainer, "Calo.Performance.root");
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  //mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  
  return task;
}


