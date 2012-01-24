//
// Wagon contacts: EMCAL Gustavo.Conesa.Balbastre@cern.ch
//                
//
AliAnalysisTaskParticleCorrelation *AddTaskCalorimeterQA(TString data, 
                                                         Int_t year = 2011, 
                                                         Bool_t kPrintSettings = kFALSE,
                                                         Bool_t kSimulation = kFALSE,
                                                         TString outputFile = "", 
                                                         const char *suffix="default")
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
  
  Bool_t kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
  
  // Configure analysis
  //===========================================================================
  
  //Reader
  //For this particular analysis few things done by the reader.
  //Nothing else needs to be set.
  
  AliCaloTrackReader * reader = 0x0;
  if     (data.Contains("AOD")) reader = new AliCaloTrackAODReader();
  else if(data.Contains("ESD")) reader = new AliCaloTrackESDReader();
  //reader->SetDebug(10);//10 for lots of messages
  
  reader->SwitchOnEMCALCells(); 
  reader->SwitchOnEMCAL();
  reader->SwitchOnPHOSCells(); // For correlation plots
  reader->SwitchOnPHOS();      // For correlation plots
  reader->SetEMCALPtMin(0.); 
  reader->SwitchOnCTS();
  reader->SetCTSPtMin  (0.);
  reader->SetZvertexCut(10.);
  
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
  reader->SwitchOffWriteDeltaAOD()  ;
  if(kPrintSettings) reader->Print("");
  
  // *** Calorimeters Utils	***
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  if(year==2010) cu->SetEMCALGeometryName("EMCAL_FIRSTYEARV1");
  else           cu->SetEMCALGeometryName("EMCAL_COMPLETEV1"); 

  
  cu->SetDebug(-1);
  if(kPrintSettings) cu->Print("");	
  
  // ##### Analysis algorithm settings ####
  
  AliAnaCalorimeterQA *emcalQA = new AliAnaCalorimeterQA();
  //emcalQA->SetDebug(10); //10 for lots of messages
  emcalQA->SetCalorimeter("EMCAL");
  if(kUseKinematics) emcalQA->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else               emcalQA->SwitchOffDataMC() ;
  emcalQA->AddToHistogramsName("EMCAL_"); //Begining of histograms name
  emcalQA->SwitchOffFiducialCut();
  emcalQA->SwitchOnCorrelation();
  emcalQA->SwitchOffFillAllTH3Histogram();
  emcalQA->SwitchOffFillAllPositionHistogram();

  //Set Histrograms bins and ranges
  emcalQA->SetHistoPtRangeAndNBins(0, 100, 200) ;
  emcalQA->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  emcalQA->SetHistoEtaRangeAndNBins(-0.71, 0.71, 200) ;
  
  if(year==2010){  
    emcalQA->SetNumberOfModules(4); 
    emcalQA->SetHistoPhiRangeAndNBins(79*TMath::DegToRad(), 121*TMath::DegToRad(), 200) ;
    emcalQA->SetHistoXRangeAndNBins(-230,90,120);
    emcalQA->SetHistoYRangeAndNBins(370,450,40);
  }
  else{            
    emcalQA->SetNumberOfModules(10); 
    emcalQA->SetHistoPhiRangeAndNBins(79*TMath::DegToRad(), 181*TMath::DegToRad(), 200) ;
    emcalQA->SetHistoXRangeAndNBins(-600,90,200);
    emcalQA->SetHistoYRangeAndNBins(100,450,100);
  }
  
  emcalQA->SetHistoMassRangeAndNBins(0., 1, 400) ;
  emcalQA->SetHistoAsymmetryRangeAndNBins(0., 1. , 10 );
  emcalQA->SetHistoPOverERangeAndNBins(0,10.,100);
  emcalQA->SetHistodEdxRangeAndNBins(0.,200.,200);
  emcalQA->SetHistodRRangeAndNBins(0.,TMath::Pi(),150);
  emcalQA->SetHistoTimeRangeAndNBins(300.,900,300);
  emcalQA->SetHistoRatioRangeAndNBins(0.,2.,100);
  emcalQA->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  emcalQA->SetHistoNClusterCellRangeAndNBins(0,500,500);
  emcalQA->SetHistoZRangeAndNBins(-400,400,200);
  emcalQA->SetHistoRRangeAndNBins(400,450,25);
  emcalQA->SetHistoV0SignalRangeAndNBins(0,5000,500);
  emcalQA->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  emcalQA->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
	
  if(kPrintSettings) emcalQA->Print("");	
  
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->SetCaloUtils(cu); //pointer to calorimeter utils
  maker->AddAnalysis(emcalQA,0);
  maker->SetAnaDebug(-1)  ; // 0 to at least print the event number
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOffAODsMaker()  ;
  if(kPrintSettings) maker->Print("");
  
  printf("======================== \n");
  printf(" End Configuration of Calorimeter QA \n");
  printf("======================== \n");
  
  // Create task
  //===========================================================================
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation (Form("CalorimeterPerformance_%s",suffix));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  task->SetAnalysisMaker(maker);	
  task->SetBranches("ESD:AliESDRun.,AliESDHeader"); //just a trick to get Constantin's analysis to work
  mgr->AddTask(task);
  
  //Create containers
  //  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer("Calo.Performance",TList::Class(),
  //							   AliAnalysisManager::kOutputContainer, "Calo.Performance.root");
	
  TString cname;
  if(outputFile.Length()==0)outputFile = AliAnalysisManager::GetCommonFileName(); 
  
  
  cname = Form("CaloQA_%s", suffix);
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(cname, TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s:%s",outputFile.Data(),cname.Data()));
  
  cname = Form("CaloQACuts_%s", suffix);
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(cname, TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             Form("%s:%s",outputFile.Data(),cname.Data()));
	//Form("%s:PartCorrCuts",outputfile.Data()));	
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  
  return task;
}


