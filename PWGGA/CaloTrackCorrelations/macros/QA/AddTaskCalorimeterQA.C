//
// Wagon contacts: EMCAL Gustavo.Conesa.Balbastre@cern.ch
//                
//
AliAnalysisTaskCaloTrackCorrelation *AddTaskCalorimeterQA(TString data, 
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
  
  if(kUseKinematics)
  {
    if(inputDataType == "ESD")
    {
      reader->SwitchOnStack();          
      reader->SwitchOffAODMCParticles(); 
    }
    else if(inputDataType == "AOD")
    {
      reader->SwitchOffStack();          
      reader->SwitchOnAODMCParticles(); 
    }
  }
  
  reader->SetDeltaAODFileName(""); //Do not create deltaAOD file, this analysis do not create branches.
  reader->SwitchOffWriteDeltaAOD()  ;
  reader->SetImportGeometryFromFile(kFALSE);
  
  if(kPrintSettings) reader->Print("");
  
  // *** Calorimeters Utils	***
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);

  AliEMCALRecoUtils* reco = cu->GetEMCALRecoUtils();
  reco->SwitchOnRejectExoticCell() ; // reject exotic cells, fill different histograms for exotic clusters and good clusters
  reco->SetExoticCellDiffTimeCut(10000); // Open  
  reco->SetExoticCellFractionCut(0.95);  // 1-Ecross/Ecell > 0.95 -> out
  reco->SetExoticCellMinAmplitudeCut(0.75); // 750 MeV    
  
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
  emcalQA->SwitchOffFillAllPositionHistogram2();
  emcalQA->SwitchOnStudyBadClusters();

  //Set Histrograms bins and ranges
  emcalQA->GetHistogramRanges()->SetHistoPtRangeAndNBins(0, 100, 200) ;
  emcalQA->GetHistogramRanges()->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  emcalQA->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-0.71, 0.71, 200) ;
  
  if     (year==2010)
  {  
    emcalQA->SetNumberOfModules(4); 
    emcalQA->GetHistogramRanges()->SetHistoPhiRangeAndNBins(79*TMath::DegToRad(), 121*TMath::DegToRad(), 200) ;
    emcalQA->GetHistogramRanges()->SetHistoXRangeAndNBins(-230,90,120);
    emcalQA->GetHistogramRanges()->SetHistoYRangeAndNBins(370,450,40);
  }
  else if(year==2011)
  {            
    emcalQA->SetNumberOfModules(10); 
    emcalQA->GetHistogramRanges()->SetHistoPhiRangeAndNBins(79*TMath::DegToRad(), 191*TMath::DegToRad(), 200) ;
    emcalQA->GetHistogramRanges()->SetHistoXRangeAndNBins(-600,90,200);
    emcalQA->GetHistogramRanges()->SetHistoYRangeAndNBins(100,450,100);
  }
  else 
  {
    emcalQA->SetNumberOfModules(12); 
    emcalQA->GetHistogramRanges()->SetHistoPhiRangeAndNBins(79*TMath::DegToRad(), 181*TMath::DegToRad(), 200) ; // revise
    emcalQA->GetHistogramRanges()->SetHistoXRangeAndNBins(-700,90,200); // revise
    emcalQA->GetHistogramRanges()->SetHistoYRangeAndNBins(50,450,100);  // revise     
  }
  
  emcalQA->GetHistogramRanges()->SetHistoMassRangeAndNBins(0., 1, 400) ;
  emcalQA->GetHistogramRanges()->SetHistoAsymmetryRangeAndNBins(0., 1. , 10 );
  emcalQA->GetHistogramRanges()->SetHistoPOverERangeAndNBins(0,3.,90);
  emcalQA->GetHistogramRanges()->SetHistodEdxRangeAndNBins(0.,200.,200);
  emcalQA->GetHistogramRanges()->SetHistodRRangeAndNBins(0.,0.15,150);
  emcalQA->GetHistogramRanges()->SetHistoTimeRangeAndNBins(300.,900,300);
  emcalQA->GetHistogramRanges()->SetHistoRatioRangeAndNBins(0.,2.,100);
  emcalQA->GetHistogramRanges()->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  emcalQA->GetHistogramRanges()->SetHistoNClusterCellRangeAndNBins(0,500,500);
  emcalQA->GetHistogramRanges()->SetHistoZRangeAndNBins(-400,400,200);
  emcalQA->GetHistogramRanges()->SetHistoRRangeAndNBins(400,450,25);
  emcalQA->GetHistogramRanges()->SetHistoV0SignalRangeAndNBins(0,5000,500);
  emcalQA->GetHistogramRanges()->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  emcalQA->GetHistogramRanges()->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
	
  if(kPrintSettings) emcalQA->Print("");	
  
  // #### Configure Maker ####
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
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
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("CalorimeterPerformance_%s",suffix));
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


