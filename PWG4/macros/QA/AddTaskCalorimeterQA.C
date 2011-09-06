//
// Wagon contacts: EMCAL Gustavo.Conesa.Balbastre@cern.ch
//                 PHOS  Yuri.Kharlov@cern.ch
//
AliAnalysisTaskParticleCorrelation *AddTaskCalorimeterQA(TString data, Int_t year = 2011, Bool_t kPrintSettings = kFALSE,Bool_t kSimulation = kFALSE,TString outputFile = "", const char *suffix="default")
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
	
  cout<<"********* ACCESS KINE? "<<kUseKinematics<<endl;
	
  // Configure analysis
  //===========================================================================
  
  //Reader
  //For this particular analysis few things done by the reader.
  //Nothing else needs to be set.
  Bool_t doEmcal = 1;
  Bool_t doPhos  = 1;
  if (data.Contains("PHOS"))
    doEmcal = 0;
  else if (data.Contains("EMCAL"))
    doPhos = 0;

  AliCaloTrackReader * reader = 0x0;
  if (data.Contains("AOD")) reader = new AliCaloTrackAODReader();
  else if(data.Contains("ESD")) reader = new AliCaloTrackESDReader();
  //reader->SetDebug(10);//10 for lots of messages
  if (doEmcal) {
    reader->SwitchOnEMCALCells(); 
    reader->SwitchOnEMCAL();
    reader->SetEMCALPtMin(0.); 
  }
  if (doPhos) {
    reader->SwitchOnPHOSCells(); 
    reader->SwitchOnPHOS();
    reader->SetPHOSPtMin (0.);
  }
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
  //if(!kSimulation) reader->SetFiredTriggerClassName("CINT1B-ABCE-NOPF-ALL");
  reader->SetDeltaAODFileName(""); //Do not create deltaAOD file, this analysis do not create branches.
  reader->SwitchOffWriteDeltaAOD()  ;
  if(kPrintSettings) reader->Print("");
  
  // *** Calorimeters Utils	***
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  if(year==2010){
    cu->SetEMCALGeometryName("EMCAL_FIRSTYEARV1");
    //cu->SwitchOnBadChannelsRemoval(); // just a trick
  }
  else {
    cu->SetEMCALGeometryName("EMCAL_COMPLETEV1"); 
    // Remove EMCAL hottest channels for LHC11 periods 	 
    cu->SwitchOnBadChannelsRemoval();
    Int_t badAbsID[]={74, 103, 152, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 368, 369, 370, 371, 372, 373, 374,375, 376, 377, 378, 379, 380, 381, 382, 383, 917, 1275, 1288, 1519, 1595, 1860, 1967, 2022, 2026, 2047, 2117, 2298, 2540, 2776, 3135, 3764, 6095, 6111, 6592, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 7371, 7425, 7430, 7457, 7491, 7709, 8352, 8353, 8356, 8357, 8808, 8810, 8812, 8814, 9056, 9769, 9815, 9837};
    
    Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
    AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
    for(Int_t i=0;i < 89; i++){
      geom->GetCellIndex(badAbsID[i],iSM,iMod,iIphi,iIeta); 
      // Gives SuperModule and Tower numbers
      geom->GetCellPhiEtaIndexInSModule(iSM,iMod,
                                        iIphi, iIeta,iRow,iCol);
      //printf("bad ID %d, col %d, row %d, sm %d\n",badAbsID[i],iCol,iRow,iSM);
      cu->SetEMCALChannelStatus(iSM , iCol, iRow,1);
    }
    
  }
  
  cu->SetDebug(-1);
  if(kPrintSettings) cu->Print("");	
  
  // ##### Analysis algorithm settings ####
  
  AliAnaCalorimeterQA *emcalQA = new AliAnaCalorimeterQA();
  //emcalQA->SetDebug(10); //10 for lots of messages
  emcalQA->SetCalorimeter("EMCAL");
  if(kUseKinematics) emcalQA->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else  emcalQA->SwitchOffDataMC() ;
  emcalQA->AddToHistogramsName("EMCAL_"); //Begining of histograms name
  //emcalQA->SetFiducialCut(fidCut);
  emcalQA->SwitchOffFiducialCut();
  emcalQA->SwitchOffPlotsMaking();
  emcalQA->SwitchOnCorrelation();
  emcalQA->SwitchOffFillAllTH3Histogram();
  //emcalQA->SwitchOffFillAllPositionHistogram2();
//  if(!kUseKinematics)emcalQA->SetTimeCut(400,850);//Open for the moment
  //Set Histrograms bins and ranges
  emcalQA->SetHistoPtRangeAndNBins(0, 50, 200) ;
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
  //emcalQA->GetMCAnalysisUtils()->SetDebug(10);
  if (!doPhos)
    emcalQA->SwitchOffCorrelation();
	
  if(kPrintSettings&&doEmcal) emcalQA->Print("");	
  
  AliAnaCalorimeterQA *phosQA = new AliAnaCalorimeterQA();
  //phosQA->SetDebug(2); //10 for lots of messages
  phosQA->SetCalorimeter("PHOS");
  if(kUseKinematics) phosQA->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else  phosQA->SwitchOffDataMC() ;  
  phosQA->AddToHistogramsName("PHOS_");//Begining of histograms name
  //phosQA->SetFiducialCut(fidCut);
  phosQA->SwitchOffFiducialCut();
  //phosQA->GetMCAnalysisUtils()->SetDebug(10);
  phosQA->SwitchOffPlotsMaking();
  phosQA->SwitchOffFillAllTH3Histogram();
  phosQA->SwitchOffFillAllPositionHistogram2();
  //Set Histrograms bins and ranges
  phosQA->SetHistoPtRangeAndNBins(0, 50, 200) ;
  phosQA->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  phosQA->SetHistoPhiRangeAndNBins(259*TMath::DegToRad(), 321*TMath::DegToRad(), 130) ;
  phosQA->SetHistoEtaRangeAndNBins(-0.125, 0.125, 57) ;
  phosQA->SetNumberOfModules(3); //PHOS first year
  phosQA->SetHistoMassRangeAndNBins(0., 1., 400) ;
  phosQA->SetHistoAsymmetryRangeAndNBins(0., 1. , 10) ;
  phosQA->SetHistoPOverERangeAndNBins(0,10.,100);
  phosQA->SetHistodEdxRangeAndNBins(0.,200.,200);
  phosQA->SetHistodRRangeAndNBins(0.,TMath::Pi(),150);
  phosQA->SetHistoTimeRangeAndNBins(0.,300,300);
  phosQA->SetHistoRatioRangeAndNBins(0.,2.,100);
  phosQA->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  phosQA->SetHistoNClusterCellRangeAndNBins(0,500,500);
  phosQA->SetHistoXRangeAndNBins(-100,400,100);
  phosQA->SetHistoYRangeAndNBins(-490,-290,100);
  phosQA->SetHistoZRangeAndNBins(-80,80,100);
  phosQA->SetHistoRRangeAndNBins(440,480,80);  
  phosQA->SetHistoV0SignalRangeAndNBins(0,5000,500);
  phosQA->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  phosQA->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
  if (!doEmcal)
    phosQA->SwitchOffCorrelation();
  if(kPrintSettings&&doPhos) phosQA->Print("");	

  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->SetCaloUtils(cu); //pointer to calorimeter utils
  Int_t n=0;//analysis counter
  if (doEmcal)
    maker->AddAnalysis(emcalQA,n++);
  if (doPhos)
    maker->AddAnalysis(phosQA,n);
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


