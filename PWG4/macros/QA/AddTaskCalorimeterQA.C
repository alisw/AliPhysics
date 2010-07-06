//
// Wagon contacts: EMCAL Gustavo.Conesa.Balbastre@cern.ch
//                 PHOS  Yuri.Kharlov@cern.ch
//

AliAnalysisTaskParticleCorrelation *AddTaskCalorimeterQA(TString data, Bool_t kPrintSettings = kFALSE,Bool_t kSimulation = kFALSE)
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
  AliCaloTrackReader * reader = 0x0;
  if(data=="AOD")      reader = new AliCaloTrackAODReader();
  else if(data=="ESD") reader = new AliCaloTrackESDReader();
  //reader->SetDebug(10);//10 for lots of messages
  reader->SwitchOffEMCAL();
  reader->SwitchOffPHOS();
  reader->SwitchOffCTS();
	
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
   if(kPrintSettings) reader->Print("");
	
  // *** Calorimeters Utils	***
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  cu->SwitchOnNoFiducialBorderInEMCALEta0();
	
  // Remove EMCAL hottest channels for first LHC10 periods 	
  cu->SwitchOnBadChannelsRemoval();
  // SM0
  cu->SetEMCALChannelStatus(0,3,13);  cu->SetEMCALChannelStatus(0,44,1); cu->SetEMCALChannelStatus(0,3,13); 
  cu->SetEMCALChannelStatus(0,20,7);  cu->SetEMCALChannelStatus(0,38,2);   
  // SM1
  cu->SetEMCALChannelStatus(1,4,7);   cu->SetEMCALChannelStatus(1,4,13);  cu->SetEMCALChannelStatus(1,9,20); 
  cu->SetEMCALChannelStatus(1,14,15); cu->SetEMCALChannelStatus(1,23,16); cu->SetEMCALChannelStatus(1,32,23); 
  cu->SetEMCALChannelStatus(1,37,5);  cu->SetEMCALChannelStatus(1,40,1);  cu->SetEMCALChannelStatus(1,40,2);
  cu->SetEMCALChannelStatus(1,40,5);  cu->SetEMCALChannelStatus(1,41,0);  cu->SetEMCALChannelStatus(1,41,1);
  cu->SetEMCALChannelStatus(1,41,2);  cu->SetEMCALChannelStatus(1,41,4);
  // SM2 	
  cu->SetEMCALChannelStatus(2,14,15); cu->SetEMCALChannelStatus(2,18,16); cu->SetEMCALChannelStatus(2,18,17); 
  cu->SetEMCALChannelStatus(2,18,18); cu->SetEMCALChannelStatus(2,18,20); cu->SetEMCALChannelStatus(2,18,21); 
  cu->SetEMCALChannelStatus(2,18,23); cu->SetEMCALChannelStatus(2,19,16); cu->SetEMCALChannelStatus(2,19,17); 
  cu->SetEMCALChannelStatus(2,19,19); cu->SetEMCALChannelStatus(2,19,20); cu->SetEMCALChannelStatus(2,19,21); 
  cu->SetEMCALChannelStatus(2,19,22);
  //SM3
  cu->SetEMCALChannelStatus(3,4,7);
	
  //Recalibration
  //cu->SwitchOnRecalibration();
  //TFile * f = new TFile("RecalibrationFactors.root","read");
  //cu->SetEMCALChannelRecalibrationFactors(0,(TH2F*)f->Get("EMCALRecalFactors_SM0"));
  //cu->SetEMCALChannelRecalibrationFactors(1,(TH2F*)f->Get("EMCALRecalFactors_SM1"));
  //cu->SetEMCALChannelRecalibrationFactors(2,(TH2F*)f->Get("EMCALRecalFactors_SM2"));
  //cu->SetEMCALChannelRecalibrationFactors(3,(TH2F*)f->Get("EMCALRecalFactors_SM3"));
	
  cu->SetDebug(-1);
  if(kPrintSettings) cu->Print("");	
	
  // ##### Analysis algorithm settings ####
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
  emcalQA->SwitchOffPlotsMaking();
  emcalQA->SwitchOnCalorimetersCorrelation();
  emcalQA->SetTimeCut(420,825);//Open for the moment
  //Set Histrograms bins and ranges
  emcalQA->SetHistoPtRangeAndNBins(0, 10, 100) ;
  emcalQA->SetHistoPhiRangeAndNBins(75*TMath::DegToRad(), 125*TMath::DegToRad(), 100) ;
  emcalQA->SetHistoEtaRangeAndNBins(-0.8, 0.8, 80) ;
  emcalQA->SetNumberOfModules(4); //EMCAL first year
  emcalQA->SetHistoMassRangeAndNBins(0., 0.6, 200) ;
  emcalQA->SetHistoAsymmetryRangeAndNBins(0., 1. , 10) ;
  emcalQA->SetHistoPOverERangeAndNBins(0,10.,100);
  emcalQA->SetHistodEdxRangeAndNBins(0.,400.,200);
  emcalQA->SetHistodRRangeAndNBins(0.,TMath::Pi(),150);
  emcalQA->SetHistoTimeRangeAndNBins(400.,800,200);
  emcalQA->SetHistoRatioRangeAndNBins(0.,2.,100);
  emcalQA->SetHistoVertexDistRangeAndNBins(0.,500.,100);
  emcalQA->SetHistoNClusterCellRangeAndNBins(0,300,300);
  emcalQA->SetHistoXRangeAndNBins(-230,90,160);
  emcalQA->SetHistoYRangeAndNBins(370,450,35);
  emcalQA->SetHistoZRangeAndNBins(-400,400,100);
  emcalQA->SetHistoRRangeAndNBins(400,450,25);

  //emcalQA->GetMCAnalysisUtils()->SetDebug(10);
	
  if(kPrintSettings) emcalQA->Print("");	

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
   //Set Histrograms bins and ranges
   phosQA->SetHistoPtRangeAndNBins(0, 10, 100) ;
   phosQA->SetHistoPhiRangeAndNBins(255*TMath::DegToRad(), 325*TMath::DegToRad(), 100) ;
   phosQA->SetHistoEtaRangeAndNBins(-0.13, 0.13, 80) ;
   phosQA->SetNumberOfModules(3); //PHOS first year
   phosQA->SetHistoMassRangeAndNBins(0., 0.6, 200) ;
   phosQA->SetHistoAsymmetryRangeAndNBins(0., 1. , 25) ;
   phosQA->SetHistoPOverERangeAndNBins(0,10.,100);
   phosQA->SetHistodEdxRangeAndNBins(0.,400.,200);
   phosQA->SetHistodRRangeAndNBins(0.,TMath::Pi(),150);
   phosQA->SetHistoTimeRangeAndNBins(0.,400,200);
   phosQA->SetHistoRatioRangeAndNBins(0.,2.,100);
   phosQA->SetHistoVertexDistRangeAndNBins(0.,500.,100);
   phosQA->SetHistoNClusterCellRangeAndNBins(0,300,300);
   phosQA->SetHistoXRangeAndNBins(-100,400,150);
   phosQA->SetHistoYRangeAndNBins(-490,-290,100);
   phosQA->SetHistoZRangeAndNBins(-80,80,80);
   phosQA->SetHistoRRangeAndNBins(440,480,40);
	
  //if(kPrintSettings)phosQA->Print("");	
	
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->SetCaloUtils(cu); //pointer to calorimeter utils
  maker->AddAnalysis(emcalQA,0);
  maker->AddAnalysis(phosQA,1);
  maker->SetAnaDebug(-1)  ; // 0 to at least print the event number
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
  //task->SetDebugLevel(-1);
  task->SelectCollisionCandidates();
  task->SetAnalysisMaker(maker);	
  mgr->AddTask(task);
  
  //Create containers
//  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer("Calo.Performance",TList::Class(),
//							   AliAnalysisManager::kOutputContainer, "Calo.Performance.root");
	
  TString outputfile = AliAnalysisManager::GetCommonFileName(); 
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer("CaloQA", TList::Class(), 
															 AliAnalysisManager::kOutputContainer, 
															 Form("%s:CaloQA",outputfile.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer("CaloQACuts", TList::Class(), 
															  AliAnalysisManager::kParamContainer, 
															  Form("%s:CaloQACuts",outputfile.Data()));
	//Form("%s:PartCorrCuts",outputfile.Data()));	
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);

  return task;
}


