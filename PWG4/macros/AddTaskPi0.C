AliAnalysisTaskParticleCorrelation *AddTaskPi0(TString data, TString calorimeter, Bool_t kPrintSettings = kFALSE,Bool_t kSimulation = kFALSE, Bool_t outputAOD=kFALSE, TString outputfile = "", Int_t year = 2010,TString col = "pp",Bool_t oldAOD=kFALSE)
{
  // Creates a PartCorr task, configures it and adds it to the analysis manager.
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPi0", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPi0", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = "AOD";
  if(!data.Contains("delta"))
    inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  //cout<<"DATA TYPE :: "<<inputDataType<<endl;
  // inputDataType: data managed by the input handler
  // data: can be same as one managed by input handler, or the output AOD created by the filter. By default use AOD
  
  Bool_t kUseKinematics = kFALSE; 
  if(kSimulation) { 
    kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE; 
    if (!kUseKinematics && data=="AOD" && inputDataType != "ESD") kUseKinematics = kTRUE; //AOD primary should be available ... 
  } 
  
  cout<<"********* ACCESS KINE? "<<kUseKinematics<<endl;
  
  // Configure analysis
  //===========================================================================
  
  // *** Reader ***
  AliCaloTrackReader * reader = ;
  if(data.Contains("AOD")) reader = new AliCaloTrackAODReader();
  else if(data=="ESD") reader = new AliCaloTrackESDReader();
  else if(data=="MC" && inputDataType == "ESD") reader = new AliCaloTrackMCReader();
  reader->SetDebug(-1);//10 for lots of messages
  reader->SwitchOnCTS();
  //reader->SetDeltaAODFileName("");
  //if(!kSimulation) reader->SetFiredTriggerClassName("CINT1B-ABCE-NOPF-ALL");
  if(calorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();  
    reader->SwitchOnEMCAL();
  }
  if(calorimeter == "PHOS") { 
    reader->SwitchOnPHOSCells();  
    reader->SwitchOnPHOS();
  }

   reader->SwitchOffSuspiciousClustersRemoval();  //EMCAL
  
  // for case data="deltaAOD", no need to fill the EMCAL/PHOS cluster lists
  if(data.Contains("delta")){
    reader->SwitchOffEMCAL();
    reader->SwitchOffPHOS();
    reader->SwitchOffEMCALCells(); 
    reader->SwitchOffPHOSCells(); 
  }
  
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

  //In case of AODs created only for calorimeters, and track information filtered
  //CTS is off when calling this method
  //reader->SwitchOnCaloFilterPatch();
  //Event selection
  if     (col=="pp"  ) {
    reader->SwitchOnEventSelection(); //remove pileup by default
    reader->SwitchOffV0ANDSelection() ; // and besides v0 AND
    reader->SwitchOffPrimaryVertexSelection(); // and besides primary vertex
  }
  else if(col=="PbPb") {
    reader->SwitchOffEventSelection(); //remove pileup by default
    reader->SwitchOffV0ANDSelection() ; // and besides v0 AND
    reader->SwitchOffPrimaryVertexSelection(); // and besides primary vertex
  }
  
  if     (col=="pp"  )   reader->SetZvertexCut(50.);  //Open cut
  else if(col=="PbPb")   reader->SetZvertexCut(10.);  //Centrality defined in this range.
  
  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.);
  reader->SetCTSPtMin(0.);
  if(outputAOD)  reader->SwitchOnWriteDeltaAOD()  ;
  if(oldAOD) reader->SwitchOnOldAODs();
  if(kPrintSettings) reader->Print("");
  
  // *** Calorimeters Utils	***
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  if(year==2010) cu->SetEMCALGeometryName("EMCAL_FIRSTYEARV1");
  else           cu->SetEMCALGeometryName("EMCAL_COMPLETEV1");
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  // Remove EMCAL hottest channels for first LHC10 periods 	
  //  cu->SwitchOnBadChannelsRemoval();
  //  cu->SwitchOnDistToBadChannelRecalculation();

//   TFile * fbad = new TFile("BadChannels.root","read");
//   TH2I * hbad0 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod0");
//   TH2I * hbad1 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod1");
//   TH2I * hbad2 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod2");
//   TH2I * hbad3 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod3");
//   cu->SetEMCALChannelStatusMap(0,hbad0);
//   cu->SetEMCALChannelStatusMap(1,hbad1);
//   cu->SetEMCALChannelStatusMap(2,hbad2);
//   cu->SetEMCALChannelStatusMap(3,hbad3);
  
  
  //Recalibration
  //cu->SwitchOnRecalibration();
  //TFile * f = new TFile("RecalibrationFactors.root","read");
  //cu->SetEMCALChannelRecalibrationFactors(0,(TH2F*)f->Get("EMCALRecalFactors_SM0"));
  //cu->SetEMCALChannelRecalibrationFactors(1,(TH2F*)f->Get("EMCALRecalFactors_SM1"));
  //cu->SetEMCALChannelRecalibrationFactors(2,(TH2F*)f->Get("EMCALRecalFactors_SM2"));
  //cu->SetEMCALChannelRecalibrationFactors(3,(TH2F*)f->Get("EMCALRecalFactors_SM3"));
  //f->Close();	
  
  cu->SetDebug(-1);
  if(kPrintSettings) cu->Print("");
  
  
  // ##### Analysis algorithm settings ####
  
  // -------------------------------------------------
  // --- Photon/Pi0/Omega/Electron Analysis ---
  // -------------------------------------------------
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  if(calorimeter == "PHOS"){
    anaphoton->SetNCellCut(2);// At least 2 cells
    anaphoton->SetMinPt(3.);
    anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
  }
  else {//EMCAL
    anaphoton->SetNCellCut(1);// At least 2 cells
    anaphoton->SetMinPt(0.5); // no effect minium EMCAL cut.
    //if(!kUseKinematics) anaphoton->SetTimeCut(400,900);// Time window of [400-900] ns
    //anaphoton->SetMinDistanceToBadChannel(6, 12, 18);//For officially produced ESDs/AODs
    anaphoton->SetMinDistanceToBadChannel(1, 2, 3);//For filtered AODs, new releases.
  }
  anaphoton->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaphoton->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaphoton->SwitchOffDataMC() ;
  anaphoton->SwitchOffCaloPID();
  anaphoton->SwitchOffFiducialCut();
  if(kSimulation){
    anaphoton->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anaphoton->GetFiducialCut();
    fidCut1stYear->DoCTSFiducialCut(kFALSE) ;
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
  }
  
  if(!data.Contains("delta")) {
    anaphoton->SetOutputAODName(Form("Photons%s",calorimeter.Data()));
    anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else anaphoton->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anaphoton->AddToHistogramsName("AnaPhotonCorr_");
  //Set Histograms bins and ranges
  anaphoton->SetHistoPtRangeAndNBins(0, 30, 150) ;
  if(year==2010)anaphoton->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
  else          anaphoton->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
  anaphoton->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  if(kPrintSettings) anaphoton->Print("");
  
  // -----------------------------------
  // --- Pi0 Invariant Mass Analysis ---
  // -----------------------------------
  
  AliAnaPi0 *anapi0 = new AliAnaPi0();
  anapi0->SetDebug(-1);//10 for lots of messages
  anapi0->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anapi0->SetCalorimeter(calorimeter);
  anapi0->SwitchOffMultipleCutAnalysis(); 
  //anapi0->SetNPtCuts(2);
  //anapi0->SetNAsymCuts(2);
  //anapi0->SetNNCellCuts(2);
  //anapi0->SetNPIDBits(1);
  
  //anapi0->SetPtCutsAt(0,0.3); anapi0->SetPtCutsAt(1,0.5);
  //anapi0->SetAsymCutsAt(0,0.1);anapi0->SetAsymCutsAt(1,0.5);
  //anapi0->SetNCellCutsAt(0,1); anapi0->SetNCellCutsAt(1,2);
  //anapi0->SetPIDBitsAt(0,0); //No Cut
  //anapi0->SetPIDBitsAt(1,2); //Dispersion Cut

  
  if(kSimulation){
    anapi0->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anapi0->GetFiducialCut();
    fidCut1stYear->DoCTSFiducialCut(kFALSE) ;
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
  }  
	
  //settings for pp collision mixing
  anapi0->SwitchOnOwnMix(); //Off when mixing done with general mixing frame
  if     (col=="pp"  ) {
    anapi0->SetNCentrBin(1);
    anapi0->SetNZvertBin(50);
    anapi0->SwitchOnSMCombinations();
  }
  else if(col=="PbPb") {
    anapi0->SetNCentrBin(10);
    anapi0->SetNZvertBin(10);
    anapi0->SwitchOffSMCombinations();
  }
  anapi0->SetNRPBin(1);
  anapi0->SetNMaxEvMix(50);
  
  if(kUseKinematics)anapi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else              anapi0->SwitchOffDataMC() ;
  if(calorimeter=="PHOS") anapi0->SetNumberOfModules(3); //PHOS first year
  else {                   
     if(year==2010) anapi0->SetNumberOfModules(4); //EMCAL first year
    else            anapi0->SetNumberOfModules(10);
  }
  anapi0->SetHistoPtRangeAndNBins(0, 30, 150) ;    
  anapi0->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  if(year==2010)anaphoton->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
  else          anaphoton->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;

  anapi0->SetHistoMassRangeAndNBins(0., 1., 200) ;
  anapi0->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  anapi0->SetHistoTrackMultiplicityRangeAndNBins(0, 200, 20); 

  if(kPrintSettings) anapi0->Print("");

  
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->SetCaloUtils(cu); //pointer to calorimeter utils
  Int_t n = 0;//Analysis number, order is important
  // Particle selection analysis
  maker->AddAnalysis(anaphoton,n++);
  maker->AddAnalysis(anapi0,n++);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(data.Contains("delta")) maker->SwitchOffAODsMaker()  ;
  else                       maker->SwitchOnAODsMaker()  ;
	
  if(kPrintSettings) maker->Print("");
  
  printf("======================== \n");
  printf(" End Configuration of Pi0 analysis with detector %s \n",calorimeter.Data());
  printf("======================== \n");
  
  // Create task
  //===========================================================================
  AliAnalysisTaskParticleCorrelation * task = new AliAnalysisTaskParticleCorrelation (Form("PartCorr%s",calorimeter.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  if(data=="ESD")task->SelectCollisionCandidates();
  task->SetAnalysisMaker(maker);
  //if(!kSimulation)task->SelectCollisionCandidates(); //AliPhysicsSelection has to be attached before.
  mgr->AddTask(task);
  
  //Create containers
  char name[128];
  sprintf(name,"PartCorr_%s",calorimeter.Data());
  cout<<"Name of task "<<name<<endl;
  
  if(outputfile.Length()==0)outputfile = AliAnalysisManager::GetCommonFileName(); 

  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(calorimeter.Data(), TList::Class(), 
                                                              AliAnalysisManager::kOutputContainer, 
                                                              Form("%s:Pi0",outputfile.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("%sCuts",calorimeter.Data()), TList::Class(), 
                                                              AliAnalysisManager::kParamContainer, 
                                                              Form("%s:Pi0Cuts",outputfile.Data()));
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  if(!data.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  
  return task;
}


