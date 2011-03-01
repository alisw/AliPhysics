AliAnalysisTaskParticleCorrelation *AddTaskPartCorr(TString inputDataType, TString calorimeter, Bool_t kPrintSettings = kFALSE,Bool_t kSimulation = kFALSE, Bool_t outputAOD=kFALSE, Bool_t oldAOD=kFALSE)
{
  // Creates a PartCorr task, configures it and adds it to the analysis manager.
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPartCorr", "No analysis manager to connect to.");
    return NULL;
  }  
 
  Bool_t kUseKinematics = kFALSE; 
  if(kSimulation) { 
    kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE; 
    if (!kUseKinematics && inputDataType == "AOD") kUseKinematics = kTRUE; //AOD primary should be available ... 
  } 
  
  cout<<"********* ACCESS KINE? "<<kUseKinematics<<endl;
  
  // Configure analysis
  //===========================================================================
  
  // *** Reader ***
  AliCaloTrackReader * reader =0x0 ;
  if(inputDataType.Contains("AOD")) reader = new AliCaloTrackAODReader();
  else if(inputDataType=="ESD") reader = new AliCaloTrackESDReader();
  else if(inputDataType=="MC" && inputDataType == "ESD") reader = new AliCaloTrackMCReader();
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
  
   reader->SwitchOnSuspiciousClustersRemoval();  //EMCAL

  // for case inputDataType="deltaAOD", no need to fill the EMCAL/PHOS cluster lists
  if(inputDataType.Contains("delta")){
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
  
  reader->SetZvertexCut(10.);
  
  //Min particle pT
  reader->SetEMCALPtMin(0.3); 
  reader->SetPHOSPtMin(0.3);
  reader->SetCTSPtMin(0.1);
  if(outputAOD)  reader->SwitchOnWriteDeltaAOD()  ;
  if(oldAOD) reader->SwitchOnOldAODs();
  if(kPrintSettings) reader->Print("");
  
  // *** Calorimeters Utils	***
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
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
  //f->Close();	
  
  cu->SetDebug(-1);
  if(kPrintSettings) cu->Print("");
  
  
  // ##### Analysis algorithm settings ####
  
  // -------------------------------------------------
  // --- Photon/Pi0/Omega/Electron Analysis ---
  // -------------------------------------------------
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  //settings for different multiplicity analysis
  anaphoton->SwitchOffEventSelection() ;
  anaphoton->SetMultiplicity(80, 120);

  if(calorimeter == "PHOS"){
    anaphoton->SetNCellCut(2);// At least 3 cells
    anaphoton->SetMinPt(0.3);
    anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
  }
  else {//EMCAL
    anaphoton->SetNCellCut(1);// At least 2 cells
    anaphoton->SetMinPt(0.3); 
    //if(!kUseKinematics) anaphoton->SetTimeCut(400,900);// Time window of [400-900] ns
    //anaphoton->SetMinDistanceToBadChannel(6, 12, 18);
    anaphoton->SetMinDistanceToBadChannel(1, 2, 3);//For new releases.
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
  
  if(!inputDataType.Contains("delta")) {
    anaphoton->SetOutputAODName(Form("Photons%s",calorimeter.Data()));
    anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else anaphoton->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anaphoton->AddToHistogramsName("AnaPhotonCorr_");
  //Set Histograms bins and ranges
  anaphoton->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anaphoton->Print("");
  
  // -----------------------------------
  // --- Pi0 Invariant Mass Analysis ---
  // -----------------------------------
  
  AliAnaPi0 *anapi0 = new AliAnaPi0();
  anapi0->SetDebug(-1);//10 for lots of messages
  anapi0->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anapi0->SetCalorimeter(calorimeter);
  anapi0->SwitchOnMultipleCutAnalysis(); 
  if(kSimulation){
    anapi0->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anapi0->GetFiducialCut();
    fidCut1stYear->DoCTSFiducialCut(kFALSE) ;
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
  }  	

  //settings for pp collision
  anapi0->SwitchOnOwnMix();
  anapi0->SwitchOnEventSelection() ;
  anapi0->SetNCentrBin(1);
  //anapi0->SetMultiplicity(80, 120);
  anapi0->SetMultiBin(1);  
  if(kUseKinematics)anapi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else anapi0->SwitchOffDataMC() ;
  if(calorimeter=="PHOS") anapi0->SetNumberOfModules(3); //PHOS first year
  else  anapi0->SetNumberOfModules(4); //EMCAL first year
  anapi0->SetHistoPtRangeAndNBins(0, 20, 200) ;
  //anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //anapi0->SetHistoEtaRangeAndNBins(-0.8, 0.8, 200) ;
  anapi0->SetHistoMassRangeAndNBins(0., 0.9, 300) ;
  anapi0->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  anapi0->SetHistoTrackMultiplicityRangeAndNBins(0, 200, 20); 

  if(kPrintSettings) anapi0->Print("");
	
  //---------------------------  
  //Pi0, event by event
  //---------------------------  
  
  AliAnaPi0EbE *anapi0ebe = new AliAnaPi0EbE();
  anapi0ebe->SwitchOffEventSelection() ;
  anapi0ebe->SetMultiplicity(80, 120);
  anapi0ebe->SetMultiBin(1);  
  anapi0ebe->SetDebug(-1);//10 for lots of messages
  anapi0ebe->SetAnalysisType(AliAnaPi0EbE::kIMCalo);
  anapi0ebe->SetMinPt(0);
  anapi0ebe->SetCalorimeter(calorimeter);
  anapi0ebe->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  if(!inputDataType.Contains("delta")) {
    anapi0ebe->SetOutputAODName(Form("Pi0s%s",calorimeter.Data()));
    anapi0ebe->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else  anapi0ebe->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  
  if(kUseKinematics) anapi0ebe->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anapi0ebe->SwitchOffDataMC() ;	
  
  AliNeutralMesonSelection *nms = anapi0ebe->GetNeutralMesonSelection();
  nms->SetInvMassCutRange(0.08, 0.18)     ;
  nms->KeepNeutralMesonSelectionHistos(kTRUE);
  //Set Histrograms bins and ranges
  if(calorimeter=="EMCAL" ){
    nms->SetHistoERangeAndNBins(0, 15, 150) ;  
    anapi0ebe->SetHistoPtRangeAndNBins(0, 15, 75) ;
  }
  else{
    nms->SetHistoERangeAndNBins(0, 30, 200) ;  
    anapi0ebe->SetHistoPtRangeAndNBins(0, 30, 100) ;
  }
  //      nms->SetHistoPtRangeAndNBins(0, 50, 100) ;
  //      nms->SetHistoAngleRangeAndNBins(0, 0.3, 100) ;
  //      nsm->SetHistoIMRangeAndNBins(0, 0.4, 100) ;  
  //Set Histrograms bins and ranges
  //      anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      anapi0->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anapi0ebe->Print("");
	
  //-------------------------------------
  //*** analysis the omega->pi0+gamma ***
  //------------------------------------
  AliAnaOmegaToPi0Gamma *anaomegaToPi0Gamma = new AliAnaOmegaToPi0Gamma();
  anaomegaToPi0Gamma->SetDebug(-1);//10 for lots of messages
  anaomegaToPi0Gamma->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  anaomegaToPi0Gamma->SetInputAODPhotonName(Form("Photons%s",calorimeter.Data()));
  anaomegaToPi0Gamma->SetNPID(1);
  anaomegaToPi0Gamma->SetNVtxZ(1);
  anaomegaToPi0Gamma->SetNEventsMixed(4);
  if(calorimeter=="PHOS")
    anaomegaToPi0Gamma->SetPi0MassPeakWidthCut(0.008); // PHOS
  else if(calorimeter=="EMCAL")
    anaomegaToPi0Gamma->SetPi0MassPeakWidthCut(0.012); // EMCAL 
  anaomegaToPi0Gamma->SetHistoPtRangeAndNBins(0, 20, 100) ;
  anaomegaToPi0Gamma->SetHistoMassRangeAndNBins(0, 1, 100) ;
  anaomegaToPi0Gamma->SetPi0OverOmegaPtCut(0.8);
  anaomegaToPi0Gamma->SetGammaOverOmegaPtCut(0.2);
  if(kUseKinematics) anaomegaToPi0Gamma->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else anaomegaToPi0Gamma->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anaomegaToPi0Gamma->AddToHistogramsName(Form("AnaOmegaToPi0Gamma%s_",calorimeter.Data()));
  if(kPrintSettings)   anaomegaToPi0Gamma->Print("");
	
	
//  //---------------------------------------------------------------------
//  // Electron/btag
//  //---------------------------------------------------------------------
//  if(calorimeter=="EMCAL"){
//    
//    AliAnaBtag *anabtag = new AliAnaBtag();
//    anabtag->SetDebug(-1); //10 for lots of messages
//    if(kUseKinematics){
//      anabtag->SwitchOnDataMC();
//      anabtag->SetMinPt(1.);
//    }
//    anabtag->SetOutputAODName("ElectronsEMCAL");
//    anabtag->SetOutputAODClassName("AliAODPWG4Particle");
//    //anabtag->SetHistoPtRangeAndNBins(0, 100, 100) ;
//    //anabtag->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//    //anabtag->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
//    if(kPrintSettings)anabtag->Print("");
//  }
  
  //==================================
  // ### Isolation analysis ###	
  //=================================
  //Photon
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  anaisol->SetMinPt(0);
  anaisol->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anaisol->SetAODObjArrayName("ICPhoton"); 
  anaisol->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaisol->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaisol->SwitchOffDataMC() ;
  //Select clusters with no pair, if both clusters with pi0 mass
  anaisol->SwitchOffInvariantMass();
  //Do isolation cut
  AliIsolationCut * ic =  anaisol->GetIsolationCut();	
  ic->SetConeSize(0.4);
  ic->SetPtThreshold(0.7);
  ic->SetPtFraction(0.1);
  ic->SetSumPtThreshold(1.0) ;
  ic->SetParticleTypeInCone(AliIsolationCut::kOnlyCharged);
  ic->SetICMethod(AliIsolationCut::kSumPtFracIC);
  if(kPrintSettings) ic->Print("");
  
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisol->SwitchOffReIsolation();
  //Multiple IC
  anaisol->SwitchOffSeveralIsolation() ;
  //Set Histograms bins and ranges
  anaisol->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  anaisol->AddToHistogramsName("AnaIsolPhoton_");
  if(kPrintSettings) anaisol->Print("");
  
  //Pi0
  AliAnaParticleIsolation *anaisolpi0 = new AliAnaParticleIsolation();
  anaisolpi0->SetDebug(-1);
  anaisolpi0->SetMinPt(0);
  anaisolpi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  anaisolpi0->AddToHistogramsName("AnaIsolPi0_");
  anaisolpi0->SetAODObjArrayName("ICPi0"); 
  anaisolpi0->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaisolpi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaisolpi0->SwitchOffDataMC() ;
  //Select clusters with no pair, if both clusters with pi0 mass
  anaisolpi0->SwitchOffInvariantMass();
  //Do isolation cut
  AliIsolationCut * ic2 =  anaisolpi0->GetIsolationCut();	
  ic2->SetConeSize(0.4);
  ic2->SetPtThreshold(0.7);
  ic2->SetPtFraction(0.1);
  ic2->SetSumPtThreshold(1.0) ;
  ic2->SetICMethod(AliIsolationCut::kSumPtFracIC);
  ic2->SetParticleTypeInCone(AliIsolationCut::kOnlyCharged);
  if(kPrintSettings) ic2->Print("");
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisolpi0->SwitchOffReIsolation();
  //Multiple IC
  anaisolpi0->SwitchOffSeveralIsolation() ;
  //Set Histograms bins and ranges
  anaisolpi0->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anaisolpi0->Print("");
	
  //===========================
  //Correlation analysis
  //===========================
	
  // ### Correlation with Jet Finder AOD output
  AliAnaParticleJetFinderCorrelation *anacorrjet = new AliAnaParticleJetFinderCorrelation();
  anacorrjet->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anacorrjet->SwitchOffFiducialCut();
  anacorrjet->SetDebug(-1);
  anacorrjet->SetConeSize(1);  
  anacorrjet->SelectIsolated(kTRUE); // do correlation with isolated photons
  anacorrjet->SetPtThresholdInCone(0.2);
  anacorrjet->SetDeltaPhiCutRange(0.5,5.5);//Mostly Open Cuts 
  anacorrjet->SetRatioCutRange(0.01,3); //Mostly Open Cuts
  anacorrjet->UseJetRefTracks(kFALSE); //Not working now
  //Set Histograms bins and ranges
  anacorrjet->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrjet->Print("");
  
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorrhadron = new AliAnaParticleHadronCorrelation();
  anacorrhadron->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anacorrhadron->AddToHistogramsName("AnaHadronCorrPhoton_");
  anacorrhadron->SetAODObjArrayName("PhotonHadronCorr"); 
  anacorrhadron->SetDebug(-1);
  anacorrhadron->SwitchOffCaloPID();
  if(kSimulation){
    anacorrhadron->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anacorrhadron->GetFiducialCut();
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
    fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
  }
  anacorrhadron->SwitchOnDecayCorr();
  anacorrhadron->SetMultiBin(1);
  anacorrhadron->SwitchOffNeutralCorr();
  anacorrhadron->SwitchOffEventSelection();
  anacorrhadron->SetPtCutRange(0.1,100);
  anacorrhadron->SetDeltaPhiCutRange(1.5,4.5);
  anacorrhadron->SwitchOnSeveralUECalculation();
  anacorrhadron->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  anacorrhadron->SelectIsolated(kFALSE); // do correlation with isolated photons
  if(kUseKinematics) anacorrhadron->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anacorrhadron->SwitchOffDataMC() ;
  //if(calorimeter=="PHOS"){
  //Correlate with particles in EMCAL
  //anacorrhadron->SwitchOnCaloPID();
  //anacorrhadron->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  //}
  //Set Histograms bins and ranges
  anacorrhadron->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrhadron->Print("");
  
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorrisohadron = new AliAnaParticleHadronCorrelation();
  anacorrisohadron->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anacorrisohadron->AddToHistogramsName("AnaHadronCorrIsoPhoton_");
  anacorrisohadron->SetAODObjArrayName("IsoPhotonHadronCorr"); 
  anacorrisohadron->SetDebug(-1);
  anacorrisohadron->SwitchOffCaloPID();
    if(kSimulation){
    anacorrisohadron->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anacorrisohadron->GetFiducialCut();
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
    fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
  }
  anacorrisohadron->SwitchOnDecayCorr();
  anacorrisohadron->SetMultiBin(1);
  anacorrisohadron->SwitchOffNeutralCorr();
  anacorrisohadron->SwitchOffEventSelection();
  anacorrisohadron->SetPtCutRange(0.1,100);
  anacorrisohadron->SetDeltaPhiCutRange(1.5,4.5);
  anacorrisohadron->SwitchOnSeveralUECalculation();
  anacorrisohadron->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  anacorrisohadron->SelectIsolated(kTRUE); // do correlation with isolated photons
  if(kUseKinematics) anacorrisohadron->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anacorrisohadron->SwitchOffDataMC() ;
  //if(calorimeter=="PHOS"){
  //Correlate with particles in EMCAL
  //anacorrhadron->SwitchOnCaloPID();
  //anacorrhadron->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  //}
  //Set Histograms bins and ranges
  anacorrisohadron->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrisohadron->Print("");
  
  
  // ### Pi0 Correlation with hadrons, not isolated
  AliAnaParticleHadronCorrelation *anacorrhadronpi0 = new AliAnaParticleHadronCorrelation();
  anacorrhadronpi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  anacorrhadronpi0->AddToHistogramsName("AnaHadronCorrPi0_");
  anacorrhadronpi0->SetAODObjArrayName("Pi0HadronCorr"); 
  anacorrhadronpi0->SetDebug(-1);
  anacorrhadronpi0->SwitchOffCaloPID();
  if(kSimulation){
    anacorrhadronpi0->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anacorrhadronpi0->GetFiducialCut();
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
    fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
  }
  anacorrhadronpi0->SwitchOnDecayCorr();
  anacorrhadronpi0->SetMultiBin(1);
  anacorrhadronpi0->SwitchOffNeutralCorr();
  anacorrhadronpi0->SwitchOffEventSelection();
  anacorrhadronpi0->SetPtCutRange(0.1,100);
  anacorrhadronpi0->SetDeltaPhiCutRange(1.5,4.5);
  anacorrhadronpi0->SelectIsolated(kFALSE); // do correlation with non isolated pi0
  anacorrhadronpi0->SwitchOnSeveralUECalculation();
  anacorrhadronpi0->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  if(kUseKinematics) anacorrhadronpi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anacorrhadronpi0->SwitchOffDataMC() ;
  //if(calorimeter=="PHOS"){
  //	//Correlate with particles in EMCAL
  //	anacorrhadronpi0->SwitchOnCaloPID();
  //	anacorrhadronpi0->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  //}
  //Set Histograms bins and ranges
  anacorrhadronpi0->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrhadronpi0->Print("");
  
  // ### Pi0 Correlation with hadrons, isolated
  AliAnaParticleHadronCorrelation *anacorrhadronisopi0 = new AliAnaParticleHadronCorrelation();
  anacorrhadronisopi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  anacorrhadronisopi0->AddToHistogramsName("AnaHadronCorrIsoPi0_");
  anacorrhadronisopi0->SetAODObjArrayName("IsoPi0HadronCorr"); 
  anacorrhadronisopi0->SetDebug(-1);
  anacorrhadronisopi0->SwitchOffCaloPID();
    if(kSimulation){
    anacorrhadronisopi0->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anacorrhadronisopi0->GetFiducialCut();
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
    fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
  }
  anacorrhadronisopi0->SwitchOnDecayCorr();
  anacorrhadronisopi0->SetMultiBin(1);
  anacorrhadronisopi0->SwitchOffNeutralCorr();
  anacorrhadronisopi0->SwitchOffEventSelection();
  anacorrhadronisopi0->SetPtCutRange(0.1,100);
  anacorrhadronisopi0->SetDeltaPhiCutRange(1.5,4.5);
  anacorrhadronisopi0->SelectIsolated(kTRUE); // do correlation with isolated pi0
  anacorrhadronisopi0->SwitchOnSeveralUECalculation();
  anacorrhadronisopi0->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  if(kUseKinematics) anacorrhadronisopi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anacorrhadronisopi0->SwitchOffDataMC() ;
  //if(calorimeter=="PHOS"){
  //	//Correlate with particles in EMCAL
  //	anacorrhadronpi0->SwitchOnCaloPID();
  //	anacorrhadronpi0->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
  //}
  //Set Histograms bins and ranges
  anacorrhadronisopi0->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrhadronisopi0->Print("");
  
  
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->SetCaloUtils(cu); //pointer to calorimeter utils
  Int_t n = 0;//Analysis number, order is important
  // Particle selection analysis
  maker->AddAnalysis(anaphoton,n++);
  maker->AddAnalysis(anapi0,n++);
  maker->AddAnalysis(anapi0ebe,n++);
  maker->AddAnalysis(anaomegaToPi0Gamma,n++);  
  //if(calorimeter=="EMCAL")maker->AddAnalysis(anabtag,n++);   
  // Isolation analysis
  maker->AddAnalysis(anaisol,n++);
  maker->AddAnalysis(anaisolpi0,n++);
  // Correlation analysis
  //maker->AddAnalysis(anacorrjet,n++);
  //maker->AddAnalysis(anacorrhadron,n++);
  //maker->AddAnalysis(anacorrhadronpi0,n++);
  //maker->AddAnalysis(anacorrisohadron,n++);
  //maker->AddAnalysis(anacorrhadronisopi0,n);
  maker->SetAnaDebug(0)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(inputDataType.Contains("delta")) maker->SwitchOffAODsMaker()  ;
  else                       maker->SwitchOnAODsMaker()  ;
	
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
  if(inputDataType=="ESD" && !kSimulation) task->SelectCollisionCandidates(); //AliPhysicsSelection has to be attached before.
  mgr->AddTask(task);
  
  //Create containers
  char name[128];
  sprintf(name,"PartCorr_%s",calorimeter.Data());
  cout<<"Name of task "<<name<<endl;
  //AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form(name),TList::Class(),
  //					   AliAnalysisManager::kOutputContainer, Form("PartCorr_%s.root",calorimeter.Data()));
  
  TString outputfile = AliAnalysisManager::GetCommonFileName(); 
  //  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form("PartCorr_%s",calorimeter.Data()),  TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PartCorr_%s",outputfile.Data(),calorimeter.Data()));
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(calorimeter.Data(), TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s:PartCorr",outputfile.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("%sCuts",calorimeter.Data()), TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             Form("%s:PartCorrCuts",outputfile.Data()));
	
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  if(!inputDataType.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  
  return task;
}


