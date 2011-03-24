AliAnalysisTaskParticleCorrelation *AddTaskPartCorr
(
 TString inputDataType, 
 TString calorimeter, 
 Bool_t kPrintSettings = kTRUE,
 Bool_t kSimulation = kFALSE, 
 Bool_t outputAOD=kFALSE, 
 Bool_t oldAOD=kFALSE,
 TString period
 ) {

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

  //-----------------------------------------------------------------
  // Switch on cluster energy smearing 
  //  -> SIMULATION
  //  -> EMCAL
  //-----------------------------------------------------------------

  if (kSimulation && calorimeter == "EMCAL") {
    //switch on cluster energy smearing
    reader->SwitchOnClusterEnergySmearing();
    reader->SetSmearingParameters(0,0.07);
    reader->SetSmearingParameters(0,0.00);
    reader->SetSmearingParameters(0,0.00);
  }

  
  //-----------------------------------------------------------------
  // Z vertex cut
  reader->SetZvertexCut(10.);
  //-----------------------------------------------------------------

  //-----------------------------------------------------------------
  // Min particle pT
  //-----------------------------------------------------------------
  reader->SetEMCALPtMin(0.3); 
  reader->SetPHOSPtMin(0.3);
  reader->SetCTSPtMin(0.1);

  if(outputAOD)      reader->SwitchOnWriteDeltaAOD()  ;
  if(oldAOD)         reader->SwitchOnOldAODs();
  if(kPrintSettings) reader->Print("");
  

  //-----------------------------------------------------------------
  // Bad cluster removal
  //  -> REAL DATA ONLY
  //-----------------------------------------------------------------
  reader->SwitchOnSuspiciousClustersRemoval();


  // *** Calorimeters Utils	***
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;

  //-----------------------------------------------------------------
  // Non-linearity corrections 
  //  -> REAL DATA AND SIMULATION
  //  -> EMCAL
  //-----------------------------------------------------------------
  if (calorimeter == "EMCAL") {
    cu->SwitchOnCorrectClusterLinearity();
    if (!kSimulation) {
      cu->GetEMCALRecoUtils()->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrected);
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(0,0.976       ) ;
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(1,9.83529e-01 ) ;
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(2,-1.84235e+02) ; 
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(3,-2.05019e+00) ;
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(4,-5.89423e+00) ;
    }
    else              {
      cu->GetEMCALRecoUtils()->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MC);
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(0,1.001   ) ;
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(1,-0.01264) ;
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(2,-0.03632) ;
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(3,0.1798  ) ;
      cu->GetEMCALRecoUtils()->SetNonLinearityParam(4,-0.522  ) ;
    }
  }
  
  //-----------------------------------------------------------------
  //  Remove clusters close to borders, 
  //  at least max energy cell is 1 cell away 
  //-----------------------------------------------------------------
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  
  //-----------------------------------------------------------------
  // Remove EMCAL hottest channels 
  //  -> REAL DATA AND SIMULATION
  //  -> EMCAL
  // Recover the file from alien  
  //   /alice/cern.ch/user/g/gconesab/BadChannelsDB
  //-----------------------------------------------------------------
  if (calorimeter == "EMCAL") {
    cu->SwitchOnBadChannelsRemoval();
    cu->SwitchOnDistToBadChannelRecalculation();
    TFile * fbad = new TFile("BadChannels.root","read");
    TH2I * hbad0 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod0");
    TH2I * hbad1 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod1");
    TH2I * hbad2 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod2");
    TH2I * hbad3 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod3");
    cu->SetEMCALChannelStatusMap(0,hbad0);
    cu->SetEMCALChannelStatusMap(1,hbad1);
    cu->SetEMCALChannelStatusMap(2,hbad2);
    cu->SetEMCALChannelStatusMap(3,hbad3);
  }
  
  
  //-----------------------------------------------------------------
  // Misalignment + recalculate position 
  //  -> REAL DATA and SIMULATION
  //  -> EMCAL
  //-----------------------------------------------------------------
  if (calorimeter == "EMCAL") {
    cu->SetEMCALGeometryName("EMCAL_FIRSTYEARV1");
    cu->GetEMCALRecoUtils()->SetParticleType(AliEMCALRecoUtils::kPhoton);
    cu->GetEMCALRecoUtils()->SetW0(4.5);      
    cu->GetEMCALRecoUtils()->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
    TGeoHMatrix *matrix[4];
    
    double rotationMatrix[4][9] = {-0.014587, -0.999892, -0.002031, 0.999892, -0.014591,  0.001979, -0.002009, -0.002002,  0.999996,
				   -0.014587,  0.999892,  0.002031, 0.999892,  0.014591, -0.001979, -0.002009,  0.002002, -0.999996,
				   -0.345864, -0.938278, -0.003412, 0.938276, -0.345874,  0.003010, -0.004004, -0.002161,  0.999990,
				   -0.345861,  0.938280,  0.003412, 0.938276,  0.345874, -0.003010, -0.004004,  0.002161, -0.999990};
    
    double translationMatrix[4][3] = {0.351659,    447.576446,  176.269742,
				      1.062577,    446.893974, -173.728870,
				      -154.213287, 419.306156,  176.753692,
				      -153.018950, 418.623681, -173.243605};
    for(int j=0; j<4; j++)
      {
	matrix[j] = new TGeoHMatrix();
	matrix[j]->SetRotation(rotationMatrix[j]);
	matrix[j]->SetTranslation(translationMatrix[j]);
	matrix[j]->Print();
	cu->SetEMCALGeometryMatrixInSM(matrix[j],j);
      }
    cu->SwitchOnRecalculateClusterTrackMatching();
  }
  
  //-----------------------------------------------------------------
  // Time dependent corrections 
  //  -> REAL_DATA ONLY
  //  -> EMCAL
  //  Recover file from alien  
  //  /alice/cern.ch/user/g/gconesab/TimeDepCorrectionDB
  //-----------------------------------------------------------------
  if (!kSimulation && calorimeter == "EMCAL") {
    cu->GetEMCALRecoUtils()->SwitchOnTimeDepCorrection();
    char cmd[200] ;
    sprintf(cmd, ".!tar xvfz CorrectionFiles.tgz >& /dev/null") ;
    gROOT->ProcessLine(cmd) ;
  }
  
  
  //------------------------------------------------------------------------
  //     Recalibration factors 
  //       -> REAL DATA ONLY
  //       -> EMCAL
  //     Recover the file from alien for LHC10d pass2
  //     /alice/cern.ch/user/g/gconesab/RecalDB/december2010        -> LHC10d pass2
  //     /alice/cern.ch/user/g/gconesab/RecalDB/summer_december2010 -> LHC10e pass1
  //  ******
  //     For other periods/passes, see
  //     https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalOffline#Summary_of_Calibration_and_Align
  //------------------------------------------------------------------------
  
  if (calorimeter == "EMCAL") {
    cu->SwitchOnRecalibration();
    TFile* f = 0x0 ;
    if (!kSimulation) {
      if      (period == "LHC10d") f = new TFile("RecalibrationFactors_LHC10d.root","read");
      else if (period == "LHC10e") f = new TFile("RecalibrationFactors_LHC10e.root","read");
      else                         Fatal("AddTaskPartCorr","run period not supported");
    }
    else {
      f = new TFile("DecalibrationFactors.root","read");
    }
    if (!f || !f->IsOpen()) Fatal("AddTaskPartCorr","Re(De)-calibration file not found");

    Info("AddTaskPartCorr",Form("Using calibration files for period %s",period.Data()));

    TH2F * h0 = (TH2F*)f->Get("EMCALRecalFactors_SM0");
    TH2F * h1 = (TH2F*)f->Get("EMCALRecalFactors_SM1");
    TH2F * h2 = (TH2F*)f->Get("EMCALRecalFactors_SM2");
    TH2F * h3 = (TH2F*)f->Get("EMCALRecalFactors_SM3");
    cu->SetEMCALChannelRecalibrationFactors(0,h0);
    cu->SetEMCALChannelRecalibrationFactors(1,h1);
    cu->SetEMCALChannelRecalibrationFactors(2,h2);
    cu->SetEMCALChannelRecalibrationFactors(3,h3);
  }

  cu->SetDebug(-1);


  // ##### Analysis algorithm settings ####
  
  // -------------------------------------------------
  // --- Photon/Pi0/Omega/Electron Analysis ---
  // -------------------------------------------------
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  //settings for different multiplicity analysis
  anaphoton->SwitchOffEventSelection() ;
  //anaphoton->SetMultiplicity(80, 120);

  if(calorimeter == "PHOS"){
    anaphoton->SetNCellCut(2);
    anaphoton->SetMinPt(0.);
    anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
  }
  else {//EMCAL
    anaphoton->SetNCellCut(1);
    anaphoton->SetMinPt(0.); 
    if(!kUseKinematics) anaphoton->SetTimeCut(400,900);// Time window of [400-900] ns
    anaphoton->SetMinDistanceToBadChannel(4, 5, 10);
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
  anaphoton->SetHistoPtRangeAndNBins(0, 20, 200) ;
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
  anapi0->SetHistoPtRangeAndNBins(0, 20, 40) ;
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
  //anapi0ebe->SetMultiplicity(80, 120);
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
    anapi0ebe->SetHistoPtRangeAndNBins(0, 30, 60) ;
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
  //maker->AddAnalysis(anapi0,n++);
  maker->AddAnalysis(anapi0ebe,n++);
//   maker->AddAnalysis(anaomegaToPi0Gamma,n++);  
  //if(calorimeter=="EMCAL")maker->AddAnalysis(anabtag,n++);   
  // Isolation analysis
  maker->AddAnalysis(anaisol,n++);
//   maker->AddAnalysis(anaisolpi0,n++);
  // Correlation analysis
//   maker->AddAnalysis(anacorrjet,n++);
  maker->AddAnalysis(anacorrhadron,n++);
//   maker->AddAnalysis(anacorrhadronpi0,n++);
//   maker->AddAnalysis(anacorrisohadron,n++);
//   maker->AddAnalysis(anacorrhadronisopi0,n);
  maker->SetAnaDebug(0)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(inputDataType.Contains("delta")) maker->SwitchOffAODsMaker()  ;
  else                                maker->SwitchOnAODsMaker()  ;
	
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
