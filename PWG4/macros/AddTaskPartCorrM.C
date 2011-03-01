AliAnalysisTaskParticleCorrelationM *AddTaskPartCorrM(TString data, TString calorimeter, Bool_t kPrintSettings = kFALSE)
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
  kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;    
	
  cout<<"********* ACCESS KINE? "<<kUseKinematics<<endl;
  
  // Configure analysis
  //===========================================================================
  
  //Reader
  AliCaloTrackReader * reader = 0x0;
  if(data.Contains("AOD")) reader = new AliCaloTrackAODReader();
  else if(data=="ESD") reader = new AliCaloTrackESDReader();
  else if(data=="MC" || (kUseKinematics && data == "ESD")) reader = new AliCaloTrackMCReader();
  reader->SetDebug(-1);//10 for lots of messages
  reader->SwitchOnCTS();
  
  if(calorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();  
    reader->SwitchOnEMCAL(); 
  }
  if(calorimeter == "PHOS") { 
    reader->SwitchOnPHOSCells();  
    reader->SwitchOnPHOS();
  }
  
   reader->SwitchOnSuspiciousClustersRemoval();  //EMCAL

  // for case data="deltaAOD", no need to fill the EMCAL/PHOS cluster lists
  if(data.Contains("delta")){
    reader->SwitchOffEMCAL();
    reader->SwitchOffPHOS();
    reader->SwitchOffEMCALCells(); 
    reader->SwitchOffPHOSCells(); 
  }
	
  if(kUseKinematics){
    if(data == "ESD"){
      reader->SwitchOnStack();          
      reader->SwitchOffAODMCParticles(); 
    }
    else if(data == "AOD"){
      reader->SwitchOffStack();          
      reader->SwitchOnAODMCParticles(); 
    }
  }
  
  reader->SetZvertexCut(10.);
  
  //Min particle pT
  reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  reader->SetCTSPtMin(0.1);
  reader->SwitchOffWriteDeltaAOD()  ;
  if(kPrintSettings) reader->Print("");
  
  // ##### Analysis algorithm settings ####
  AliCaloPID * pid = new AliCaloPID();
  pid->SetDispersionCut(1.5);
  pid->SetTOFCut(5.e-9);
  pid->SetDebug(-1);
  if(kPrintSettings) pid->Print("");
  //Fiducial cut  
//  AliFiducialCut * fidCut1stYear = new AliFiducialCut();
//  fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
//  fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);
//  fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
//  fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
//  fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
//  fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,300.);
  // -------------------------------------------------
  // --- Isolation Cut ---
  // -------------------------------------------------
  AliIsolationCut * ic = new AliIsolationCut();
  ic->SetConeSize(0.4);
  //ic->SetPtThreshold(0.7);
  ic->SetPtFraction(0.1);
  ic->SetPtThreshold(0.7) ;
  ic->SetSumPtThreshold(1.0) ;
  //choose different method for IC:
  //kPtThresIC, kSumPtIC, kPtFracIC, kSumPtFracIC
  //  ic->SetICMethod(AliIsolationCut::kPtThresIC);
  ic->SetICMethod(AliIsolationCut::kSumPtFracIC);
  //particle in cone: kNeutralAndCharged=0, kOnlyNeutral=1, kOnlyCharged=2
  ic->SetParticleTypeInCone(AliIsolationCut::kOnlyCharged);
  if(kPrintSettings) ic->Print("");
  
  //analysis with calorimeter triggers  
  if(calorimeter=="PHOS" || calorimeter=="EMCAL") {  
    AliCalorimeterUtils * cu = new AliCalorimeterUtils();
    cu->SwitchOnBadChannelsRemoval();
    
    cu->SetNumberOfCellsFromEMCALBorder(1) ; //nEMCAL);
    cu->SetNumberOfCellsFromPHOSBorder(2) ; //nPHOS);
   // cu->SwitchOnNoFiducialBorderInEMCALEta0();
    cu->SetPHOSChannelStatus(1,48, 8); //PHOS new hot channel
    
    // // Remove EMCAL hottest channels from Gustavo list
    // SM0
    cu->SetEMCALChannelStatus(0,3,13);  cu->SetEMCALChannelStatus(0,44,1); cu->SetEMCALChannelStatus(0,3,13); //warm
    cu->SetEMCALChannelStatus(0,20,7);  cu->SetEMCALChannelStatus(0,38,2);   //hot
    // SM1 warm channels
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
    
    // -----------------------------------
    // --- Photon and Pi0  Analysis ---
    // -----------------------------------
    
    AliAnaPhoton *anaphoton = new AliAnaPhoton();
    anaphoton->SetDebug(-1); //10 for lots of messages
    anaphoton->SetMinPt(1.0);
    anaphoton->SetCaloPID(pid);
    anaphoton->SetCalorimeter(calorimeter);
    if(kUseKinematics) anaphoton->SwitchOnDataMC() ;//Access MC stack and fill more histograms
    else  anaphoton->SwitchOffDataMC() ;
    anaphoton->SwitchOffCaloPID();
    anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
    anaphoton->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anaphoton->GetFiducialCut();
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
    anaphoton->SwitchOnTrackMatchRejection();
    //settings for different multiplicity analysis
    anaphoton->SwitchOffEventSelection() ;
    anaphoton->SetMultiplicity(80, 120);
    
    if(calorimeter == "EMCAL"){
      anaphoton->SetNCellCut(1);
      anaphoton->SetMinPt(0.3); 
      //if(!kUseKinematics) anaphoton->SetTimeCut(525, 725);
      anaphoton->SetMinDistanceToBadChannel(1, 2, 3);
    }
    else{
      anaphoton->SetMinPt(0.3);       
      anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
      anaphoton->SetNCellCut(2);
    }
    if(!data.Contains("delta")) {
      anaphoton->SetOutputAODName(Form("Triggers%s",calorimeter.Data()));
      anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    }
    else anaphoton->SetOutputAODName(Form("Triggers%s",calorimeter.Data()));
    anaphoton->AddToHistogramsName("AnaPhoton_");
    //Set Histograms bins and ranges
    anaphoton->SetHistoPtRangeAndNBins(0, 50, 100) ;
    //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
    //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
    
    if(kPrintSettings) anaphoton->Print("");
    AliAnaPi0 *anapi0 = new AliAnaPi0();
    anapi0->SetDebug(-1);//10 for lots of messages
    anapi0->SetInputAODName(Form("Triggers%s",calorimeter.Data()));
    anapi0->AddToHistogramsName("AnaPi0_");
    anapi0->SetCaloPID(pid);
    anapi0->SetCalorimeter(calorimeter);
    anapi0->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anapi0->GetFiducialCut();
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
    anapi0->SetNPID(1); 
    //settings for different multiplicity analysis
    anapi0->SwitchOnEventSelection() ;
    anapi0->SetNCentrBin(1);
    anapi0->SetMultiplicity(80, 120);
    anapi0->SetMultiBin(1);  
    anapi0->SwitchOffDataMC() ;//Access MC stack and fill more histograms
    if(calorimeter=="PHOS") anapi0->SetNumberOfModules(3); //PHOS first year
    else if(calorimeter=="EMCAL") anapi0->SetNumberOfModules(4); //EMCAL first year
    anapi0->SetHistoPtRangeAndNBins(0, 50, 100) ;
    //anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
    //anapi0->SetHistoEtaRangeAndNBins(-0.8, 0.8, 200) ;
    anapi0->SetHistoMassRangeAndNBins(0., 1.0, 100) ;
    anapi0->SetHistoAsymmetryRangeAndNBins(0., 1. , 10) ;
    if(kPrintSettings) anapi0->Print("");
    
    // -------------------------------------------------
    // --- Pi0 EbE Analysis ---
    // -------------------------------------------------
    
    AliNeutralMesonSelection *nms = new AliNeutralMesonSelection();
    nms->SetInvMassCutRange(0.05, 0.16)     ;
    nms->SwitchOffAngleSelection() ;
    nms->KeepNeutralMesonSelectionHistos(kTRUE);
    //Set Histrograms bins and ranges
    nms->SetHistoERangeAndNBins(0, 50, 100) ;
    nms->SetHistoPtRangeAndNBins(0, 50, 100) ;
    //      nms->SetHistoAngleRangeAndNBins(0, 0.3, 100) ;
    //      nsm->SetHistoIMRangeAndNBins(0, 0.4, 100) ;  
    
    AliAnaPi0EbE *anapi0ebe = new AliAnaPi0EbE();
    anapi0ebe->SetDebug(-1);//10 for lots of messages
    anapi0ebe->SetAnalysisType(AliAnaPi0EbE::kIMCalo);
    anapi0ebe->SetMinPt(0.5);
    anapi0ebe->SetCalorimeter(calorimeter);
    anapi0ebe->SetInputAODName(Form("Triggers%s",calorimeter.Data()));
    if(!data.Contains("delta")) {
      anapi0ebe->SetOutputAODName(Form("Pi0s%s",calorimeter.Data()));
      anapi0ebe->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    }
    else  anapi0ebe->SetOutputAODName(Form("Pi0s%s",calorimeter.Data()));
    
    if(kUseKinematics) anapi0ebe->SwitchOnDataMC() ;//Access MC stack and fill more histograms
    else  anapi0ebe->SwitchOffDataMC() ;	
    anapi0ebe->SwitchOffEventSelection() ;
    anapi0ebe->SetMultiplicity(80, 120);
    anapi0ebe->SetMultiBin(1);  
    anapi0ebe->SetNeutralMesonSelection(nms);
    //Set Histrograms bins and ranges
    anapi0ebe->SetHistoPtRangeAndNBins(0, 50, 100) ;
    //      anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
    //      anapi0->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
    if(kPrintSettings) anapi0ebe->Print("");
    
    // ### Pi0 Correlation with hadrons, not isolated
    AliAnaParticleHadronCorrelation *anacorrhadronpi0 = new AliAnaParticleHadronCorrelation();
    anacorrhadronpi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
    anacorrhadronpi0->AddToHistogramsName("AnaHadronCorrPi0_");
    anacorrhadronpi0->SetDebug(-1);
    anacorrhadronpi0->SwitchOffCaloPID();
    anacorrhadronpi0->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anacorrhadronpi0->GetFiducialCut();
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
    fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
    anacorrhadronpi0->SwitchOnDecayCorr();
    anacorrhadronpi0->SetPtCutRange(0.5,50);
    anacorrhadronpi0->SetDeltaPhiCutRange(1.5,4.5);
    anacorrhadronpi0->SelectIsolated(kFALSE); // do correlation with non isolated pi0
    anacorrhadronpi0->SetMultiplicity(80, 100);
    anacorrhadronpi0->SetMultiBin(1);
    anacorrhadronpi0->SwitchOffNeutralCorr();
    anacorrhadronpi0->SwitchOffEventSelection();
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
    anacorrhadronpi0->SetHistoPtRangeAndNBins(0, 50, 100) ;
    //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
    //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
    if(kPrintSettings) anacorrhadronpi0->Print("");
    
    AliAnaParticleIsolation *anaisolpi0 = new AliAnaParticleIsolation();
    anaisolpi0->SetDebug(-1);
    anaisolpi0->SetMinPt(2.);
    anaisolpi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
    anaisolpi0->AddToHistogramsName("AnaIsolPi0_");
    anaisolpi0->SetCalorimeter(calorimeter);
    if(kUseKinematics) anaisolpi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
    else  anaisolpi0->SwitchOffDataMC() ;
    //Select clusters with no pair, if both clusters with pi0 mass
    anaisolpi0->SwitchOffInvariantMass();
    //anaisol->SetNeutralMesonSelection(nms);
    //Do isolation cut
    anaisolpi0->SetIsolationCut(ic);	
    //Do or not do isolation with previously produced AODs.
    //No effect if use of SwitchOnSeveralIsolation()
    anaisolpi0->SwitchOffReIsolation();
    //Multiple IC
    anaisolpi0->SwitchOffSeveralIsolation() ;
    //Set Histograms bins and ranges
    anaisolpi0->SetHistoPtRangeAndNBins(0, 50, 100) ;
    //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
    //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
    if(kPrintSettings) anaisol->Print("");
    
    // ### Pi0 Correlation with hadrons, isolated
    AliAnaParticleHadronCorrelation *anacorrhadronisopi0 = new AliAnaParticleHadronCorrelation();
    anacorrhadronisopi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
    anacorrhadronisopi0->AddToHistogramsName("AnaHadronCorrIsoPi0_");
    anacorrhadronisopi0->SetDebug(-1);
    anacorrhadronisopi0->SwitchOffCaloPID();
    anacorrhadronisopi0->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anacorrhadronisopi0->GetFiducialCut();
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
    fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
    anacorrhadronisopi0->SwitchOnDecayCorr();
    anacorrhadronisopi0->SetPtCutRange(0.5,50);
    anacorrhadronisopi0->SetDeltaPhiCutRange(1.5,4.5);
    anacorrhadronisopi0->SelectIsolated(kTRUE); // do correlation with isolated pi0
    anacorrhadronisopi0->SetMultiplicity(80, 100);
    anacorrhadronisopi0->SetMultiBin(1);
    anacorrhadronisopi0->SwitchOffNeutralCorr();
    anacorrhadronisopi0->SwitchOffEventSelection();  
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
    anacorrhadronisopi0->SetHistoPtRangeAndNBins(0, 50, 100) ;
    //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
    //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
    if(kPrintSettings) anacorrhadronisopi0->Print("");
    
  } //analysis in calorimeter
  else {
    //---charge particle trigger---------------------
    //------------------------------------------------
    
    AliAnaChargedParticles *anacharge = new AliAnaChargedParticles();
    anacharge->SetDebug(-1); //10 for lots of messages
    anacharge->SetMinPt(1.0);
    // anacharge->SetCaloPID(pid);
    if(kUseKinematics) anacharge->SwitchOnDataMC() ;//Access MC stack and fill more histograms
    else  anacharge->SwitchOffDataMC() ;
    anacharge->SwitchOffCaloPID();
    anacharge->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
    anacharge->SwitchOnFiducialCut();
    AliFiducialCut * fidCut1stYear = anacharge->GetFiducialCut();
    fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);       
    if(!data.Contains("delta")) {
      anacharge->SetOutputAODName(Form("Triggers%s",calorimeter.Data()));
      anacharge->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    }
    else anacharge->SetOutputAODName(Form("Triggers%s",calorimeter.Data()));
    anacharge->AddToHistogramsName("AnaCharge_");
    //Set Histograms bins and ranges
    anacharge->SetHistoPtRangeAndNBins(0, 50, 100) ;
    //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
    //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
    
    if(kPrintSettings) anacharge->Print("");    
  }
  
  // -------------------------------------------------
  // ---  Correlation Analysis with non-isolated triggers ---
  // -------------------------------------------------
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorrhadron = new AliAnaParticleHadronCorrelation();
  anacorrhadron->SetInputAODName(Form("Triggers%s",calorimeter.Data()));
  anacorrhadron->AddToHistogramsName("AnaHadronCorrTrig_");
  anacorrhadron->SetDebug(-1);
  anacorrhadron->SwitchOnCaloPID();
  anacorrhadron->SwitchOnFiducialCut();
  AliFiducialCut * fidCut1stYear = anacorrhadron->GetFiducialCut();
  fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
  fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
  fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
  fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
  fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
  fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
  anacorrhadron->SetPtCutRange(0.5,50);
  anacorrhadron->SetDeltaPhiCutRange(1.5,4.5);
  anacorrhadron->SetMultiplicity(80, 120);
  anacorrhadron->SetMultiBin(1);
  anacorrhadron->SwitchOffNeutralCorr();
  anacorrhadron->SwitchOffEventSelection();
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
  anacorrhadron->SetHistoPtRangeAndNBins(0, 50, 100) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrhadron->Print("");
  
  
  
  // ### Isolation analysis ###	
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  anaisol->SetMinPt(2.0);
  anaisol->SetInputAODName(Form("Triggers%s",calorimeter.Data()));
  anaisol->AddToHistogramsName("AnaIsolTrig_");
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
  //Set Histograms bins and ranges
  anaisol->SetHistoPtRangeAndNBins(0, 50, 100) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  
  if(kPrintSettings) anaisol->Print("");
  
  // -------------------------------------------------
  // ---  Correlation Analysis with isolated triggers ---
  // -------------------------------------------------
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorrisohadron = new AliAnaParticleHadronCorrelation();
  anacorrisohadron->SetInputAODName(Form("Triggers%s",calorimeter.Data()));
  anacorrisohadron->AddToHistogramsName("AnaHadronCorrIsoTrig_");
  anacorrisohadron->SetDebug(-1);
  anacorrisohadron->SwitchOffCaloPID();
  anacorrisohadron->SwitchOnFiducialCut();
  AliFiducialCut * fidCut1stYear = anacorrisohadron->GetFiducialCut();
  fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
  fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
  fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
  fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
  fidCut1stYear->DoCTSFiducialCut(kTRUE) ;
  fidCut1stYear->SetSimpleCTSFiducialCut(0.8,0.,360.);    
  anacorrisohadron->SetPtCutRange(0.5,50);
  anacorrisohadron->SetDeltaPhiCutRange(1.5,4.5);
  anacorrisohadron->SetMultiplicity(80, 100);
  anacorrisohadron->SetMultiBin(1);
  anacorrisohadron->SwitchOffNeutralCorr();
  anacorrisohadron->SwitchOffEventSelection();
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
  anacorrisohadron->SetHistoPtRangeAndNBins(0, 50, 100) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrisohadron->Print("");
  
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  Int_t n = 0;//Analysis number, order is important
  // Particle selection analysis
  if(calorimeter=="PHOS" || calorimeter=="EMCAL") {   
    maker->SetCaloUtils(cu);  //pointer to calorimeter utils
    maker->AddAnalysis(anaphoton,n++);
    maker->AddAnalysis(anapi0,n++);
    maker->AddAnalysis(anapi0ebe,n++);
    maker->AddAnalysis(anacorrhadronpi0,n++);
    maker->AddAnalysis(anaisolpi0,n++);
    maker->AddAnalysis(anacorrhadronisopi0,n++);
  }
   if(calorimeter=="CTS") 
     maker->AddAnalysis(anacharge,n++);
  // Correlation analysis
  maker->AddAnalysis(anacorrhadron,n++);
  // Isolation analysis
  maker->AddAnalysis(anaisol,n++);
  // Correlation analysis with isolated triggers
  maker->AddAnalysis(anacorrisohadron,n);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(data.Contains("delta")) maker->SwitchOffAODsMaker()  ;
  else                       maker->SwitchOnAODsMaker()  ;
	
  if(kPrintSettings) maker->Print("");
  
  printf("======================== \n");
  printf(" End Configuration of PartCorr analysis with detector %s \n",calorimeter.Data());
  printf("======================== \n");
  
  // Create task
  //===========================================================================
 // AliAnalysisTaskParticleCorrelationM * task = new AliAnalysisTaskParticleCorrelationM(Form("PartCorr%s",calorimeter.Data()));
  AliAnalysisTaskParticleCorrelationM * task = new AliAnalysisTaskParticleCorrelationM("PartCorr");
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  //task->SetDebugLevel(-1);
  task->SetAnalysisMaker(maker);
  //if(!kSimulation)task->SelectCollisionCandidates(); //AliPhysicsSelection has to be attached before.
  mgr->AddTask(task);
  
//  char name[128];
//  sprintf(name,"PartCorr_%s",calorimeter.Data());
//  cout<<"Name of task "<<name<<endl;
  //AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form(name),TList::Class(),
  //					   AliAnalysisManager::kOutputContainer, Form("PartCorr_%s.root",calorimeter.Data()));
  
  TString outputfile = AliAnalysisManager::GetCommonFileName(); 
  outputfile.ReplaceAll(".root","") ;
  outputfile.Append("M.root") ;  
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
 //if(!data.Contains("delta")) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
 // mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);  
  
  return task;
}


