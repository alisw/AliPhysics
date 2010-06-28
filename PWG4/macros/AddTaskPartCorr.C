AliAnalysisTaskParticleCorrelation *AddTaskPartCorr(TString data, TString calorimeter, Bool_t kPrintSettings = kFALSE,Bool_t kSimulation = kFALSE)
{
  // Creates a PartCorr task, configures it and adds it to the analysis manager.
  
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
  
  //Min particle pT
  reader->SetEMCALPtMin(0.1); 
  reader->SetPHOSPtMin(0.);
  reader->SetCTSPtMin(0.);
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
  // cu->SetEMCALChannelStatus(0,3,13);  cu->SetEMCALChannelStatus(0,44,1); cu->SetEMCALChannelStatus(0,3,13); 
  cu->SetEMCALChannelStatus(0,20,7);  cu->SetEMCALChannelStatus(0,38,2);   
  // SM1
  // cu->SetEMCALChannelStatus(1,4,7);   cu->SetEMCALChannelStatus(1,4,13);  cu->SetEMCALChannelStatus(1,9,20); 
  // cu->SetEMCALChannelStatus(1,14,15); cu->SetEMCALChannelStatus(1,23,16); cu->SetEMCALChannelStatus(1,32,23); 
  // cu->SetEMCALChannelStatus(1,37,5);  cu->SetEMCALChannelStatus(1,40,1);  cu->SetEMCALChannelStatus(1,40,2);
  // cu->SetEMCALChannelStatus(1,40,5);  cu->SetEMCALChannelStatus(1,41,0);  cu->SetEMCALChannelStatus(1,41,1);
  // cu->SetEMCALChannelStatus(1,41,2);  cu->SetEMCALChannelStatus(1,41,4);
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
  
  
  AliFiducialCut * fidCut1stYear = new AliFiducialCut();
  fidCut1stYear->DoCTSFiducialCut(kFALSE) ;
  if(kSimulation){
    fidCut1stYear->DoEMCALFiducialCut(kTRUE) ;
    fidCut1stYear->DoPHOSFiducialCut(kTRUE) ;
    fidCut1stYear->SetSimpleEMCALFiducialCut(0.7,80.,120.);
    fidCut1stYear->SetSimplePHOSFiducialCut(0.12,260.,320.);
  } 
  else{
    fidCut1stYear->DoEMCALFiducialCut(kFALSE) ;
    fidCut1stYear->DoPHOSFiducialCut(kFALSE) ;
  }	
  
  
  // -------------------------------------------------
  // --- Photon Isolation and Correlation Analysis ---
  // -------------------------------------------------
  
  AliAnaPhoton *anaphoton = new AliAnaPhoton();
  anaphoton->SetDebug(-1); //10 for lots of messages
  if(calorimeter == "PHOS"){
    anaphoton->SetNCellCut(1);// At least 2 cells
    anaphoton->SetMinPt(0.2);
    anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
  }
  else {//EMCAL
    //anaphoton->SetNCellCut(0);// At least 2 cells
    anaphoton->SetMinPt(0.1); // no effect minium EMCAL cut.
    anaphoton->SetTimeCut(525,725);// Time window of [550-750] ns
    anaphoton->SetMinDistanceToBadChannel(6, 12, 18);
  }
  anaphoton->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaphoton->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaphoton->SwitchOffDataMC() ;
  anaphoton->SwitchOffCaloPID();
  anaphoton->SwitchOffFiducialCut();
  if(kSimulation){
    anaphoton->SwitchOnFiducialCut();
    anaphoton->SetFiducialCut(fidCut1stYear);
  }
  
  if(!data.Contains("delta")) {
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
  if(kSimulation){
    anapi0->SwitchOnFiducialCut();
    anapi0->SetFiducialCut(fidCut1stYear);
  }  
  anapi0->SetNPID(1); //Available from tag AliRoot::v4-18-15-AN
  //settings for pp collision
  anapi0->SetNCentrBin(1);
  anapi0->SetNZvertBin(1);
  anapi0->SetNRPBin(1);
  anapi0->SetNMaxEvMix(10);
  anapi0->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  if(calorimeter=="PHOS") anapi0->SetNumberOfModules(3); //PHOS first year
  else  anapi0->SetNumberOfModules(4); //EMCAL first year
  anapi0->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //anapi0->SetHistoEtaRangeAndNBins(-0.8, 0.8, 200) ;
  anapi0->SetHistoMassRangeAndNBins(0., 0.6, 200) ;
  anapi0->SetHistoAsymmetryRangeAndNBins(0., 1. , 10) ;
  if(kPrintSettings) anapi0->Print("");
  
  // ### Isolation analysis ###	
  
  AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
  anaisol->SetDebug(-1);
  anaisol->SetMinPt(0);
  anaisol->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anaisol->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaisol->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaisol->SwitchOffDataMC() ;
  //Select clusters with no pair, if both clusters with pi0 mass
  anaisol->SwitchOffInvariantMass();
  //anaisol->SetNeutralMesonSelection(nms);
  //Do isolation cut
  AliIsolationCut * ic =  anaisol->GetIsolationCut();	
  ic->SetConeSize(0.4);
  ic->SetPtThreshold(0.2);
  ic->SetICMethod(AliIsolationCut::kPtThresIC);
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
  anacorrhadron->SetDebug(-1);
  anacorrhadron->SwitchOffCaloPID();
  anacorrhadron->SwitchOffFiducialCut();
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
  anacorrisohadron->SetDebug(-1);
  anacorrisohadron->SwitchOffCaloPID();
  anacorrisohadron->SwitchOffFiducialCut();
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
  
  // -------------------------------------------------
  // --- Pi0    Isolation and Correlation Analysis ---
  // -------------------------------------------------
	
  AliNeutralMesonSelection *nms = new AliNeutralMesonSelection();
  nms->SetInvMassCutRange(0.05, 0.2)     ;
  nms->KeepNeutralMesonSelectionHistos(kTRUE);
  //Set Histrograms bins and ranges
  nms->SetHistoERangeAndNBins(0, 50, 200) ;
  //      nms->SetHistoPtRangeAndNBins(0, 50, 100) ;
  //      nms->SetHistoAngleRangeAndNBins(0, 0.3, 100) ;
  //      nsm->SetHistoIMRangeAndNBins(0, 0.4, 100) ;  
  
  AliAnaPi0EbE *anapi0ebe = new AliAnaPi0EbE();
  anapi0ebe->SetDebug(-1);//10 for lots of messages
  anapi0ebe->SetAnalysisType(AliAnaPi0EbE::kIMCalo);
  anapi0ebe->SetMinPt(0);
  anapi0ebe->SetCalorimeter(calorimeter);
  anapi0ebe->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  if(!data.Contains("delta")) {
    anapi0ebe->SetOutputAODName(Form("Pi0s%s",calorimeter.Data()));
    anapi0ebe->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else  anapi0ebe->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  
  if(kUseKinematics) anapi0ebe->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anapi0ebe->SwitchOffDataMC() ;	
  anapi0ebe->SetNeutralMesonSelection(nms);
  //Set Histrograms bins and ranges
  anapi0ebe->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      anapi0->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anapi0ebe->Print("");
  
  AliAnaParticleIsolation *anaisolpi0 = new AliAnaParticleIsolation();
  anaisolpi0->SetDebug(-1);
  anaisolpi0->SetMinPt(0);
  anaisolpi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  anaisolpi0->AddToHistogramsName("AnaIsolPi0_");
  anaisolpi0->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaisolpi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaisolpi0->SwitchOffDataMC() ;
  //Select clusters with no pair, if both clusters with pi0 mass
  anaisolpi0->SwitchOffInvariantMass();
  //anaisol->SetNeutralMesonSelection(nms);
  //Do isolation cut
  AliIsolationCut * ic2 =  anaisolpi0->GetIsolationCut();	
  ic2->SetConeSize(0.4);
  ic2->SetPtThreshold(0.2);
  ic2->SetICMethod(AliIsolationCut::kPtThresIC);
  if(kPrintSettings) ic->Print("");
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisolpi0->SwitchOffReIsolation();
  //Multiple IC
  anaisolpi0->SwitchOffSeveralIsolation() ;
  //Set Histograms bins and ranges
  anaisolpi0->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anaisol->Print("");
  
  
  // ### Pi0 Correlation with hadrons, not isolated
  AliAnaParticleHadronCorrelation *anacorrhadronpi0 = new AliAnaParticleHadronCorrelation();
  anacorrhadronpi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  anacorrhadronpi0->AddToHistogramsName("AnaHadronCorrPi0_");
  anacorrhadronpi0->SetDebug(-1);
  anacorrhadronpi0->SwitchOffCaloPID();
  anacorrhadronpi0->SwitchOffFiducialCut();
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
  anacorrhadronisopi0->SetDebug(-1);
  anacorrhadronisopi0->SwitchOffCaloPID();
  anacorrhadronisopi0->SwitchOffFiducialCut();
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
  
  //*** analysis the omega->pi0+gamma ***
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
  // Isolation analysis
  maker->AddAnalysis(anaisol,n++);
  maker->AddAnalysis(anacorrisohadron,n++);
  maker->AddAnalysis(anaisolpi0,n++);
  // Correlation analysis
  maker->AddAnalysis(anacorrjet,n++);
  maker->AddAnalysis(anacorrhadron,n++);
  maker->AddAnalysis(anacorrhadronpi0,n++);
  maker->AddAnalysis(anacorrhadronisopi0,n);
  maker->SetAnaDebug(0)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(data.Contains("delta")) maker->SwitchOffAODsMaker()  ;
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
  task->SelectCollisionCandidates();
  task->SetAnalysisMaker(maker);
  //if(!kSimulation)task->SelectCollisionCandidates(); //AliPhysicsSelection has to be attached before.
  mgr->AddTask(task);
  
  char name[128];
  sprintf(name,"PartCorr_%s",calorimeter.Data());
  cout<<"Name of task "<<name<<endl;
  //AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form(name),TList::Class(),
  //					   AliAnalysisManager::kOutputContainer, Form("PartCorr_%s.root",calorimeter.Data()));
  
  TString outputfile = AliAnalysisManager::GetCommonFileName(); 
  //  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(Form("PartCorr_%s",calorimeter.Data()),  TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PartCorr_%s",outputfile.Data(),calorimeter.Data()));
  AliAnalysisDataContainer *cout_pc = mgr->CreateContainer(calorimeter.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PartCorr",outputfile.Data()));
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  if(!data.Contains("delta")) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  
  return task;
}


