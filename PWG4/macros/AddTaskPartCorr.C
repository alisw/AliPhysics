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
   
   //Reader
   AliCaloTrackReader * reader = 0x0;
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
  
  // ##### Analysis algorithm settings ####
  
  
  AliCaloPID * pid = new AliCaloPID();
  pid->SetDispersionCut(1.5);
  pid->SetTOFCut(5.e-9);
  pid->SetDebug(-1);
  if(kPrintSettings) pid->Print("");
	
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
  
	// --------------------
	// --- QA Analysis ---
	// --------------------
	
//  AliAnaCalorimeterQA *qa = new AliAnaCalorimeterQA();
//  //qa->SetDebug(10); //10 for lots of messages
//  qa->SetCalorimeter(calorimeter);
//  if(kUseKinematics) qa->SwitchOnDataMC() ;//Access MC stack or AODMCOParticles
//  else  qa->SwitchOffDataMC() ;
//  qa->AddToHistogramsName("AnaCaloQA_");
//  if(kSimulation){
//    qa->SetFiducialCut(fidCut1stYear);
//    qa->SwitchOnFiducialCut();
//  }
//  if(qa=="PHOS") qa->SetNumberOfModules(3); //PHOS first year
//  else  qa->SetNumberOfModules(4); //EMCAL first year
//  //Set Histograms bins and ranges
//  qa->SetHistoPtRangeAndNBins(0, 50, 500) ;
//  qa->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//  if(calorimeter == "PHOS"){
//    qa->SetHistoEtaRangeAndNBins(-0.13, 0.13, 100) ;
//    qa->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 330*TMath::DegToRad() , 200) ;
//  }
//  else if(calorimeter == "EMCAL"){
//    qa->SetHistoEtaRangeAndNBins(-0.8, 0.8, 200) ;
//    qa->SetHistoPhiRangeAndNBins(70*TMath::DegToRad(), 130*TMath::DegToRad(), 200) ;
//  }
//  qa->SetHistoMassRangeAndNBins(0., 0.6, 300) ;
//  qa->SetHistoAsymmetryRangeAndNBins(0., 1. , 25) ;
//  qa->SetHistoPOverERangeAndNBins(0,10.,100);
//  qa->SetHistodEdxRangeAndNBins(0.,400.,200);
//  qa->SetHistodRRangeAndNBins(0.,TMath::Pi(),300);
//  qa->SetHistoTimeRangeAndNBins(0.,1000,1000);
//  qa->SetHistoRatioRangeAndNBins(0.,2.,100);
//  qa->SetHistoVertexDistRangeAndNBins(0.,500.,100);
//  qa->SetHistoNClusterCellRangeAndNBins(0,300,300);
//	
//  if(kPrintSettings) qa->Print("");	
  
	
  // -----------------------------------
  // --- Pi0 Invariant Mass Analysis ---
  // -----------------------------------

  AliAnaPhoton *anaphoton1 = new AliAnaPhoton();
  anaphoton1->SetDebug(-1); //10 for lots of messages
  if(calorimeter == "PHOS"){
	  anaphoton1->SetNCellCut(1);// At least 2 cells
	  anaphoton1->SetMinPt(0.2);
  }
  else {//EMCAL
	  //anaphoton1->SetNCellCut(0);// At least 2 cells
	  anaphoton1->SetMinPt(0.1); // no effect minium EMCAL cut.
	  anaphoton1->SetTimeCut(550,750);// Time window of [550-750] ns 
  }

  anaphoton1->SetMinDistanceToBadChannel(2, 4, 5);
  anaphoton1->SetCaloPID(pid);
  anaphoton1->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaphoton1->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaphoton1->SwitchOffDataMC() ;
  anaphoton1->SwitchOffCaloPID();
  anaphoton1->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
  if(kSimulation){
    anaphoton1->SwitchOnFiducialCut();
    anaphoton1->SetFiducialCut(fidCut1stYear);
  }
  //anaphoton1->SwitchOnTrackMatchRejection();
  if(!data.Contains("delta")) anaphoton1->SetOutputAODName(Form("PhotonsForIM%s",calorimeter.Data()));
  else                        anaphoton1->SetInputAODName (Form("PhotonsForIM%s",calorimeter.Data()));
  //Set Histograms bins and ranges
  anaphoton1->SetHistoPtRangeAndNBins(0, 50, 500) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  
  if(kPrintSettings) anaphoton1->Print("");

  AliAnaPi0 *anapi0 = new AliAnaPi0();
  anapi0->SetDebug(-1);//10 for lots of messages
  anapi0->SetInputAODName(Form("PhotonsForIM%s",calorimeter.Data()));
  anapi0->SetCaloPID(pid);
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
  anapi0->SetHistoPtRangeAndNBins(0, 50, 500) ;
  //anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //anapi0->SetHistoEtaRangeAndNBins(-0.8, 0.8, 200) ;
  anapi0->SetHistoMassRangeAndNBins(0., 0.6, 300) ;
  anapi0->SetHistoAsymmetryRangeAndNBins(0., 1. , 10) ;
  if(kPrintSettings) anapi0->Print("");
  
  
  // -------------------------------------------------
  // --- Photon Isolation and Correlation Analysis ---
  // -------------------------------------------------
  
  AliAnaPhoton *anaphoton2 = new AliAnaPhoton();
  anaphoton2->SetDebug(-1); //10 for lots of messages
  if(calorimeter == "PHOS"){
		anaphoton2->SetNCellCut(1);// At least 2 cells
		anaphoton2->SetMinPt(0.2);
  }
  else {//EMCAL
		//anaphoton2->SetNCellCut(0);// At least 2 cells
		anaphoton2->SetMinPt(0.1); // no effect minium EMCAL cut.
		anaphoton2->SetTimeCut(550,750);// Time window of [550-750] ns 
  }
  anaphoton2->SetCaloPID(pid);
  anaphoton2->SetCalorimeter(calorimeter);
  if(kUseKinematics) anaphoton2->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else  anaphoton2->SwitchOffDataMC() ;
  anaphoton2->SwitchOffCaloPID();
  anaphoton2->SwitchOffFiducialCut();
  if(kSimulation){
		anaphoton2->SwitchOnFiducialCut();
		anaphoton2->SetFiducialCut(fidCut1stYear);
  }
	
  if(!data.Contains("delta")) {
		anaphoton2->SetOutputAODName(Form("Photons%s",calorimeter.Data()));
		anaphoton2->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else anaphoton2->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  anaphoton2->AddToHistogramsName("AnaPhotonCorr_");
  //Set Histograms bins and ranges
  anaphoton2->SetHistoPtRangeAndNBins(0, 50, 500) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anaphoton2->Print("");
  // ### Isolation analysis ###	
  
  AliIsolationCut * ic = new AliIsolationCut();
  ic->SetConeSize(0.4);
  ic->SetPtThreshold(0.2);
  ic->SetICMethod(AliIsolationCut::kPtThresIC);
  if(kPrintSettings) ic->Print("");
  
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
  anaisol->SetIsolationCut(ic);	
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisol->SwitchOffReIsolation();
  //Multiple IC
  anaisol->SwitchOffSeveralIsolation() ;
  //Set Histograms bins and ranges
  anaisol->SetHistoPtRangeAndNBins(0, 50, 500) ;
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
  anacorrjet->SetHistoPtRangeAndNBins(0, 50, 500) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrjet->Print("");
  
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorrhadron = new AliAnaParticleHadronCorrelation();
  anacorrhadron->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  if(!data.Contains("delta")) {
	  anacorrhadron->SetOutputAODName(Form("CorrGammaHadrons%s",calorimeter.Data()));
	  anacorrhadron->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else anacorrhadron->SetInputAODName(Form("CorrGammaHadrons%s",calorimeter.Data()));
	
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
  anacorrhadron->SetHistoPtRangeAndNBins(0, 50, 500) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrhadron->Print("");
  
  // ### Correlation with hadrons
  AliAnaParticleHadronCorrelation *anacorrisohadron = new AliAnaParticleHadronCorrelation();
  anacorrisohadron->SetInputAODName(Form("Photons%s",calorimeter.Data()));
  if(!data.Contains("delta")) {
	  anacorrisohadron->SetOutputAODName(Form("CorrIsoGammaHadrons%s",calorimeter.Data()));
	  anacorrisohadron->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else  anacorrisohadron->SetInputAODName(Form("CorrIsoGammaHadrons%s",calorimeter.Data()));
	
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
  anacorrisohadron->SetHistoPtRangeAndNBins(0, 50, 500) ;
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
  nms->SetHistoERangeAndNBins(0, 50, 500) ;
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
  anapi0ebe->SetHistoPtRangeAndNBins(0, 50, 500) ;
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
  anaisolpi0->SetIsolationCut(ic);	
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  anaisolpi0->SwitchOffReIsolation();
  //Multiple IC
  anaisolpi0->SwitchOffSeveralIsolation() ;
  //Set Histograms bins and ranges
  anaisolpi0->SetHistoPtRangeAndNBins(0, 50, 500) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anaisol->Print("");
  
  
  // ### Pi0 Correlation with hadrons, not isolated
  AliAnaParticleHadronCorrelation *anacorrhadronpi0 = new AliAnaParticleHadronCorrelation();
  anacorrhadronpi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  if(!data.Contains("delta")){ 
	  anacorrhadronpi0->SetOutputAODName(Form("CorrPi0Hadrons%s",calorimeter.Data()));
	  anacorrhadronpi0->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else anacorrhadronpi0->SetInputAODName(Form("CorrPi0Hadrons%s",calorimeter.Data()));
	
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
  anacorrhadronpi0->SetHistoPtRangeAndNBins(0, 50, 500) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrhadronpi0->Print("");
  
  // ### Pi0 Correlation with hadrons, isolated
  AliAnaParticleHadronCorrelation *anacorrhadronisopi0 = new AliAnaParticleHadronCorrelation();
  anacorrhadronisopi0->SetInputAODName(Form("Pi0s%s",calorimeter.Data()));
  if(!data.Contains("delta")) {
		anacorrhadronisopi0->SetOutputAODName(Form("CorrIsoPi0Hadrons%s",calorimeter.Data()));
		anacorrhadronisopi0->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else  anacorrhadronisopi0->SetInputAODName(Form("CorrIsoPi0Hadrons%s",calorimeter.Data()));
	
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
  anacorrhadronisopi0->SetHistoPtRangeAndNBins(0, 50, 500) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  if(kPrintSettings) anacorrhadronisopi0->Print("");
 
  //analysis the omega->pi0+gamma
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
  anaomegaToPi0Gamma->SetHistoPtRangeAndNBins(0, 20, 200) ;
  anaomegaToPi0Gamma->SetHistoMassRangeAndNBins(0, 1, 200) ;
  anaomegaToPi0Gamma->SetPi0OverOmegaPtCut(0.8);
  anaomegaToPi0Gamma->SetGammaOverOmegaPtCut(0.2);
  if(kUseKinematics) anaomegaToPi0Gamma->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  else anaomegaToPi0Gamma->SwitchOffDataMC() ;//Access MC stack and fill more histograms
  anaomegaToPi0Gamma->AddToHistogramsName(Form("AnaOmegaToPi0Gamma%s_",calorimeter.Data()));
 if(kPrintSettings)   anaomegaToPi0Gamma->Print("");
 
  // #### Configure Maker ####
  AliAnaPartCorrMaker * maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  //if(!data.Contains("delta")) maker->AddAnalysis(qa,0);
  maker->AddAnalysis(anaphoton1,0);
  maker->AddAnalysis(anapi0,1);
  maker->AddAnalysis(anaphoton2,2);
  maker->AddAnalysis(anaisol,3);
  maker->AddAnalysis(anacorrjet,4);
  maker->AddAnalysis(anacorrhadron,5);
  maker->AddAnalysis(anacorrisohadron,6);
  maker->AddAnalysis(anapi0ebe,7);
  maker->AddAnalysis(anaisolpi0,8);
  maker->AddAnalysis(anacorrhadronpi0,9);
  maker->AddAnalysis(anacorrhadronisopi0,10);
  maker->AddAnalysis(anaomegaToPi0Gamma,11);   
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


