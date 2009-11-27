/* $Id: $ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon - jet correlation analysis with ESDs
// First find photons with AliAnaPhoton, then
// isolate them with AliAnaParticleIsolation and finally correlate           
// them with AliAnaParticleJetLeadingConeCorrelation
//
// Author : Gustavo Conesa Balbastre (INFN-LNF)
//------------------------------------

AliAnaPartCorrMaker*  ConfigAnalysis()
{
	//
	// Configuration goes here
	// 
	printf("======================== \n");
	printf("ConfigAnalysis() \n");
	printf("======================== \n");
	
	
	//Detector Fiducial Cuts
	AliFiducialCut * fidCut = new AliFiducialCut();
	fidCut->DoCTSFiducialCut(kTRUE) ;
	fidCut->DoEMCALFiducialCut(kTRUE) ;
	fidCut->DoPHOSFiducialCut(kTRUE) ;
	
	fidCut->SetSimpleCTSFiducialCut(0.9,0.,360.);
	fidCut->SetSimpleEMCALFiducialCut(0.7,80.,190.);
	fidCut->SetSimplePHOSFiducialCut(0.13,220.,320.);
	
	fidCut->Print("");
	
	//-----------------------------------------------------------  
	// Reader
	//-----------------------------------------------------------
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);
	
	//Switch on or off the detectors information that you want
	reader->SwitchOnEMCAL();
	reader->SwitchOnCTS();
	reader->SwitchOnPHOS();
	
	//Min particle pT
	reader->SetEMCALPtMin(0.5); 
	reader->SetPHOSPtMin(0.5);
	reader->SetCTSPtMin(0.2);
	
	reader->SetFiducialCut(fidCut);
	
	//     //We want tracks fitted in the detectors:
	//     ULong_t status=AliAODTrack::kTPCrefit;
	//     status|=AliAODTrack::kITSrefit; //(default settings)
	
	//     We want tracks whose PID bit is set:
	//     ULong_t status =AliAODTrack::kITSpid;
	//     status|=AliAODTrack::kTPCpid;	
	
	//	reader->SetTrackStatus(status);
	
	//Remove the temporal AODs we create.	
	reader->SwitchOnCleanStdAOD();	
	
	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	//<<<< first analysis >>> select the photons
	//Detector Fiducial Cuts for analysis part
	AliFiducialCut * fidCut2 = new AliFiducialCut();
	fidCut2->DoCTSFiducialCut(kFALSE) ;
	fidCut2->DoEMCALFiducialCut(kTRUE) ;
	fidCut2->DoPHOSFiducialCut(kFALSE) ;
	
	fidCut2->SetSimpleCTSFiducialCut(0.9,0.,360.);
	fidCut2->SetSimpleEMCALFiducialCut(0.7,80.,190.);
	fidCut2->SetSimplePHOSFiducialCut(0.13,220.,320.);
	fidCut2->Print("");
	
	AliCaloPID * pid = new AliCaloPID();
	// use selection with simple weights
	pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
	pid->SetEMCALPhotonWeight(0.7);    pid->SetEMCALPi0Weight(0.7);
	
	pid->Print("");
	
  	AliAnaPhoton *anaphoton = new AliAnaPhoton();
	anaphoton->SetDebug(-1); //10 for lots of messages
	anaphoton->SetMinPt(0.2);
	anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
	anaphoton->SetCaloPID(pid);
	anaphoton->SetFiducialCut(fidCut2); //More acceptance selections if needed at this level
	anaphoton->SetCalorimeter("PHOS");
	anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anaphoton->SwitchOnCaloPID();
	anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
	anaphoton->SwitchOnFiducialCut();
	anaphoton->SetOutputAODName("Photons");
	anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
	//Set Histrograms bins and ranges
	//	anaphoton->SetHistoPtRangeAndNBins(0, 50, 100) ;
	//	anaphoton->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
	//	anaphoton->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
	anaphoton->Print("");
	
	// >>>> Second Analysis <<<< Isolate the photons
	AliIsolationCut * ic = new AliIsolationCut();
	ic->SetConeSize(0.5);
	ic->SetPtThreshold(1.);
	ic->SetICMethod(AliIsolationCut::kPtThresIC);
	ic->Print("");
	
	AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
	anaisol->SetDebug(-1);
	anaisol->SetMinPt(5);
	anaisol->SetInputAODName("Photons");
	anaisol->SetCalorimeter("PHOS");
	anaisol->SwitchOffDataMC() ;//Access MC stack and fill more histograms
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
	
	anaisol->Print("");
	
	//<<<Third analysis>>> Isolated Photon- Jet Leading in opposite Cone correlation
	AliAnaParticleJetLeadingConeCorrelation *anacorr = new AliAnaParticleJetLeadingConeCorrelation();
	anacorr->SetDebug(-1);
	anacorr->SetInputAODName("Photons");
	anacorr->SwitchOnCaloPID();
	anacorr->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
	anacorr->SwitchOffFiducialCut();
	anacorr->SwitchOffJetsOnlyInCTS();
	anacorr->SwitchOffJetsRecalculation();
	//Analysis cuts for leading particle selection
	anacorr->SetDeltaPhiCutRange(1.5,4.5); //Back-Leading particle angular cut
	anacorr->SetLeadingRatioCutRange(0.,3);//Cut for the momentum of leading
	//Analysis cuts for jet selection
	anacorr->SetppCollisions(); //Jet particles Pt threshold for different collisions
	anacorr->SetCone(0.7); //Jet cone size
	anacorr->SetJetPtThreshold(0.2); //Jet particle threshold 
	anacorr->SetJetRatioCutRange(0.7, 1.3);//Only if SwitchOffJetsOnlyInCTS(); and SetJetSelectionMode(2)
	anacorr->SetJetCTSRatioCutRange(0.3,1.3); //Only if SwitchOnJetsOnlyInCTS(); and SetJetSelectionMode(2)

	anacorr->Print("");	
	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->AddAnalysis(anaphoton,0);
	maker->AddAnalysis(anaisol,1);
	maker->AddAnalysis(anacorr,2);
	maker->SetAnaDebug(-1)  ;
	maker->SwitchOnHistogramsMaker()  ;
	maker->SwitchOnAODsMaker()  ;
	
	maker->Print("");
	//
	printf("======================== \n");
	printf("END ConfigAnalysis() \n");
	printf("======================== \n");
	return maker ;
}

