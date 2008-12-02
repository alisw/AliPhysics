/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do photon identification analysis with ESDs
// in EMCAL, after pi0 identification event by event
// in photon paris invariant mass and aperture angle widow.
// After, study the isolation of the selected pi0.
// For PHOS, change EMCAL by PHOS where necessary
// Do photon identification analysis with ESDs
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
	
	
	//Detector Fidutial Cuts
	AliFidutialCut * fidCut = new AliFidutialCut();
	fidCut->DoCTSFidutialCut(kFALSE) ;
	fidCut->DoEMCALFidutialCut(kFALSE) ;
	fidCut->DoPHOSFidutialCut(kTRUE) ;
	
	//fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
	//fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
	//fidCut->SetSimplePHOSFidutialCut(0.12,220.,320.);
	
	fidCut->Print("");
	
	//-----------------------------------------------------------  
	// Reader
	//-----------------------------------------------------------
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);//10 for lots of messages
	
	//Switch on or off the detectors information that you want
	reader->SwitchOnEMCAL();
	reader->SwitchOnCTS();
	reader->SwitchOffPHOS();
	
	//Min particle pT
	reader->SetEMCALPtMin(0.5); 
	//reader->SetPHOSPtMin(0.5);
	reader->SetCTSPtMin(0.2);
	
	reader->SetFidutialCut(fidCut);
	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	
	AliCaloPID * pid = new AliCaloPID();
	// use selection with simple weights
	pid->SetEMCALPhotonWeight(0.7);    pid->SetEMCALPi0Weight(0.7);
	pid->SetDispersionCut(1.5);
	pid->SetTOFCut(5.e-9);
	pid->SetDebug(-1);
	pid->Print("");

	AliAnaPhoton *anaphoton = new AliAnaPhoton();
	anaphoton->SetDebug(-1); //10 for lots of messages
	anaphoton->SetMinPt(1);
	anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
	anaphoton->SetCaloPID(pid);
	//anaphoton->SetFidutialCut(fidCut2); //More acceptance selections if needed at this level
	anaphoton->SetCalorimeter("EMCAL");
	anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anaphoton->SwitchOffCaloPID();
	anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
	anaphoton->SwitchOffFidutialCut();
	anaphoton->SetOutputAODName("Photons");
	anaphoton->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
	//Set Histrograms bins and ranges
//	anaphoton->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	anaphoton->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	anaphoton->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
	anaphoton->Print("");

	// >>> Second Analysis <<<
	
	AliNeutralMesonSelection *nms = new AliNeutralMesonSelection();
	nms->SetInvMassCutRange(0.10, 0.17)	;
	nms->KeepNeutralMesonSelectionHistos(kTRUE);
	//Set Histrograms bins and ranges
//	nms->SetHistoERangeAndNBins(0, 50, 100) ;
//	nms->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	nms->SetHistoAngleRangeAndNBins(0, 0.3, 100) ;
//	nsm->SetHistoIMRangeAndNBins(0, 0.4, 100) ;  
	
	AliAnaPi0EbE *anapi0 = new AliAnaPi0EbE();
	anapi0->SetDebug(-1);//10 for lots of messages
	anapi0->SetAnalysisType(AliAnaPi0EbE::kIMCalo);
	anapi0->SetInputAODName("Photons");
	anapi0->SetOutputAODName("Pi0s");
	anapi0->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
	anapi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
	anapi0->SetNeutralMesonSelection(nms);
	//Set Histrograms bins and ranges
//	anapi0->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	anapi0->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
	anapi0->Print("");
	
	// >>>> Third Analysis <<<<
	AliIsolationCut * ic = new AliIsolationCut();
	ic->SetConeSize(0.4);
	ic->SetPtThreshold(1.);
	ic->SetICMethod(AliIsolationCut::kPtThresIC);
	ic->Print("");
	
	AliAnaParticleIsolation *anaisol = new AliAnaParticleIsolation();
	anaisol->SetDebug(-1);
	anaisol->SetInputAODName("Pi0s");
	anaisol->SetMinPt(5);
	anaisol->SetCalorimeter("EMCAL"); //Must be the same calorimeter where pi0/photon falls
	anaisol->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	//Select clusters with no pair, if both clusters with pi0 mass
	anaisol->SwitchOffInvariantMass(); //It is a pi0 not needed
	//anaisol->SetNeutralMesonSelection(nms);
	//Do isolation cut
	anaisol->SetIsolationCut(ic);	
	//Do or not do isolation with previously produced AODs.
	//No effect if use of SwitchOnSeveralIsolation()
	anaisol->SwitchOffReIsolation();
	
	//Multiple IC
	anaisol->SwitchOffSeveralIsolation() ;
	//Make sure that the initial isolation cut parameters set in AliIsolationCut 
	//pointer are relativelly open.
//	anaisol->SetNCones(2) ;     
//	anaisol->SetNPtThresFrac(2) ;     
//	anaisol->SetConeSizes(0, 0.3) ;		anaisol->SetConeSizes(1, 0.4) ;    
//	anaisol->SetPtThresholds(0, 0.5) ;	anaisol->SetPtThresholds(1, 1.) ;   
//	anaisol->SetPtFractions(0, 1.) ;	anaisol->SetPtFractions(1, 1.5) ;  
	//Set Histrograms bins and ranges
//	anaisol->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	anaisol->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	anaisol->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;	
	anaisol->SetHistoPtSumRangeAndNBins(0, 10, 100) ;
	anaisol->SetHistoPtInConeRangeAndNBins(0, 50, 100) ;	
	anaisol->Print("");
	
	
	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->AddAnalysis(anaphoton,0);
	maker->AddAnalysis(anapi0,1);
	maker->AddAnalysis(anaisol,2);
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
