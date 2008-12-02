/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do prompt photon - hadron correlation analysis with ESDs
// First find photons with AliAnaPhoton, then
// isolate them with AliAnaParticleIsolation and finally correlate 
// them with AliAnaParticleHadronCorrelation
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
	fidCut->DoCTSFidutialCut(kTRUE) ;
	fidCut->DoEMCALFidutialCut(kTRUE) ;
	fidCut->DoPHOSFidutialCut(kTRUE) ;
	
	fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
	fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
	fidCut->SetSimplePHOSFidutialCut(0.13,220.,320.);
	
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
	reader->SwitchOffEMCALCells();
	reader->SwitchOffPHOSCells();
	
	
	//Min particle pT
	reader->SetEMCALPtMin(0.5); 
	reader->SetPHOSPtMin(0.5);
	reader->SetCTSPtMin(0.2);
	
	reader->SetFidutialCut(fidCut);
	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	//<<<< first analysis >>> select the photons
	//Detector Fidutial Cuts for analysis part
	AliFidutialCut * fidCut2 = new AliFidutialCut();
	fidCut2->DoCTSFidutialCut(kFALSE) ;
	fidCut2->DoEMCALFidutialCut(kTRUE) ;
	fidCut2->DoPHOSFidutialCut(kFALSE) ;
	
	fidCut2->SetSimpleCTSFidutialCut(0.9,0.,360.);
	fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,190.);
	fidCut2->SetSimplePHOSFidutialCut(0.13,220.,320.);
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
	anaphoton->SetFidutialCut(fidCut2); //More acceptance selections if needed at this level
	anaphoton->SetCalorimeter("PHOS");
	anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anaphoton->SwitchOnCaloPID();
	anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
	anaphoton->SwitchOnFidutialCut();
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
	
	//<<<Third analysis>>> Isolated Photon- hadron correlation
	AliAnaParticleHadronCorrelation *anacorr = new AliAnaParticleHadronCorrelation();
	anacorr->SetInputAODName("Photons");
	anacorr->SetDebug(-1);
	anacorr->SetCaloPID(pid);
	anacorr->SwitchOnCaloPID();
	anacorr->SwitchOnCaloPIDRecalculation(); //recommended for EMCAL
	anacorr->SwitchOffFidutialCut();
	anacorr->SetPtCutRange(1,100);
	anacorr->SetDeltaPhiCutRange(1.5,4.5);
	//Set Histrograms bins and ranges
	//	anaphoton->SetHistoPtRangeAndNBins(0, 50, 100) ;
	//	anaphoton->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
	//	anaphoton->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;	
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

