/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do photon identification analysis with ESDs
// in EMCAL, after pi0 identification event by event
// in photon paris invariant mass and aperture angle widow
// For PHOS, change EMCAL by PHOS where necessary
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
	fidCut->DoCTSFiducialCut(kFALSE) ;
	fidCut->DoEMCALFiducialCut(kFALSE) ;
	fidCut->DoPHOSFiducialCut(kTRUE) ;
	
	//fidCut->SetSimpleCTSFiducialCut(0.9,0.,360.);
	//fidCut->SetSimpleEMCALFiducialCut(0.7,80.,190.);
	fidCut->SetSimplePHOSFiducialCut(0.12,220.,320.);
	
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
	//reader->SetEMCALPtMin(0.5); 
	reader->SetPHOSPtMin(0.5);
	//reader->SetCTSPtMin(0.2);
	
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
	//anaphoton->SetFiducialCut(fidCut2); //More acceptance selections if needed at this level
	anaphoton->SetCalorimeter("EMCAL");
	anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anaphoton->SwitchOffCaloPID();
	anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
	anaphoton->SwitchOffFiducialCut();
	anaphoton->SetOutputAODName("Photons");
	anaphoton->SetOutputAODClassName("AliAODPWG4Particle");
	//Set Histrograms bins and ranges
//	anaphoton->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	anaphoton->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	anaphoton->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;	anaphoton->Print("");

	//Detector Fiducial Cuts
	AliFiducialCut * fidCut2 = new AliFiducialCut();
	fidCut2->DoPHOSFiducialCut(kTRUE) ;
	fidCut2->SetSimplePHOSFiducialCut(0.12,220.,320.);
	
	fidCut->Print("");
	
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
	anapi0->SetOutputAODClassName("AliAODPWG4Particle");
	anapi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
	anapi0->SetNeutralMesonSelection(nms);
	//Set Histrograms bins and ranges
//	anapi0->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	anapi0->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	anapi0->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
	anapi0->Print("");
	
	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->AddAnalysis(anaphoton,0);
	maker->AddAnalysis(anapi0,1);
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
