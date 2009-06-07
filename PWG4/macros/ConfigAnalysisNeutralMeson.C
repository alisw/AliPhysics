/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do photon identification analysis with ESDs
// in PHOS, after, pi0 identification with Invariant mass. 
// For EMCAL, change PHOS by EMCAL where necessary
//
// Author : Renzhuo Wan
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
	fidCut->DoEMCALFidutialCut(kTRUE) ;
	fidCut->DoPHOSFidutialCut(kTRUE) ;
	
	//fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
	fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
	fidCut->SetSimplePHOSFidutialCut(0.12,220.,320.);
	
	fidCut->Print("");
	
	
	//-----------------------------------------------------------  
	// Reader
	//-----------------------------------------------------------
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);//10 for lots of messages
	
	//Switch on or off the detectors information that you want
	reader->SwitchOnEMCAL();  //switch off the EMCal
	reader->SwitchOffCTS();
	reader->SwitchOnPHOS();
	
	//Min particle pT
	reader->SetEMCALPtMin(0.2); 
	reader->SetPHOSPtMin(0.2);
	//reader->SetCTSPtMin(0.2);
	
	reader->SetFidutialCut(fidCut);
	reader->Print("");
	
	//PID settings
	AliCaloPID * pid = new AliCaloPID();
	// use selection with simple weights
	pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
	// use more complicated selection, particle weight depending on cluster energy
	//	   pid->UsePHOSPIDWeightFormula(kTRUE);
	//	   TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
	//	   TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
	//	   pid->SetPHOSPhotonWeightFormula(photonF);
	//	   pid->SetPHOSPi0WeightFormula(pi0F);
	
	pid->SetDispersionCut(2);
	pid->SetTOFCut(5.e-9);
	pid->SetDebug(-1);
	pid->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm PHOS
	//---------------------------------------------------------------------
		
	AliAnaPhoton *anaphoton = new AliAnaPhoton();
	anaphoton->SetDebug(-1); //10 for lots of messages
	anaphoton->SetMinPt(0.2);
	anaphoton->SetMinDistanceToBadChannel(2, 4, 5);
	anaphoton->SetCaloPID(pid);
	anaphoton->SetCalorimeter("PHOS");
	anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anaphoton->SwitchOffCaloPID();
	anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
	anaphoton->SwitchOffFidutialCut();
	anaphoton->SetOutputAODName("PhotonsPHOS");
	anaphoton->SetOutputAODClassName("AliAODPWG4Particle");
	anaphoton->AddToHistogramsName("AnaPhotonPHOS_");
	//Set Histrograms bins and ranges
	anaphoton->SetHistoPtRangeAndNBins(0, 20, 200) ;
	anaphoton->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
	anaphoton->SetHistoEtaRangeAndNBins(-0.5, 0.5, 200) ;	
	anaphoton->Print("");
		
	AliAnaNeutralMeson *ananeutral = new AliAnaNeutralMeson();
	ananeutral->SetDebug(-1);//10 for lots of messages
	ananeutral->SetInputAODName("PhotonsPHOS");
	ananeutral->SetCaloPID(pid);
//	ananeutral->SetNCentrBin(5); //number of bins in centrality 
//	ananeutral->SetNZvertBin(5); //number of bins for vertex position
//	ananeutral->SetNRPBin(6); //number of bins in reaction plain
	ananeutral->SetAnaPi0Eta(kTRUE); //whether analysis pi0 and eta
	ananeutral->SetAnaOmega(kTRUE);   //whether analysis omega
	ananeutral->SetNPID(3);
	ananeutral->SetInvMassCut(1.);
	ananeutral->SetNEventsMixed(6);
	ananeutral->SetNAsyBinsMinMax(200,0,1.);
	ananeutral->SetNPtBinsMinMax(200,0,20.);
	ananeutral->SetNMassBinsMinMas(200,0,1.);
	ananeutral->SetPi0MassPeakWidthCut(0.015);
	ananeutral->SetCalorimeter("PHOS");
	ananeutral->SwitchOnFidutialCut();
	ananeutral->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	ananeutral->AddToHistogramsName("AnaNeuMesonPHOS_");
	ananeutral->Print("");
	
	//---------------------------------------------------------------------
	// Analysis algorithm EMCAL
	//---------------------------------------------------------------------
	
	AliAnaPhoton *anaphoton2 = new AliAnaPhoton();
	anaphoton2->SetDebug(-1); //10 for lots of messages
	anaphoton2->SetMinPt(0.2);
	anaphoton2->SetMinDistanceToBadChannel(2, 4, 5);
	anaphoton2->SetCaloPID(pid);
	anaphoton2->SetCalorimeter("EMCAL");
	anaphoton2->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anaphoton2->SwitchOffCaloPID();
	anaphoton2->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
	anaphoton2->SwitchOffFidutialCut();
	anaphoton2->SetOutputAODName("PhotonsEMCAL");
	anaphoton2->SetOutputAODClassName("AliAODPWG4Particle");
	anaphoton2->AddToHistogramsName("AnaPhotonEMCAL_");
	//Set Histrograms bins and ranges
	anaphoton2->SetHistoPtRangeAndNBins(0, 20, 200) ;
	anaphoton2->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
	anaphoton2->SetHistoEtaRangeAndNBins(-0.5, 0.5, 200) ;	
	anaphoton2->Print("");
	
	AliAnaNeutralMeson *ananeutral2 = new AliAnaNeutralMeson();
	ananeutral2->SetDebug(-1);//10 for lots of messages
	ananeutral2->SetInputAODName("PhotonsEMCAL");
	ananeutral2->SetCaloPID(pid);
	//	ananeutral2->SetNCentrBin(5); //number of bins in centrality 
	//	ananeutral2->SetNZvertBin(5); //number of bins for vertex position
	//	ananeutral2->SetNRPBin(6); //number of bins in reaction plain
	ananeutral2->SetAnaPi0Eta(kTRUE); //whether analysis pi0 and eta
	ananeutral2->SetAnaOmega(kTRUE);   //whether analysis omega
	ananeutral2->SetNPID(3);
	ananeutral2->SetInvMassCut(1.);
	ananeutral2->SetNEventsMixed(6);
	ananeutral2->SetNAsyBinsMinMax(200,0,1.);
	ananeutral2->SetNPtBinsMinMax(200,0,20.);
	ananeutral2->SetNMassBinsMinMas(200,0,1.);
	ananeutral2->SetPi0MassPeakWidthCut(0.015);
	ananeutral2->SetCalorimeter("EMCAL");
	ananeutral2->SwitchOnFidutialCut();
	ananeutral2->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	ananeutral2->AddToHistogramsName("AnaNeutralMesonEMCAL_");
	ananeutral2->Print("");
	
	
	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->AddAnalysis(anaphoton,0);
	maker->AddAnalysis(ananeutral,1);
	maker->AddAnalysis(anaphoton2,2);
	maker->AddAnalysis(ananeutral2,3);
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
