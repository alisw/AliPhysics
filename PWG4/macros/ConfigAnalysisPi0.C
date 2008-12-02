/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do photon identification analysis with ESDs
// in PHOS, after, pi0 identification with Invariant mass. 
// For EMCAL, change PHOS by EMCAL where necessary
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
	fidCut->SetSimplePHOSFidutialCut(0.12,220.,320.);
	
	fidCut->Print("");
	
	
	//-----------------------------------------------------------  
	// Reader
	//-----------------------------------------------------------
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);//10 for lots of messages
	
	//Switch on or off the detectors information that you want
	reader->SwitchOnEMCAL();
	reader->SwitchOnCTS();
	reader->SwitchOnPHOS();
	
	//Min particle pT
	//reader->SetEMCALPtMin(0.5); 
	reader->SetPHOSPtMin(0.5);
	//reader->SetCTSPtMin(0.2);
	
	reader->SetFidutialCut(fidCut);
	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	
	AliCaloPID * pid = new AliCaloPID();
	// use selection with simple weights
	pid->SetPHOSPhotonWeight(0.7);    pid->SetPHOSPi0Weight(0.7); 
	// use more complicated selection, particle weight depending on cluster energy
	//   pid->UsePHOSPIDWeightFormula(kTRUE);
	//   TFormula * photonF = new TFormula("photonWeight","0.98*(x<40)+ 0.68*(x>=100)+(x>=40 && x<100)*(0.98+x*(6e-3)-x*x*(2e-04)+x*x*x*(1.1e-06))");
	//   TFormula * pi0F = new TFormula("pi0Weight","0.98*(x<65)+ 0.915*(x>=100)+(x>=65 && x-x*(1.95e-3)-x*x*(4.31e-05)+x*x*x*(3.61e-07))");
	//   pid->SetPHOSPhotonWeightFormula(photonF);
	//   pid->SetPHOSPi0WeightFormula(pi0F);
	
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
	anaphoton->SetCalorimeter("PHOS");
	anaphoton->SwitchOffDataMC() ;//Access MC stack and fill more histograms
	anaphoton->SwitchOffCaloPID();
	anaphoton->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL
	anaphoton->SwitchOffFidutialCut();
	anaphoton->SetOutputAODName("Photons");
	anaphoton->SetOutputAODClassName("AliAODPWG4Particle");
	//Set Histrograms bins and ranges
//	anaphoton->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	anaphoton->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	anaphoton->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;	anaphoton->Print("");
	
	//Detector Fidutial Cuts
	AliFidutialCut * fidCut2 = new AliFidutialCut();
	fidCut2->DoPHOSFidutialCut(kTRUE) ;
	fidCut2->SetSimplePHOSFidutialCut(0.12,220.,320.);
	
	fidCut->Print("");
		
	AliAnaPi0 *anapi0 = new AliAnaPi0();
	anapi0->SetDebug(-1);//10 for lots of messages
        anapi0->SetInputAODName("Photons");
	anapi0->SetCaloPID(pid);
//	anapi0->SetNCentrBin(5); //number of bins in centrality 
//	anapi0->SetNZvertBin(5); //number of bins for vertex position
//	anapi0->SetNRPBin(6); //number of bins in reaction plain
//	anapi0->SetNMaxEvMix(20);//Maximal number of events for mixing
	anapi0->SetCalorimeter("PHOS");
	anapi0->SetFidutialCut(fidCut2); //More acceptance selections if needed at this level. Not used if real geometry
	anapi0->SwitchOnFidutialCut();
	anapi0->SwitchOnDataMC() ;//Access MC stack and fill more histograms
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
