/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do some track selection, for input of correlation analysis.
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
	fidCut->DoEMCALFidutialCut(kFALSE) ;
	fidCut->DoPHOSFidutialCut(kFALSE) ;
	
	fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
	fidCut->Print("");
	
	
	//-----------------------------------------------------------  
	// Reader
	//-----------------------------------------------------------
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);//10 for lots of messages
	
	//Switch on or off the detectors information that you want
	reader->SwitchOffEMCAL();
	reader->SwitchOffPHOS();
	reader->SwitchOnCTS(); //We want only tracks
	
	//Min particle pT
	reader->SetCTSPtMin(0.2);
	
	//Uncomment only with MCReader
//	TArrayI statusArray(1) ;
//	statusArray.SetAt(1,0); 
//	reader->AddStatusArray(statusArray)  ;
//	reader->SwitchOnStatusSelection() ;
	
	reader->SetFidutialCut(fidCut);

	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------

	AliAnaChargedParticles *anatrack = new AliAnaChargedParticles();
	anatrack->SetDebug(-1);//10 for lots of messages
	anatrack->SetMinPt(5.);
	anatrack->SetOutputAODName("Charged");
	anatrack->SwitchOffFidutialCut();
	anatrack->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
	anatrack->SwitchOnDataMC() ;//Access MC stack and fill more histograms
	//Set Histrograms bins and ranges
//	anatrack->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	anatrack->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	anatrack->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
	anatrack->Print("");
	
	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->AddAnalysis(anatrack,0);
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
