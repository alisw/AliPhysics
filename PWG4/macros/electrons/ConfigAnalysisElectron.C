/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do EMCal electron analysis with ESDs
//
//------------------------------------

AliAnaPartCorrMaker*  ConfigAnalysis()
{
	//
	// Configuration goes here
	// 
	printf("======================== \n");
	printf("ConfigAnalysisElectron() \n");
	printf("======================== \n");
	
	
	//Detector Fidutial Cuts
	//AliFidutialCut * fidCut = new AliFidutialCut();
	//fidCut->DoEMCALFidutialCut(kTRUE) ;
	
	//fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
	
	//fidCut->Print("");
	
	
	//-----------------------------------------------------------  
	// Reader
	//-----------------------------------------------------------
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(10);//10 for lots of messages
	
	//Switch on or off the detectors information that you want
	reader->SwitchOffCTS();
	reader->SwitchOffEMCAL();
	reader->SwitchOffEMCALCells();	
	reader->SwitchOffPHOS();
	reader->SwitchOffPHOSCells();	
	//Min particle pT
	//reader->SetEMCALPtMin(1.); 
	
	//reader->SetFidutialCut(fidCut);
	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	
	AliAnaElectron *anaelectron = new AliAnaElectron();
	anaelectron->SetDebug(10); //10 for lots of messages
	anaelectron->SetCalorimeter("EMCAL");
	anaelectron->SetpOverEmin(0.8);
	anaelectron->SetpOverEmax(1.2);
	anaelectron->SetResidualCut(0.05);
	anaelectron->SwitchOnDataMC();
	anaelectron->SetMinPt(1.);
	//Set Histrograms bins and ranges
	anaelectron->SetHistoPtRangeAndNBins(0, 50, 100) ;
	anaelectron->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
	anaelectron->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;	anaelectron->Print("");
	
	//Detector Fidutial Cuts
	AliFidutialCut * fidCut2 = new AliFidutialCut();
	fidCut2->DoEMCALFidutialCut(kTRUE) ;
	fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,190.);
	fidCut2->Print("");
		
	
	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->AddAnalysis(anaelectron,0);
	maker->SetAnaDebug(10)  ;
	maker->SwitchOnHistogramsMaker()  ;
	//maker->SwitchOnAODsMaker()  ;
	
	maker->Print("");
	//
	printf("======================== \n");
	printf("END ConfigAnalysisElectron() \n");
	printf("======================== \n");
	return maker ;
}
