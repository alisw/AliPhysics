/* $Id: $ */

//------------------------------------
// Configuration macro example:
//
// Calorimeters QA
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
	
	
	
	//-----------------------------------------------------------  
	// Reader
	//-----------------------------------------------------------
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);
	
	//Switch on or off the detectors information that you want
	reader->SwitchOnEMCAL();
	reader->SwitchOffCTS();
	reader->SwitchOnPHOS();
	reader->SwitchOnEMCALCells();
	reader->SwitchOnPHOSCells();
	
	//Min particle pT
	reader->SetEMCALPtMin(0.); 
	reader->SetPHOSPtMin(0.);
	
	//     //We want tracks fitted in the detectors:
	//     ULong_t status=AliAODTrack::kTPCrefit;
	//     status|=AliAODTrack::kITSrefit; //(default settings)
	
	//     We want tracks whose PID bit is set:
	//     ULong_t status =AliAODTrack::kITSpid;
	//     status|=AliAODTrack::kTPCpid;	
	
	//	reader->SetTrackStatus(status);
	
	//Remove the temporal AODs we create.	
	reader->SwitchOffWriteStdAOD();	
	
	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	
	AliFidutialCut * fidCut = new AliFidutialCut();
	fidCut->DoCTSFidutialCut(kFALSE) ;
	fidCut->DoEMCALFidutialCut(kTRUE) ;
	fidCut->DoPHOSFidutialCut(kTRUE) ;
	
	
  	AliAnaCalorimeterQA *anaEMCAL = new AliAnaCalorimeterQA();
	anaEMCAL->SetDebug(-1); //10 for lots of messages
	anaEMCAL->SetCalorimeter("EMCAL");
	anaEMCAL->SwitchOnDataMC() ;//Access MC stack and fill more histograms
	//anaEMCAL->SetStyleMacro("$MACROS/style.C") ;
	anaEMCAL->AddToHistogramsName("AnaCaloQA_EMCAL_");
	anaEMCAL->SetFidutialCut(fidCut);
	anaEMCAL->SwitchOnFidutialCut();
	anaEMCAL->Print("");
	
	AliAnaCalorimeterQA *anaPHOS = new AliAnaCalorimeterQA();
	anaPHOS->SetDebug(-1); //10 for lots of messages
	anaPHOS->SetCalorimeter("PHOS");
	anaPHOS->SwitchOnDataMC() ;//Access MC stack and fill more histograms
	//anaPHOS->SetStyleMacro("$MACROS/style.C") ;
	anaPHOS->AddToHistogramsName("AnaCaloQA_PHOS_");
	anaPHOS->SetFidutialCut(fidCut);
	anaPHOS->SwitchOnFidutialCut();
	anaPHOS->Print("");
	
	
	//---------------------------------------------------------------------
	// Set  analysis algorithm and reader
	//---------------------------------------------------------------------
	maker = new AliAnaPartCorrMaker();
	maker->SetReader(reader);//pointer to reader
	maker->AddAnalysis(anaEMCAL,1);
	maker->AddAnalysis(anaPHOS,0);
	//maker->SetAnaDebug(0)  ;
	maker->SwitchOnHistogramsMaker()  ;
	maker->SwitchOffAODsMaker()  ;
	
	maker->Print("");
	//
	printf("======================== \n");
	printf("END ConfigAnalysis() \n");
	printf("======================== \n");
	return maker ;
}

