/* $Id: $ */

//------------------------------------
// Configuration macro example:
//
// Calorimeters QA: Validation of data (MC)
// Valid for ESDs, comment/uncomment below  in the reader part for AODs
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
	//For this particular analysis few things done by the reader.
	//Nothing else needs to be set.
	AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
	reader->SetDebug(-1);
	reader->SwitchOnStack();          
	reader->SwitchOffAODMCParticles(); 	
	reader->SetDeltaAODFileName(""); //Do not create deltaAOD file, this analysis do not create branches.
	reader->Print("");
	
	//For AODs:
//	AliCaloTrackAODReader *reader = new AliCaloTrackAODReader();
//	reader->SetDebug(-1);
//	reader->SwitchOffStack();          
//	reader->SwitchOnAODMCParticles(); 	
//	reader->SetDeltaAODFileName(""); //Do not create deltaAOD file, this analysis do not create branches.
//	reader->Print("");
	
	
	//---------------------------------------------------------------------
	// Analysis algorithm
	//---------------------------------------------------------------------
	
	AliFiducialCut * fidCut = new AliFiducialCut();
	fidCut->DoCTSFiducialCut(kFALSE) ;
	fidCut->DoEMCALFiducialCut(kTRUE) ;
	fidCut->DoPHOSFiducialCut(kTRUE) ;
	
	
  	AliAnaCalorimeterQA *anaEMCAL = new AliAnaCalorimeterQA();
	anaEMCAL->SetDebug(-1); //10 for lots of messages
	anaEMCAL->SetCalorimeter("EMCAL");
	anaEMCAL->SwitchOnDataMC() ;//Access MC stack and fill more histograms
	//anaEMCAL->SetStyleMacro("$MACROS/style.C") ;
	anaEMCAL->AddToHistogramsName("AnaCaloQA_EMCAL_");
	anaEMCAL->SetFiducialCut(fidCut);
	anaEMCAL->SwitchOnFiducialCut();
	anaEMCAL->Print("");
	
	AliAnaCalorimeterQA *anaPHOS = new AliAnaCalorimeterQA();
	anaPHOS->SetDebug(-1); //10 for lots of messages
	anaPHOS->SetCalorimeter("PHOS");
	anaPHOS->SwitchOnDataMC() ;//Access MC stack and fill more histograms
	//anaPHOS->SetStyleMacro("$MACROS/style.C") ;
	anaPHOS->AddToHistogramsName("AnaCaloQA_PHOS_");
	anaPHOS->SetFiducialCut(fidCut);
	anaPHOS->SwitchOnFiducialCut();
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
	maker->SwitchOffAODsMaker()  ;//No AODs created in this task.
	
	maker->Print("");
	//
	printf("======================== \n");
	printf("END ConfigAnalysis() \n");
	printf("======================== \n");
	return maker ;
}

