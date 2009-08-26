/* $Id: $ */
/* $Log$ */

//------------------------------------
// Configuration macro example:
//
// Do photon identification analysis with Generator data
// Gamma in PHOS. For EMCAL clusters change 
// PHOS by EMCAL where necessary
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
  
  //fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
  //fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  //fidCut->SetSimplePHOSFidutialCut(0.13,220.,320.);
   
  fidCut->Print("");
  
  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  AliCaloTrackMCReader *reader = new AliCaloTrackMCReader();
  reader->SetDebug(-1);

  //Switch on or off the detectors information that you want
  reader->SwitchOffEMCAL();
  reader->SwitchOffCTS();
  reader->SwitchOnPHOS();
  
  //Min particle pT
  //reader->SetEMCALPtMin(0.5); 
  reader->SetPHOSPtMin(0.5);
  //reader->SetCTSPtMin(0.2);
  reader->SwitchOnStack(); //On by default, remember to SwitchOnMCData() in analysis classes
  //In case of generating jet events (with PYTHIA), if pt hard bin is small
  //reject events with large difference between ptHard and triggered jet	
  //reader->SetPtHardAndJetPtComparison(kTRUE);
	
  reader->SetFidutialCut(fidCut);

  //Anaysis of final particles, not pi0/eta etc.
  TArrayI statusArray(1) ;
  statusArray.SetAt(1,0); 
  reader->AddStatusArray(statusArray)  ;
  reader->SwitchOnStatusSelection() ;

  //Keep pi0 in the list and not the 2 photon if decay angle is small.
  reader->SwitchOffOverlapCheck();        	
	
  reader->Print("");
  
  
  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------
  
  //Detector Fidutial Cuts for analysis part
  AliFidutialCut * fidCut2 = new AliFidutialCut();
  fidCut2->DoCTSFidutialCut(kFALSE) ;
  fidCut2->DoEMCALFidutialCut(kFALSE) ;
  fidCut2->DoPHOSFidutialCut(kFALSE) ;
  
  //fidCut2->SetSimpleCTSFidutialCut(0.9,0.,360.);
  //fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  //fidCut2->SetSimplePHOSFidutialCut(0.13,220.,320.);

  fidCut2->Print("");

  AliAnaPhoton *ana = new AliAnaPhoton();
  ana->SetDebug(-1);
  ana->SetMinPt(5.);
  ana->SetMinDistanceToBadChannel(2, 4, 5);
  ana->SetFidutialCut(fidCut2);
  ana->SetCalorimeter("PHOS");
  ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  ana->SwitchOffCaloPID(); //No need with MC reader
  ana->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL, no need with MC reader
  ana->SwitchOffTrackMatchRejection(); //Only in use when OnCaloPID
  ana->SwitchOffFidutialCut();
  ana->SetOutputAODName("Photons");
  ana->SetOutputAODClassName("AliAODPWG4Particle");
  //Set Histrograms bins and ranges
//	ana->SetHistoPtRangeAndNBins(0, 50, 100) ;
//	ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
//	ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;
  ana->Print("");

  //---------------------------------------------------------------------
  // Set  analysis algorithm and reader
  //---------------------------------------------------------------------
  maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(ana,0);
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
