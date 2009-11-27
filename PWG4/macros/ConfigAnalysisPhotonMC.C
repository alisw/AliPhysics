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
  
  
  //Detector Fiducial Cuts
  AliFiducialCut * fidCut = new AliFiducialCut();
  fidCut->DoCTSFiducialCut(kTRUE) ;
  fidCut->DoEMCALFiducialCut(kTRUE) ;
  fidCut->DoPHOSFiducialCut(kTRUE) ;
  
  //fidCut->SetSimpleCTSFiducialCut(0.9,0.,360.);
  //fidCut->SetSimpleEMCALFiducialCut(0.7,80.,190.);
  //fidCut->SetSimplePHOSFiducialCut(0.13,220.,320.);
   
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
	
  reader->SetFiducialCut(fidCut);

  //Anaysis of final particles, not pi0/eta etc.
  TArrayI statusArray(1) ;
  statusArray.SetAt(1,0); 
  reader->AddStatusArray(statusArray)  ;
  reader->SwitchOnStatusSelection() ;

  //Keep pi0 in the list and not the 2 photon if decay angle is small.
  reader->SwitchOffOverlapCheck();        	

  //Remove the temporal AODs we create.	
  reader->SwitchOnCleanStdAOD();	
	
  reader->Print("");
  
  
  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------
  
  //Detector Fiducial Cuts for analysis part
  AliFiducialCut * fidCut2 = new AliFiducialCut();
  fidCut2->DoCTSFiducialCut(kFALSE) ;
  fidCut2->DoEMCALFiducialCut(kFALSE) ;
  fidCut2->DoPHOSFiducialCut(kFALSE) ;
  
  //fidCut2->SetSimpleCTSFiducialCut(0.9,0.,360.);
  //fidCut2->SetSimpleEMCALFiducialCut(0.7,80.,190.);
  //fidCut2->SetSimplePHOSFiducialCut(0.13,220.,320.);

  fidCut2->Print("");

  AliAnaPhoton *ana = new AliAnaPhoton();
  ana->SetDebug(-1);
  ana->SetMinPt(5.);
  ana->SetMinDistanceToBadChannel(2, 4, 5);
  ana->SetFiducialCut(fidCut2);
  ana->SetCalorimeter("PHOS");
  ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  ana->SwitchOffCaloPID(); //No need with MC reader
  ana->SwitchOffCaloPIDRecalculation(); //recommended for EMCAL, no need with MC reader
  ana->SwitchOffTrackMatchRejection(); //Only in use when OnCaloPID
  ana->SwitchOffFiducialCut();
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
