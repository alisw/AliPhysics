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
  AliFidutialCut * fidCut = new AliFidutialCut();
  fidCut->DoCTSFidutialCut(kFALSE) ;
  fidCut->DoEMCALFidutialCut(kFALSE) ;
  fidCut->DoPHOSFidutialCut(kFALSE) ;
  
  //fidCut->SetSimpleCTSFidutialCut(0.9,0.,360.);
  //fidCut->SetSimpleEMCALFidutialCut(0.7,80.,190.);
  //fidCut->SetSimplePHOSFidutialCut(0.13,220.,320.);
  
  fidCut->Print("");
  
  
  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  //AliCaloTrackAODReader *reader = new AliCaloTrackAODReader();
  AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
  reader->SetDebug(-1);//10 for lots of messages
  
  //Switch on or off the detectors information that you want
  reader->SwitchOnCTS();
  reader->SwitchOnEMCAL();
  reader->SwitchOffPHOS();
  
  //Kine
  //reader->SwitchOffStack();  // On  by default, remember to SwitchOnMCData() in analysis classes
  //reader->SwitchOnAODMCParticles(); // Off by default, remember to SwitchOnMCData() in analysis classes

  //Min particle pT
  reader->SetEMCALPtMin(1.); 
  
  //In case of generating jet events (with PYTHIA), if pt hard
  //bin is small reject events with large difference between
  //ptHard and triggered jet
  //reader->SetPtHardAndJetPtComparison(kTRUE);
  
  reader->SetFidutialCut(fidCut);
  
  //Embedding/mixing/merging analysis with other events in
  //another AOD file
  //reader->SetSecondInputFileName("./aod.root");
  //Standard event loop can have less events, start mixing at
  //event N
  //reader->SetSecondInputFirstEvent(0);
  
  // Analysis with tracks, select only tracks with
  // following bits
  
  //We want tracks fitted in the detectors:
  //     ULong_t status=AliAODTrack::kTPCrefit;
  //     status|=AliAODTrack::kITSrefit;
  
  //We want tracks whose PID bit is set:
  //     ULong_t status =AliAODTrack::kITSpid;
  //     status|=AliAODTrack::kTPCpid;
  
  //    reader->SetTrackStatus(status);
  
  reader->Print("");
  
  
  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------
  
  //Detector Fidutial Cuts
  AliFidutialCut * fidCut2 = new AliFidutialCut();
  fidCut2->DoEMCALFidutialCut(kTRUE) ;
  fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,190.);

  fidCut2->DoCTSFidutialCut(kTRUE) ;
  fidCut2->SetSimpleCTSFidutialCut(0.9,0.,360.); 

  fidCut2->Print("");

  AliAnaElectron *anaelectron = new AliAnaElectron();
  anaelectron->SetDebug(-1); //10 for lots of messages
  anaelectron->SetCalorimeter("EMCAL");
  anaelectron->SetpOverEmin(0.8);
  anaelectron->SetpOverEmax(1.2);
  anaelectron->SetResidualCut(0.05);
  anaelectron->SetWriteNtuple(kTRUE);
  anaelectron->SwitchOnDataMC();  //Access MC stack and fill more histograms
  anaelectron->SetMinPt(1.);
  anaelectron->SetOutputAODName("Electrons");
  anaelectron->SetOutputAODClassName("AliAODPWG4Particle");
  //Set Histrograms bins and ranges
  anaelectron->SetHistoPtRangeAndNBins(0, 50, 100) ;
  anaelectron->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  anaelectron->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;	
  anaelectron->Print("");
		
	
  //---------------------------------------------------------------------
  // Set  analysis algorithm and reader
  //---------------------------------------------------------------------
  maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(anaelectron,0);
  maker->SetAnaDebug(-1)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  
  maker->Print("");
  //
  printf("======================== \n");
  printf("END ConfigAnalysisElectron() \n");
  printf("======================== \n");
  return maker ;
}
