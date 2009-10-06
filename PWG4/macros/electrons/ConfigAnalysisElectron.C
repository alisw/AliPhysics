//------------------------------------
// Configuration macro example:
//
// Configure EMCal electron analysis
//
// Modified by: K. Read
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
  Bool_t kInputIsESD = kTRUE;     //uncomment for input ESD
//Bool_t kInputIsESD = kFALSE;    //uncomment for input AODs
  Bool_t kFollowsFilter = kTRUE;  //uncomment if follows ESD filter task
//Bool_t kFollowsFilter = kFALSE; //uncomment if no ESD filter task

  //Alternatively, select input via anaInputData environment variable.
  if (gSystem->Getenv("anaInputData")){
    TString kInputData = gSystem->Getenv("anaInputData");
    if( kInputData == "AOD" ){
      kInputIsESD = kFALSE;
      kFollowsFilter = kFALSE;
    }
  }

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
  if(kInputIsESD && !kFollowsFilter)AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
  else           AliCaloTrackAODReader *reader = new AliCaloTrackAODReader();
  reader->SetDebug(-1);//10 for lots of messages

  //Switch on or off the detectors information that you want
  reader->SwitchOnEMCAL();
  reader->SwitchOnCTS();
  //reader->SwitchOffEMCALCells();	
  reader->SwitchOffPHOS();
  //reader->SwitchOffPHOSCells();	

  //Kine
  if(!kInputIsESD){
    reader->SwitchOffStack();          // On  by default, remember to SwitchOnMCData() in analysis classes
    reader->SwitchOnAODMCParticles();  // Off by default, remember to SwitchOnMCData() in analysis classes
  }

  //Min particle pT
  reader->SetCTSPtMin(0.0);   //new
  reader->SetEMCALPtMin(0.0); //new
  if(kFollowsFilter)reader->SetTrackStatus(0);  //to prevent automatic TPC and ITS refit

  //In case of generating jet events (with PYTHIA), if pt hard bin is small
  //reject events with large difference between ptHard and triggered jet
  //reader->SetPtHardAndJetPtComparison(kTRUE);

  reader->SetFidutialCut(fidCut);

  if(!kInputIsESD){
    // Analysis with tracks, select only tracks with
    // following bits

    //     //We want tracks fitted in the detectors:
    //     ULong_t status=AliAODTrack::kTPCrefit;
    //     status|=AliAODTrack::kITSrefit;
   
    //     We want tracks whose PID bit is set:
    //     ULong_t status =AliAODTrack::kITSpid;
    //     status|=AliAODTrack::kTPCpid;	

    //	reader->SetTrackStatus(status);
  }

  reader->Print("");


  //Detector Fidutial Cuts
  AliFidutialCut * fidCut2 = new AliFidutialCut();
  fidCut2->DoEMCALFidutialCut(kTRUE) ;
  fidCut2->SetSimpleEMCALFidutialCut(0.7,80.,190.);

  fidCut2->DoCTSFidutialCut(kTRUE) ;
  fidCut2->SetSimpleCTSFidutialCut(0.9,0.,360.); 

  fidCut2->Print("");

  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------

  AliAnaElectron *anaelectron = new AliAnaElectron();
  anaelectron->SetDebug(-1); //10 for lots of messages
  anaelectron->SetCalorimeter("EMCAL");
  anaelectron->SwitchOnDataMC();
  anaelectron->SetpOverEmin(0.8);
  anaelectron->SetpOverEmax(1.2);
  anaelectron->SetResidualCut(0.02);
  anaelectron->SetMinPt(1.);
  anaelectron->SetImpactCut(1.0);  //instead of 0.5
  anaelectron->SetSdcaCut(0.05);  //instead of 0.1
  anaelectron->SetOutputAODName("Electrons");
  anaelectron->SetOutputAODClassName("AliAODPWG4Particle");
  anaelectron->SetWriteNtuple(kFALSE);
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
  printf("============================ \n");
  printf("END ConfigAnalysisElectron() \n");
  printf("============================ \n");
  return maker ;
}
