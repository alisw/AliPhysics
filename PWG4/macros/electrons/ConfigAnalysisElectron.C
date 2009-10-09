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

  //enum for the different electron cut sets
  //defined for dR and p/E
  //kTight2 is the default standard cuts
  enum kCutTypes {kTight2, kLooseTight, kTightLoose, kLoose2};
  Int_t kCutSet = kTight2;
  Double_t pOverEmin = 0.8;  //tight
  Double_t pOverEmax = 1.1;  //tight
  Double_t dRmax     = 0.02; //tight
  if (gSystem->Getenv("ELECUTSET")){
    kCutSet = atoi(gSystem->Getenv("ELECUTSET"));
  }
  if(kCutSet == kLooseTight) {
    pOverEmin = 0.6;  //loose
    pOverEmax = 1.3;  //loose
    dRmax     = 0.02; //tight
  }
  if(kCutSet == kTightLoose) {
    pOverEmin = 0.8;  //tight
    pOverEmax = 1.1;  //tight
    dRmax     = 0.05; //loose
  }
  if(kCutSet == kLoose2) {
    pOverEmin = 0.6;  //loose
    pOverEmax = 1.3;  //loose
    dRmax     = 0.05; //loose
  }    

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
  anaelectron->SetMinPt(1.);
  anaelectron->SetOutputAODName("Electrons");
  anaelectron->SetOutputAODClassName("AliAODPWG4Particle");
  anaelectron->SetWriteNtuple(kFALSE);
  //Determine which cuts to use based on enum
  anaelectron->SetpOverEmin(pOverEmin);
  anaelectron->SetpOverEmax(pOverEmax);
  anaelectron->SetResidualCut(dRmax);
  //Set Histrograms bins and ranges
  anaelectron->SetHistoPtRangeAndNBins(0, 100, 100) ;
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
