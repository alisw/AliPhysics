//------------------------------------
// Config file for B-tagging
//
// Author T. Aronsson
//
//------------------------------------

AliAnaPartCorrMaker*  ConfigAnalysis()
{
  //
  // Configuration goes here
  // 
  printf("======================== \n");
  printf("==Preforming Btag.cxx=== \n");
  printf("======================== \n");
  Bool_t kInputIsESD = kTRUE;     //uncomment for input ESD
//Bool_t kInputIsESD = kFALSE;    //uncomment for input AODs
  Bool_t kFollowsFilter = kTRUE;  //uncomment if follows ESD filter task
//Bool_t kFollowsFilter = kFALSE; //uncomment if no ESD filter task
  Bool_t kMC = kTRUE; //set to kFALSE for data

   

  //Alternatively, select input via anaInputData environment variable.
  if (gSystem->Getenv("anaInputData")){
    TString kInputData = gSystem->Getenv("anaInputData");
    if( kInputData == "AOD" ){
      kInputIsESD = kFALSE;
      kFollowsFilter = kFALSE;
    }
  }

  //Alternatively, adjust for real data based on kMC value.
  if (gSystem->Getenv("anakMC")){
    kMC = atoi(gSystem->Getenv("anakMC"));
  }

  //Detector Fiducial Cuts
  AliFiducialCut * fidCut = new AliFiducialCut();
  fidCut->DoCTSFiducialCut(kFALSE) ;
  fidCut->DoEMCALFiducialCut(kFALSE) ;
  fidCut->DoPHOSFiducialCut(kFALSE) ;

  //fidCut->SetSimpleCTSFiducialCut(0.9,0.,360.);
  //fidCut->SetSimpleEMCALFiducialCut(0.7,80.,190.);
  //fidCut->SetSimplePHOSFiducialCut(0.13,220.,320.);

  fidCut->Print("");

  //-----------------------------------------------------------  
  // Reader
  //-----------------------------------------------------------
  if(kInputIsESD && !kFollowsFilter)AliCaloTrackESDReader *reader = new AliCaloTrackESDReader();
  else           AliCaloTrackAODReader *reader = new AliCaloTrackAODReader();
  reader->SetDebug(-1);//10 for lots of messages
  reader->SwitchOnWriteDeltaAOD();
  //Switch on or off the detectors information that you want
  reader->SwitchOnEMCAL();
  reader->SwitchOnCTS();	
  reader->SwitchOffPHOS();
	

  //Kine
  if(kMC && !kInputIsESD){
    reader->SwitchOffStack();          // On  by default, remember to SwitchOnMCData() in analysis classes
    reader->SwitchOnAODMCParticles();  // Off by default, remember to SwitchOnMCData() in analysis classes
  }

  //Select Trigger Class for real data
  if(!kMC) reader->SetFiredTriggerClassName("CINT1B-ABCE-NOPF-ALL");

  //Min particle pT
  reader->SetCTSPtMin(0.0);   //new
  reader->SetEMCALPtMin(0.0); //new
  if(kFollowsFilter)reader->SetTrackStatus(0);  //to prevent automatic TPC and ITS refit



  reader->SetFiducialCut(fidCut);

  if(!kInputIsESD){

  }

  reader->Print("");


  //Detector Fiducial Cuts
  AliFiducialCut * fidCut2 = new AliFiducialCut();
  fidCut2->DoEMCALFiducialCut(kTRUE) ;
  fidCut2->SetSimpleEMCALFiducialCut(0.7,80.,190.);

  fidCut2->DoCTSFiducialCut(kTRUE) ;
  fidCut2->SetSimpleCTSFiducialCut(0.9,0.,360.); 

  fidCut2->Print("");





  //---------------------------------------------------------------------
  // Analysis algorithm
  //---------------------------------------------------------------------

  AliAnaBtag *btag = new AliAnaBtag();

 //Base class
  btag->SetDebug(-1); //10 for lots of messages
  //btag->SetWriteNtuple(1); //Can be used to write out NTuples for local analysis (1000 times faster than AOD analysis), default is off.
  if(kMC){
    btag->SwitchOnDataMC();
    btag->SetMinPt(1.);
  }
  btag->SetOutputAODName("Electrons");
  btag->SetOutputAODClassName("AliAODPWG4Particle");





  //---------------------------------------------------------------------
  // Set  analysis algorithm and reader
  //---------------------------------------------------------------------
  maker = new AliAnaPartCorrMaker();
  maker->SetReader(reader);//pointer to reader
  maker->AddAnalysis(btag,0);
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
