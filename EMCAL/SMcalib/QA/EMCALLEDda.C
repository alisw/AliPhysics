/*
  a very basic test macro for processing LED events with AliCaloCalibSignal
  Josh Hamblen
*/
#define AliDebugLevel() -1

void EMCALLEDda(const int runno = 476){
  
  int i, status;
  
  /* log start of process */
  printf("EMCAL DA started - %s\n",__FILE__);
  
  AliCaloCalibSignal * calibSignal = new 
    AliCaloCalibSignal(AliCaloCalibSignal::kEmCal); // pedestal and noise calibration
  
  int nevents=0;
  // Assume we are just interested in the 1st segment, _0.root below for fname*
  Char_t fname[256];
  sprintf(fname, "/local/data/Run_%09d.Seq_1A.Stream_0.root",runno);
    
  AliRawReader *rawReader = NULL;
  rawReader = new AliRawReaderRoot(fname);
  AliCaloRawStream *in = NULL; 
  in = new AliCaloRawStream(rawReader,"EMCAL");
  //in->SetOldRCUFormat(kTRUE);
  
  AliRawEventHeaderBase *aliHeader=NULL;
  int nev=0;
  /* read until EOF */
  while ( rawReader->NextEvent()) {
    nev++;
    
    aliHeader = (AliRawEventHeaderBase*) rawReader->GetEventHeader();
    calibSignal->SetRunNumber( aliHeader->Get("RunNb") ); // just for fun; keep info on last run looked at
    
    if ( aliHeader->Get("Type") == AliRawEventHeaderBase::kStartOfRun){
      // set start time
      calibSignal->SetStartTime(aliHeader->Get("Timestamp"));
      int time = aliHeader->Get("Timestamp");
    }
    
    // set parameters for led event checking
    calibSignal->SetAmpCut(50);
    // for now, fraction choice is set for testbeam files:
    // 64 channels / (a total of 24 rows * 48 cols * 12 modules) = 0.00462962963
    calibSignal->SetReqFractionAboveAmpCutVal(0.004);
    calibSignal->SetReqFractionAboveAmp(kTRUE);
    calibSignal->SetUseAverage(kTRUE);
    calibSignal->SetSecInAverage(900);
    
    // select physics and calibration events now (only calibration in future)
    if ( aliHeader->Get("Type") == AliRawEventHeaderBase::kPhysicsEvent || 
	 aliHeader->Get("Type") == AliRawEventHeaderBase::kCalibrationEvent  ) {
      
      nevents++;
      if(nevents%1000==0)cout<<"Event "<<nevents<<endl;
      //  Signal calibration
      calibSignal->ProcessEvent(in,aliHeader);
    }
  } // loop over all events in file
    /* cleanup the reading handles */
  delete in;
  delete rawReader;    
  
  //
  // write results/histograms to rootfile
  //
  
  printf ("%d physics/calibration events processed.\n",nevents);
 
  Char_t outname[256];
  sprintf(outname, "LED_%09d.root",runno);
  
  calibSignal->Save(outname);
  printf("Wrote %s.\n",outname);
  calibSignal->Analyze();
  printf("NEvents Processed %d Accepted %d\n",calibSignal->GetNEvents(),calibSignal->GetNAcceptedEvents());
  
  // see if we can delete our analysis helper also
  delete calibSignal;
  
}
