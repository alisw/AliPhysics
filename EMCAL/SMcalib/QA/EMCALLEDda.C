void EMCALLEDda(const int year=12, const int runno = 191741, const int streamno=8, const int segno=10, const int kMaxPhysEvents = 1000, const int debug = 1){
  
  int i, status;
  gROOT->SetStyle("Plain");
  
  /* log start of process */
  printf("EMCAL DA started - %s\n",__FILE__);
  
  AliCaloCalibSignal * calibSignal = new 
    AliCaloCalibSignal(AliCaloCalibSignal::kEmCal); 

  // setup
  calibSignal->SetAmpCut(1);
  calibSignal->SetReqFractionAboveAmpCutVal(0.1);
  calibSignal->SetReqFractionAboveAmp(false);

  // setup; LEDRef
  calibSignal->SetAmpCutLEDRef(1);
  calibSignal->SetReqLEDRefAboveAmpCutVal(false);

  int nevents=0;

  Char_t fname[256];
  sprintf(fname, "%02d%09d0%02d.%d.root", year, runno, streamno, segno);
  AliRawReader *rawReader = NULL;
  rawReader = new AliRawReaderRoot(fname);

  AliCaloRawStreamV3 *in = NULL; 
  in = new AliCaloRawStreamV3(rawReader,"EMCAL");
  rawReader->Select("EMCAL", 0, AliEMCALGeoParams::fgkLastAltroDDL) ; //select EMCAL DDL range
  //in->SetOldRCUFormat(kTRUE);
  
  int nev=0;
  int type = 0;
  UInt_t timestamp = 0;

  /* read until EOF */
  while ( rawReader->NextEvent() && nevents<kMaxPhysEvents) {
    
    calibSignal->SetRunNumber( rawReader->GetRunNumber() ); // just for fun; keep info on last run looked at
    //    cout << " Event " << nev << endl;

    type = rawReader->GetType();
    timestamp = rawReader->GetTimestamp();
    
    // select physics and calibration events now (only calibration in future)
    if ( type == AliRawEventHeaderBase::kCalibrationEvent  ) {

      rawReader->Reset();
      
      nevents++;
      // if(nevents%1000==0)cout<<"Event "<<nevents<<endl;
      calibSignal->ProcessEvent(in, timestamp);
    }

    nev++;
    
  } // loop over all events in file
    /* cleanup the reading handles */
  delete in;
  delete rawReader;    

  // calculate average values also, for the LED info
  calibSignal->SetUseAverage(kTRUE);
  calibSignal->Analyze();
  
  // by default, we only save the full info in debug mode  
  if (debug==0) {
    // reset the full trees, when we are not in debug mode
    calibSignal->GetTreeAmpVsTime()->Reset();
    calibSignal->GetTreeLEDAmpVsTime()->Reset();
  }

  //
  // write results/histograms to rootfile
  //
  
  printf ("%d physics/calibration events processed.\n",nevents);
  
  // create output histograms and write to file
  Char_t outname[256];
  sprintf(outname, "LED_%09d.root",runno);
  //sprintf(outname, "EMCALLED.root");

  TFile destFile(outname, "update");
  destFile.cd();
  calibSignal->Write("emcCalibSignal");
  destFile.Close();

  printf("Wrote %s.\n",outname);
  
  // see if we can delete our analysis helper also
  delete calibSignal;
  
}
