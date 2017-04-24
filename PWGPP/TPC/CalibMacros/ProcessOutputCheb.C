void PrintProcStatus(int status) {
  fflush(stdout);
  printf("extraction of TPC SP calibration status: %d\n",status);
  fflush(stdout);
}

Bool_t ProcessOutputCheb(TString filesToProcess, Int_t startRun, Int_t endRun, const char* ocdbStorage, 
			 Bool_t corr=kTRUE, Bool_t dist=kTRUE) {

  // macro that process a list of files (xml or txt) to produce then the
  // OCDB entry for the TPC SP Distortion calibration; inspired by
  // AliAnalysisAlien::MergeOutput for simplicity

  enum {kStatusFail=1, kStatusOK=0, kStatusLowStatUpdate=-1};
  if (!gSystem->AccessPathName("CreateCorrMapObj.C")) {
    gROOT->LoadMacro("CreateCorrMapObj.C");
    ::Info("ProcessOutputCheb","using local version of ProcessOutputCheb.C");
  }
  else {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/CreateCorrMapObj.C");
    ::Info("ProcessOutputCheb","using %s","$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/CreateCorrMapObj.C");
  }
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/CreateCorrMapObj.C");
  //gROOT->LoadMacro("CreateCorrMapObj.C");
    
  Bool_t isGrid = kTRUE;
  TObjArray *listoffiles = new TObjArray();

  int ntrminUser = -1;
  int ntrminUserStop = -1;
  float stopDefFrac = 0.7;
  TString ntrminUserS = gSystem->Getenv("distMinTracks");
  if (!ntrminUserS.IsNull() && (ntrminUser=ntrminUserS.Atoi())>0) {
    ::Info("ProcessOutput","User provided min tracks to validate object as good: %d",ntrminUser);
  } 
  TString ntrminUserStopS = gSystem->Getenv("distMinTracksStop");
  if (!ntrminUserStopS.IsNull() && (ntrminUserStop=ntrminUserStopS.Atoi())>0) {
    ::Info("ProcessOutput","User provided min tracks to validate object as unusable: %d",ntrminUserStop);
  }
  
  if (filesToProcess.Contains(".xml")) {
    // Merge files pointed by the xml 
    TGridCollection *coll = (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"%s\");", filesToProcess.Data()));
    if (!coll) {
      ::Error("ProcessOutput", "Input XML collection empty.");
      PrintProcStatus(kStatusFail);
      return kFALSE;
    }
    coll->Print();
    // Iterate grid collection
    while (coll->Next()) {
      TString fname = coll->GetTURL();
      Printf("fname = %s", fname.Data());
      listoffiles->Add(new TNamed(fname.Data(),""));
    }   
  }
  
  else if (filesToProcess.Contains(".txt")) {
    TString line;
    ifstream in;
    in.open(filesToProcess);
    if (in.fail()) {
      ::Error("ProcessOutput", "File %s cannot be opened. Processing stopped." ,filesToProcess.Data());
      PrintProcStatus(kStatusFail);
      return kFALSE;
    }
    Int_t nfiles = 0;
    while (in.good()) {
      in >> line;
      if (line.IsNull()) continue;
      if (!line.Contains("alien:")) isGrid = kFALSE;
      nfiles++;
      listoffiles->Add(new TNamed(line.Data(),""));
    }
    in.close();
    if (!nfiles) {
      ::Error("ProcessOutput","Input file %s contains no files to be processed\n", filesToProcess.Data());
      delete listoffiles;
      PrintProcStatus(kStatusFail);
      return kFALSE;
    }
  }
  
  if (!listoffiles->GetEntries()) {
    ::Error("ProcessOutput","No files to process\n");
    delete listoffiles;
    PrintProcStatus(kStatusFail);
    return kFALSE;
  }
  
  if (startRun != endRun) {
    Printf("The processing now is only run-level, please check again!");
    PrintProcStatus(kStatusFail);
    return kFALSE;
  }
  
  Int_t run = startRun;
  TObject *nextfile;
  TString snextfile;
  TIter next(listoffiles);   
  TObjArray* acorr = new TObjArray();
  TObjArray* adist = new TObjArray();
  acorr->SetOwner(kTRUE);
  adist->SetOwner(kTRUE);


  Int_t lowStatJobs = 0, badStatJobs = 0;
  Int_t nJobs = 0;
  while (nextfile=next()) {
    snextfile = nextfile->GetName();
    Printf("opening file %s", snextfile.Data());
    if (isGrid) TGrid::Connect("alien://");
    AliTPCDcalibRes* dcalibRes = AliTPCDcalibRes::Load(snextfile.Data());
    if (!dcalibRes) {
      ::Error("ProcessOutput","Did not find calib object in %s, job Killed",snextfile.Data());
      PrintProcStatus(kStatusFail);
      exit(1);
    }
    int ntrUse = dcalibRes->GetNTracksUsed();
    int ntrMin = dcalibRes->GetMinTrackToUse();
    if (ntrminUser>0) ntrMin = ntrminUser;
    int ntrMinBad = ntrminUserStop>0 ? ntrminUserStop : int(ntrMin*stopDefFrac);
    if (ntrMinBad>=ntrMin) ntrMinBad = ntrMin;
    if (ntrUse<ntrMin) {
      ::Warning("ProcessOutput","Low stat:%d tracks used (min: %d, unusable: %d) in %s",ntrUse,ntrMin,ntrMinBad,snextfile.Data());
      lowStatJobs++;
      if (ntrUse<ntrMinBad) badStatJobs++;
    }
    else {
      ::Info("ProcessOutput","stat is OK :%d tracks used (min: %d, unusable: %d) in %s",ntrUse,ntrMin,ntrMinBad,snextfile.Data());
    }
    if (corr) { 
      AliTPCChebCorr* c = dcalibRes->GetChebCorrObject();
      if (!c) {
	::Error("ProcessOutput","Did not find %s Cheb.parm in %s","Correction",snextfile.Data());
	PrintProcStatus(kStatusFail);
	exit(1);
      }
      acorr->Add(c);
    }
    if (dist) { // 
      AliTPCChebDist* d = dcalibRes->GetChebDistObject();
      if (!d) {
	::Error("ProcessOutput","Did not find %s Cheb.parm in %s","Distortion" ,snextfile.Data());
	//PrintProcStatus(kStatusFail);
	//exit(1);
      }
      else adist->Add(d);
    }
    //
    nJobs++;
  }
  if (badStatJobs) {
    ::Error("ProcessOutput","%d out of %d timebins have too low stat, will not update OCDB",badStatJobs,nJobs);
    PrintProcStatus(kStatusFail);
  }
  else {
    int status = kStatusOK;
    if (lowStatJobs) {
      ::Warning("ProcessOutput","%d out of %d timebins have low stat, still will update OCDB",lowStatJobs,nJobs);
      status = kStatusLowStatUpdate;
    }
    if (corr) {
      printf("Corrections\n");
      acorr->Print();
      if (!CreateCorrMapObjTime(acorr, startRun, endRun, ocdbStorage)) status = kStatusFail;      
    }
    if (dist && adist->GetEntries()) {
      printf("Distortions\n");
      adist->Print();
      if (!CreateCorrMapObjTime(adist, startRun, endRun, ocdbStorage)) status = kStatusFail;
    }
    PrintProcStatus(status);
  }
  acorr->SetOwner(kFALSE);
  adist->SetOwner(kFALSE);
  delete acorr;
  delete adist;
  delete listoffiles;
  
  return kTRUE;
  
}
