void PrintProcStatus(int status) {
  ::Info("ProcessOutput","processing ends with status %d",status);
}

Bool_t ProcessOutputCheb(TString filesToProcess, Int_t startRun, Int_t endRun, const char* ocdbStorage, Bool_t corr=kTRUE) {

  // macro that process a list of files (xml or txt) to produce then the
  // OCDB entry for the TPC SP Distortion calibration; inspired by
  // AliAnalysisAlien::MergeOutput for simplicity

  enum {kStatusFail=-1, kStatusOK=0, kStatusNoUpdate=1};

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/CreateCorrMapObj.C");
  //gROOT->LoadMacro("CreateCorrMapObj.C");
    
  Bool_t isGrid = kTRUE;
  TObjArray *listoffiles = new TObjArray();

  int ntrminUser = -1;
  TString ntrminUserS = gSystem->Getenv("distMinTracks");
  if (!ntrminUserS.IsNull() && (ntrminUser=ntrminUserS.Atoi())>0) {
    ::Info("ProcessOutput","User provided min tracks to validate object: %d",ntrminUser);
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
  TObjArray* a = new TObjArray();
  a->SetOwner(kTRUE);


  Int_t lowStatJobs = 0;
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
    if (ntrUse<ntrMin) {
      ::Error("ProcessOutput","Low stat:%d tracks used (min: %d) in %s",ntrUse,ntrMin,snextfile.Data());
      lowStatJobs++;
    }
    else {
      ::Info("ProcessOutput","stat is OK :%d tracks used (min: %d) in %s",ntrUse,ntrMin,snextfile.Data());
    }
    AliTPCChebCorr* c = corr ? dcalibRes->GetChebCorrObject() : dcalibRes->GetChebDistObject();
    if (!c) {
      ::Error("ProcessOutput","Did not find %s Cheb.parm in %s",corr ? "Correction":"Distortion" ,snextfile.Data());
      PrintProcStatus(kStatusFail);
      exit(1);
    }
    a->Add(c);
    nJobs++;
  }
  if (lowStatJobs) {
    ::Error("ProcessOutput","%d out of %d timebins have low stat, will not update OCDB",lowStatJobs,nJobs);
    PrintProcStatus(kStatusNoUpdate);
  }
  else {
    a->Print();
    CreateCorrMapObjTime(a, startRun, endRun, ocdbStorage);
    PrintProcStatus(kStatusOK);
  }
  delete a;
  delete listoffiles;
  
  return kTRUE;
  
}
